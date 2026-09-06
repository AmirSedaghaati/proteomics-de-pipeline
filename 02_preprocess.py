import re
import numpy as np
import pandas as pd

print("Loading raw protein abundance table...")
raw = pd.read_csv("data/Hipp-protGrp.batchadj.csv", index_col=0)

print("Loading sample metadata...")
meta = pd.read_csv("data/Hipp-meta.csv")

# The source data is already log2-scale normalized protein abundance
# (Level 3B = "normalized protein abundance across all batches"), so no
# additional log2 transform is applied here - re-log-transforming already
# log-scale data would distort the distribution.

# --- Step 1: drop pooled QC / reference samples ---
# These are internal quality-control pools (Brain Reference, Hipp Reference),
# not individual patient samples, and must be excluded from a disease-vs-control
# comparison.
qc_samples = meta.loc[meta["Condition"].isin(["Brain Reference", "Hipp Reference"]), "Sample Label"]
real_samples = [c for c in raw.columns if c not in qc_samples.values]
print(f"Excluding {len(qc_samples)} pooled QC/reference samples: {list(qc_samples)}")
data = raw[real_samples].copy()

meta_real = meta[meta["Sample Label"].isin(real_samples)].set_index("Sample Label")

# --- Step 2: collapse the 4-way clinical grouping into a 2-way comparison ---
# Sporadic AD + AutoDom AD (autosomal-dominant AD) -> "AD"
# Control High Path + Control Low Path -> "Control"
# (Both control subtypes are cognitively normal; "High/Low Path" refers to
# histopathologic AD burden found at autopsy despite no dementia diagnosis.)
group_map = {
    "Sporadic AD": "AD",
    "AutoDom AD": "AD",
    "Control High Path": "Control",
    "Control Low Path": "Control",
}
meta_real["Group"] = meta_real["Condition"].map(group_map)
print(meta_real["Group"].value_counts())

# --- Step 3: extract a readable gene/protein symbol from the row descriptor ---
# Descriptor format ends with "... @ sp|ACCESSION|GENE_HUMAN, sp|ACCESSION2|GENE2_HUMAN"
# Take the first UniProt entry's gene symbol as the representative label for that
# protein group (a protein group can map to more than one near-identical isoform;
# this keeps one readable label per row without discarding the group).
def extract_gene(descriptor):
    match = re.search(r"sp\|[^|]+\|([A-Za-z0-9]+)", descriptor, re.IGNORECASE)
    return match.group(1) if match else descriptor[:30]

data.index = [extract_gene(d) for d in data.index]
data.index.name = "Gene"

# --- Step 4: treat 0 as missing (not-detected), not a true zero abundance ---
# This dataset encodes "peptide/protein not detected in this run" as 0, which is
# a detection-limit / missingness signal, not a real abundance of zero.
data = data.replace(0, np.nan)

# --- Step 5: filter out protein groups with too much missingness ---
# Require detection in at least 70% of samples in at least ONE of the two groups
# (a standard proteomics filtering rule - a protein consistently missing in one
# group but present in the other is itself potentially informative, so we don't
# require 70% in both groups, only in at least one).
ad_cols = meta_real.index[meta_real["Group"] == "AD"]
ctrl_cols = meta_real.index[meta_real["Group"] == "Control"]

frac_present_ad = data[ad_cols].notna().mean(axis=1)
frac_present_ctrl = data[ctrl_cols].notna().mean(axis=1)
keep = (frac_present_ad >= 0.7) | (frac_present_ctrl >= 0.7)
print(f"Keeping {keep.sum()} of {len(keep)} protein groups after missingness filter (>=70% detected in at least one group)")
data_filtered = data.loc[keep]

# --- Step 6: impute remaining missing values ---
# Proteomics missingness is typically MNAR (missing not at random - low-abundance
# proteins are more likely to go undetected, not random dropout). A simple,
# honest left-censored imputation: replace a protein's missing values with a value
# drawn from a normal distribution centered below that protein's observed minimum,
# with a small spread. This is a defensible first-pass method - not as rigorous as
# dedicated tools (e.g., MSstats, DEP's QRILC/MinProb), which is stated plainly in
# the README's Limitations section rather than overclaiming statistical rigor.
rng = np.random.default_rng(42)

def impute_row(row):
    observed = row.dropna()
    if observed.empty or row.isna().sum() == 0:
        return row
    shift = 1.8       # shift downward, in log2 units, below the observed minimum
    scale = 0.3        # narrow spread, typical for this kind of imputation
    center = observed.min() - shift
    n_missing = row.isna().sum()
    imputed_values = rng.normal(loc=center, scale=scale, size=n_missing)
    row_filled = row.copy()
    row_filled[row_filled.isna()] = imputed_values
    return row_filled

print("Imputing remaining missing values (left-censored, per-protein)...")
data_imputed = data_filtered.apply(impute_row, axis=1)

data_imputed.to_csv("data/protein_abundance_clean.csv")
meta_real.to_csv("data/sample_groups.csv")
print("Saved data/protein_abundance_clean.csv and data/sample_groups.csv")
