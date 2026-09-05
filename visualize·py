import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs("results", exist_ok=True)

de = pd.read_csv("data/differential_expression_results.csv")
data = pd.read_csv("data/protein_abundance_clean.csv", index_col=0)
meta = pd.read_csv("data/sample_groups.csv", index_col=0)
raw = pd.read_csv("data/Hipp-protGrp.batchadj.csv", index_col=0)

# --- Volcano plot ---
print("Building volcano plot...")
de["neg_log10_p"] = -np.log10(de["p_value"])
colors = np.where(
    de["significant"], np.where(de["log2FC"] > 0, "#d62728", "#1f77b4"), "#999999"
)

plt.figure(figsize=(8, 6))
plt.scatter(de["log2FC"], de["neg_log10_p"], c=colors, s=12, alpha=0.7, edgecolors="none")
plt.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1)
plt.axvline(0, color="gray", linestyle="-", linewidth=0.5)
plt.xlabel("log2 Fold Change (AD vs Control)")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot: Hippocampal Protein Abundance, AD vs Control")
plt.tight_layout()
plt.savefig("results/volcano_plot.png", dpi=150)
plt.close()

# --- Heatmap of top differentially abundant proteins ---
print("Building heatmap of top hits...")
top_genes = de.sort_values("p_value").head(30)["Gene"].tolist()
heatmap_data = data.loc[data.index.intersection(top_genes)]

# z-score each protein across samples so the heatmap shows relative pattern,
# not absolute abundance scale (which differs protein to protein).
heatmap_z = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1), axis=0)

# order columns by group for a visually readable heatmap
ordered_cols = list(meta.index[meta["Group"] == "Control"]) + list(meta.index[meta["Group"] == "AD"])
heatmap_z = heatmap_z[ordered_cols]

plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_z, cmap="vlag", center=0, xticklabels=True, yticklabels=True,
            cbar_kws={"label": "z-score"})
plt.title("Top 30 Differentially Abundant Proteins (Hippocampus, AD vs Control)")
plt.xlabel("Sample")
plt.ylabel("Gene / Protein")
plt.tight_layout()
plt.savefig("results/heatmap_top_proteins.png", dpi=150)
plt.close()

# --- Missing-value / detection-pattern plot (proteomics-specific) ---
# This has no transcriptomics equivalent - RNA-seq/microarray platforms don't
# have the same detection-limit-driven missingness that DIA mass spec does.
print("Building missing-value detection pattern plot...")
qc_samples = meta.index  # already excludes reference/QC pools from step 02
missing_frac_per_sample = (raw[qc_samples] == 0).mean(axis=0).sort_values()

plt.figure(figsize=(10, 5))
plt.bar(range(len(missing_frac_per_sample)), missing_frac_per_sample.values, color="#4c72b0")
plt.xticks(range(len(missing_frac_per_sample)), missing_frac_per_sample.index, rotation=90, fontsize=6)
plt.ylabel("Fraction of proteins not detected")
plt.title("Per-Sample Missingness (Detection-Limit Pattern), Hippocampus DIA-MS")
plt.tight_layout()
plt.savefig("results/missingness_pattern.png", dpi=150)
plt.close()

print("Saved: results/volcano_plot.png, results/heatmap_top_proteins.png, results/missingness_pattern.png")
