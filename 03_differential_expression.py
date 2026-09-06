import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

print("Loading cleaned data...")
data = pd.read_csv("data/protein_abundance_clean.csv", index_col=0)
meta = pd.read_csv("data/sample_groups.csv", index_col=0)

ad_cols = meta.index[meta["Group"] == "AD"]
ctrl_cols = meta.index[meta["Group"] == "Control"]

print(f"AD samples: {len(ad_cols)} | Control samples: {len(ctrl_cols)}")

results = []
for protein, row in data.iterrows():
    ad_vals = row[ad_cols].values
    ctrl_vals = row[ctrl_cols].values

    # Welch's t-test: does not assume equal variance between groups, a safer
    # default than Student's t-test for real biological data.
    t_stat, p_val = stats.ttest_ind(ad_vals, ctrl_vals, equal_var=False)

    # log2FC directly from the log2-scale data: difference of means IS the log2
    # fold-change, since log2(A) - log2(B) = log2(A/B).
    log2fc = np.mean(ad_vals) - np.mean(ctrl_vals)

    results.append({"Gene": protein, "log2FC": log2fc, "p_value": p_val})

results_df = pd.DataFrame(results)

# Benjamini-Hochberg FDR correction across all tested proteins.
print("Applying Benjamini-Hochberg FDR correction...")
reject, q_values, _, _ = multipletests(results_df["p_value"], alpha=0.05, method="fdr_bh")
results_df["q_value"] = q_values
results_df["significant"] = reject

results_df = results_df.sort_values("p_value")
results_df.to_csv("data/differential_expression_results.csv", index=False)

n_sig = results_df["significant"].sum()
print(f"\n{n_sig} of {len(results_df)} proteins significant at FDR < 0.05")
print("Saved data/differential_expression_results.csv")
print("\nTop 10 hits:")
print(results_df.head(10).to_string(index=False))
