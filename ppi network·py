import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

de = pd.read_csv("data/differential_expression_results.csv")
sig_genes = de.loc[de["significant"], "Gene"].dropna().unique().tolist()

print(f"Querying STRING API for {len(sig_genes)} significant proteins...")

STRING_API_URL = "https://string-db.org/api/tsv/network"
params = {
    "identifiers": "%0d".join(sig_genes),
    "species": 9606,           # Homo sapiens
    "required_score": 700,     # high-confidence interactions only
}
response = requests.get(STRING_API_URL, params=params, timeout=60)
response.raise_for_status()

lines = [l for l in response.text.strip().split("\n") if l]
header = lines[0].split("\t")
rows = [l.split("\t") for l in lines[1:]]
edges_df = pd.DataFrame(rows, columns=header)
edges_df.to_csv("results/string_ppi_edges.csv", index=False)
print(f"Retrieved {len(edges_df)} high-confidence interactions from STRING")

# --- Build network and identify hub proteins by degree centrality ---
print("Building network and ranking hub proteins by degree...")
G = nx.Graph()
for _, row in edges_df.iterrows():
    G.add_edge(row["preferredName_A"], row["preferredName_B"])

degree = dict(G.degree())
hub_ranking = pd.DataFrame(
    {"Gene": list(degree.keys()), "Degree": list(degree.values())}
).sort_values("Degree", ascending=False)
hub_ranking.to_csv("results/hub_proteins.csv", index=False)

print("\nTop 10 hub proteins (by PPI degree):")
print(hub_ranking.head(10).to_string(index=False))

# --- Visualize network, sized by degree ---
plt.figure(figsize=(10, 10))
pos = nx.spring_layout(G, seed=42, k=0.4)
sizes = [degree[n] * 60 for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color="#4c72b0", alpha=0.8)
nx.draw_networkx_edges(G, pos, alpha=0.3)
top_hubs = set(hub_ranking.head(10)["Gene"])
labels = {n: n for n in G.nodes() if n in top_hubs}
nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)
plt.title("Protein-Protein Interaction Network (significant DE proteins, hippocampus)")
plt.axis("off")
plt.tight_layout()
plt.savefig("results/ppi_network.png", dpi=150)
plt.close()

print("\nSaved: results/string_ppi_edges.csv, results/hub_proteins.csv, results/ppi_network.png")
print(
    "\nNote: this mirrors the hub-gene identification approach from the candidate's "
    "own Scientific African paper (CAMK1G, PRKCB, MAP3K9), now applied starting "
    "from real protein-level differential abundance instead of transcript-level data."
)
