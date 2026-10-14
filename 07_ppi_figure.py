"""
07_ppi_figure.py

Redraws the STRING association network from the pipeline's own outputs
(results/string_ppi_edges.csv and results/hub_proteins.csv) with a layout chosen
for legibility. No new analysis: edges, degrees and hubs come straight from
06_ppi_network.py. Node positions carry no biological meaning.

Layout
  * hubs (degree >= 12) on an inner ring; the other proteins of the largest
    connected component on an outer ring, ordered by where their hub neighbours sit;
  * the remaining small connected components in a grid on the right.
Node size = number of connections. The three highest-degree hubs are labelled.
The layout is deterministic: the same inputs always give the same figure.

Run from the repository root:  python 07_ppi_figure.py
Output: results/ppi_network_clean.png
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

EDGES = "results/string_ppi_edges.csv"
HUBS = "results/hub_proteins.csv"
OUT = "results/ppi_network_clean.png"

SURFACE, INK, INK2, MUTED = "#ffffff", "#0b0b0b", "#52514e", "#898781"
ACCENT, NEUTRAL, EDGE = "#2a78d6", "#b9b8b0", "#c3c2b7"
HUB_MIN_DEGREE, N_LABELS = 12, 3
R_IN, R_OUT, CENTER = 0.42, 1.05, np.array([-0.45, 0.0])

e = pd.read_csv(EDGES)
G = nx.Graph()
G.add_edges_from(zip(e["preferredName_A"], e["preferredName_B"]))
hubs_df = pd.read_csv(HUBS).sort_values("Degree", ascending=False, kind="stable")
deg = dict(G.degree())
assert all(deg[g] == d for g, d in zip(hubs_df["Gene"], hubs_df["Degree"])), "degree mismatch"

# deterministic component order: largest first, ties broken alphabetically
comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
giant, small = comps[0], comps[1:]
hub_nodes = [n for n in hubs_df["Gene"] if deg[n] >= HUB_MIN_DEGREE]
assert all(h in giant for h in hub_nodes), "a hub is outside the largest component"
labels = list(hubs_df["Gene"][:N_LABELS])


def polar(angle_deg, r):
    a = math.radians(angle_deg)
    return CENTER + r * np.array([math.cos(a), math.sin(a)])


def circ_mean(angles):
    a = np.radians(angles)
    return math.degrees(math.atan2(np.sin(a).mean(), np.cos(a).mean())) % 360


# --- inner ring: the three labelled hubs at 12, 4 and 8 o'clock, other hubs in between
slots = [None] * len(hub_nodes)
for k, n in enumerate(hub_nodes[:N_LABELS]):
    slots[k * (len(hub_nodes) // N_LABELS)] = n
rest = iter(hub_nodes[N_LABELS:])
for i in range(len(slots)):
    if slots[i] is None:
        slots[i] = next(rest)
hub_angle = {n: (90 - 360 * i / len(slots)) % 360 for i, n in enumerate(slots)}

# --- outer ring: angle = circular mean of already-placed neighbours, then equal spacing
others = sorted(n for n in giant if n not in hub_angle)
raw = {}
for _ in range(6):
    for n in others:
        nb = [hub_angle[m] for m in sorted(G[n]) if m in hub_angle] or \
             [raw[m] for m in sorted(G[n]) if m in raw]
        if nb:
            raw[n] = circ_mean(nb)
for n in others:
    raw.setdefault(n, 270.0)
order = sorted(others, key=lambda n: (raw[n], n))
start_off = raw[order[0]]
out_angle = {n: (start_off + 360 * k / len(order)) % 360 for k, n in enumerate(order)}

pos = {n: polar(a, R_IN) for n, a in hub_angle.items()}
pos.update({n: polar(a, R_OUT) for n, a in out_angle.items()})

# --- small components: grid on the right, each on its own small circle
cols = 3
rows = math.ceil(len(small) / cols)
cell_w, cell_h = 0.58, 2.3 / rows
for i, comp in enumerate(small):
    r, c = divmod(i, cols)
    cx, cy = 1.0 + c * cell_w, 1.15 - (r + 0.5) * cell_h
    members = sorted(comp)
    scale = min(cell_w, cell_h) * 0.30
    for k, n in enumerate(members):
        if len(members) == 2:
            xy = np.array([-0.5 if k == 0 else 0.5, 0.0])
        else:
            a = 2 * math.pi * k / len(members) + math.pi / 2
            xy = np.array([math.cos(a), math.sin(a)])
        pos[n] = np.array([cx + xy[0] * scale, cy + xy[1] * scale])

is_hub = {n: deg[n] >= HUB_MIN_DEGREE for n in G}

fig = plt.figure(figsize=(8, 6), dpi=200, facecolor=SURFACE)
ax = fig.add_axes([0.02, 0.10, 0.96, 0.76])
ax.set_facecolor(SURFACE)
ax.axis("off")
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(-1.6, 2.61)
ax.set_ylim(-1.25, 1.25)

nx.draw_networkx_edges(G, pos, ax=ax, edge_color=EDGE, width=0.6, alpha=0.9)
for flag, colr in ((False, NEUTRAL), (True, ACCENT)):
    nodes = [n for n in G.nodes if is_hub[n] == flag]
    nx.draw_networkx_nodes(
        G, pos, nodelist=nodes, ax=ax, node_color=colr,
        node_size=[16 + 14 * deg[n] for n in nodes],
        edgecolors=SURFACE, linewidths=1.2,
    )

halo = [pe.withStroke(linewidth=3, foreground=SURFACE)]
for n in labels:
    a = hub_angle[n]
    x, y = polar(a, R_IN + 0.14)
    cos_a = math.cos(math.radians(a))
    ha = "left" if cos_a > 0.3 else "right" if cos_a < -0.3 else "center"
    ax.text(x, y, f"{n} ({deg[n]})", ha=ha, va="center", fontsize=8.6, fontweight="bold",
            color=INK, path_effects=halo, zorder=10)

fig.text(0.03, 0.945, "STRING associations among the 369 significant proteins",
         fontsize=13.5, fontweight="bold", color=INK, ha="left", va="center")
fig.text(0.03, 0.895,
         f"{G.number_of_nodes()} proteins have at least one high-confidence STRING association "
         f"(score ≥ 0.7; {G.number_of_edges()} edges). Node size = connections.",
         fontsize=8.2, color=INK2, ha="left", va="center")

lx, ly = 0.03, 0.082
for colr, txt, dx in ((ACCENT, f"Hubs (≥ {HUB_MIN_DEGREE} connections)", 0.0),
                      (NEUTRAL, "Other proteins", 0.27)):
    fig.patches.append(plt.Circle((lx + dx + 0.007, ly), 0.0075, transform=fig.transFigure,
                                  color=colr, figure=fig))
    fig.text(lx + dx + 0.022, ly, txt, fontsize=8.4, color=INK2, ha="left", va="center")
fig.text(0.03, 0.043, "Labelled: the three highest-degree hubs, with their number of connections in brackets.",
         fontsize=6.9, color=MUTED, ha="left", va="center")
fig.text(0.03, 0.017, "Data: PXD034525 (Merrihew et al. 2023), AD vs control hippocampus; STRING.  "
         "Code: github.com/AmirSedaghaati/proteomics-de-pipeline",
         fontsize=6.9, color=MUTED, ha="left", va="center")

fig.savefig(OUT, dpi=200, facecolor=SURFACE)
print(f"Saved {OUT}")
