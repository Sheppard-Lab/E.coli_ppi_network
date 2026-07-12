import igraph as ig
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------
# INPUT FILES
# ----------------------------------------------------------------------
nodes = pd.read_csv("/Users/simrangambhir/Desktop/network_perspective/version_2/nodes.tsv", sep="\t")
edges = pd.read_csv("/Users/simrangambhir/Desktop/network_perspective/version_2/edges_typed.tsv", sep="\t")

# Node fills 
CORE_FILL = "#d9d9d9"   
ACC_FILL  = "#585858"   # accessory = light grey
NODE_EDGE = "#000000"   # black outline on all nodes

#edge color 
E_CORE = "#d9d9d9"   
E_ACC  = "#585858"   # accessory = darker grey
E_MIX  = "#e8740c"   

# ----------------------------------------------------------------------
# CLEAN: keep valid edges, drop unconnected (singleton) nodes
# ----------------------------------------------------------------------
valid = set(nodes["locus_tag"])
edges = edges[edges["protein1"].isin(valid) & edges["protein2"].isin(valid)].copy()
connected = set(edges["protein1"]) | set(edges["protein2"])
nodes = nodes[nodes["locus_tag"].isin(connected)].reset_index(drop=True)
print(f"{len(nodes)} connected nodes, {len(edges)} edges")

ids = nodes["locus_tag"].tolist()
idx = {g_: i for i, g_ in enumerate(ids)}
edges = edges[edges["protein1"].isin(idx) & edges["protein2"].isin(idx)].copy()

# ----------------------------------------------------------------------
# BUILD GRAPH
# ----------------------------------------------------------------------
g = ig.Graph()
g.add_vertices(len(ids))
g.vs["gclass"] = nodes["gene_class"].tolist()
g.vs["gname"]  = nodes["gene_name"].fillna("").tolist()
g.add_edges(list(zip(edges["protein1"].map(idx), edges["protein2"].map(idx))))
g.es["weight"] = edges["score"].tolist()
g.es["itype"]  = edges["interaction_type"].tolist()

# ----------------------------------------------------------------------
# LOUVAIN community detection (reproducible: seed igraph's own RNG)
# ----------------------------------------------------------------------
random.seed(42)
ig.set_random_number_generator(random)
part = g.community_multilevel(weights="weight")
membership = np.array(part.membership)
n_comm = len(set(membership))
print(f"Louvain: {n_comm} communities, modularity {part.modularity:.3f}")

# ----------------------------------------------------------------------
# GROUP-AWARE LAYOUT: amplify within-community edges so each community
# contracts into its own visible hub
# ----------------------------------------------------------------------
base = np.array(g.es["weight"], dtype=float)
base = base / base.max()
same = np.array([membership[e.source] == membership[e.target] for e in g.es])
boost = np.where(same, 1.0, 0.05)        
lay_w = base * boost

random.seed(42)
ig.set_random_number_generator(random)
coords = np.array(g.layout_fruchterman_reingold(niter=1000, weights=lay_w.tolist()).coords)


coords[:, 0] *= 1.3

# ----------------------------------------------------------------------
# NODE COLOURS (by pangenome class) and SIZES (by degree)
# ----------------------------------------------------------------------
is_core = np.array([c == "core" for c in g.vs["gclass"]])
node_colors = np.where(is_core, CORE_FILL, ACC_FILL)

deg = np.array(g.degree())
node_sizes = 12 + (deg / deg.max()) * 180   # min size 12, scales with degree

# ----------------------------------------------------------------------
# FUNCTIONAL LABELS anchored to MARKER GENES 

# ----------------------------------------------------------------------
marker_labels = {
    "rpsL": "Translation &\ntranscription",
    "eno":  "Glycolysis /\npentose phosphate",
    "pykF": "Central carbon /\npyruvate metabolism",
    "guaA": "Amino acid &\nnucleotide biosynthesis",
    "malE": "ABC transporters /\nimport",
    "bamA": "Cell envelope &\nouter membrane",
    "recA": "DNA replication,\nrepair & stress",
    "pheA": "Aromatic amino acid\nbiosynthesis",
    "iscU": "tRNA modification &\nFe-S / Mo cofactors",
    "ftsZ": "Cell division &\npeptidoglycan",
    "fliC": "Flagella &\nchemotaxis",
    "narG": "Anaerobic respiration\n(nitrate/nitrite)",
    "nuoF": "Respiration\n(NADH / cytochrome)",
    "paaF": "Fatty acid \u03b2-oxidation /\nCoA metabolism",
    "fimA": "Fimbriae / pili\n(adhesion)",
    "mazF": "Toxin\u2013antitoxin\nsystems",
    "gspD": "Type II secretion",
    "puuA": "Polyamine\nmetabolism",
    "katG": "Oxidative stress\nresponse",
    "intE": "Prophage /\nmobile elements",
}

# resolve each marker gene -> its community id at runtime
gname_arr = np.array(g.vs["gname"])
community_labels = {}          # community_id -> label text
for marker, text in marker_labels.items():
    hits = np.where(gname_arr == marker)[0]
    if len(hits) == 0:
        print(f"  WARNING: marker '{marker}' not found, skipping label '{text.replace(chr(10),' ')}'")
        continue
    comm_id = membership[hits[0]]
    if comm_id not in community_labels:
        community_labels[comm_id] = text
    else:
        print(f"  NOTE: marker '{marker}' shares community {comm_id} with an existing label")

# ----------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(18, 13))

# edges coloured by interaction type
for itype, col, lw, a, z in [
    ("core-core", E_CORE, 0.45, 0.95, 3),
    ("accessory-accessory", E_ACC, 0.45, 0.90, 4),
    ("mixed", E_MIX, 0.55, 0.90, 6),
]:
    xs, ys = [], []
    for e in g.es:
        if e["itype"] == itype:
            s, t = e.source, e.target
            xs += [coords[s, 0], coords[t, 0], None]
            ys += [coords[s, 1], coords[t, 1], None]
    ax.plot(xs, ys, color=col, lw=lw, alpha=a, zorder=z, solid_capstyle="round")

# nodes: fill = pangenome class, size = degree, black outline
ax.scatter(coords[:, 0], coords[:, 1], s=node_sizes, c=node_colors,
           edgecolors=NODE_EDGE, linewidths=0.3, zorder=6)

# functional labels at community centroids (no box border, faint white backing)
for comm, text in community_labels.items():
    mask = membership == comm
    if mask.sum() == 0:
        continue
    cx = coords[mask, 0].mean()
    cy = coords[mask, 1].mean()
    ax.text(cx, cy, text, fontsize=14,
            ha="center", va="center", zorder=10, color="#111111",
            linespacing=0.95,
            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                      ec="none", alpha=0.7))

lo = np.percentile(coords, 2, axis=0)
hi = np.percentile(coords, 98, axis=0)
pad = (hi - lo) * 0.18
ax.set_xlim(lo[0]-pad[0], hi[0]+pad[0])
ax.set_ylim(lo[1]-pad[1], hi[1]+pad[1])

# legend
legend = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=CORE_FILL,
           markeredgecolor=NODE_EDGE, markersize=11, label='Core gene'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=ACC_FILL,
           markeredgecolor=NODE_EDGE, markersize=11, label='Accessory gene'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor="#888888",
           markeredgecolor=NODE_EDGE, markersize=6, label='Node size = number of interactions'),
    Line2D([0], [0], color=E_CORE, lw=2, label='core\u2013core interaction'),
    Line2D([0], [0], color=E_ACC, lw=2, label='accessory\u2013accessory interaction'),
    Line2D([0], [0], color=E_MIX, lw=2, label='core\u2013accessory interaction'),
]
ax.legend(handles=legend, loc='upper left', frameon=False, fontsize=16)

ax.axis("off")
plt.tight_layout()
plt.savefig("/Users/simrangambhir/Desktop/network_perspective/version_2/louvain_hubs_labelled_cc_v7.svg",
            dpi=300, bbox_inches="tight")
print("saved louvain_hubs_labelled_cc.png")
