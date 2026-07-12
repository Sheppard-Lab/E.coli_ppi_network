import igraph as ig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# INPUT FILES  (update these paths to your local folder)
# ----------------------------------------------------------------------
nodes = pd.read_csv("/path/nodes.tsv", sep="\t")
edges = pd.read_csv("/path/edges_typed.tsv", sep="\t")

CORE_C = "#d9d9d9"   # core = lighter grey
ACC_C  = "#585858"   # accessory = darker grey
MIX_C  = "#e8740c"   
LINE_C = "#000000"   

FS_LABEL = 15   # axis label font size 
FS_TICK  = 13   # tick font size
FS_VAL   = 13   # value labels on bars

# ----------------------------------------------------------------------
# GRAPH
# ----------------------------------------------------------------------
valid = set(nodes["locus_tag"])
edges = edges[edges["protein1"].isin(valid) & edges["protein2"].isin(valid)].copy()
connected = set(edges["protein1"]) | set(edges["protein2"])
nodes = nodes[nodes["locus_tag"].isin(connected)].reset_index(drop=True)
ids = nodes["locus_tag"].tolist()
idx = {g_: i for i, g_ in enumerate(ids)}
edges = edges[edges["protein1"].isin(idx) & edges["protein2"].isin(idx)].copy()

g = ig.Graph(); g.add_vertices(len(ids))
g.vs["gclass"] = nodes["gene_class"].tolist()
g.add_edges(list(zip(edges["protein1"].map(idx), edges["protein2"].map(idx))))
gclass = np.array(g.vs["gclass"]); is_core = gclass == "core"
deg = np.array(g.degree())
core_deg = deg[is_core]; acc_deg = deg[~is_core]

# interaction-type counts (panel C)
counts = edges["interaction_type"].value_counts()
cc = int(counts.get("core-core", 0))
aa = int(counts.get("accessory-accessory", 0))
mx = int(counts.get("mixed", 0))

# accessory-to-core connection (panel D)
acc = np.where(~is_core)[0]
has_core = np.array([any(gclass[n] == "core" for n in g.neighbors(v)) for v in acc])
n_conn = int(has_core.sum()); n_iso = len(acc) - n_conn

# ----------------------------------------------------------------------
# FIGURE: 3 stacked panels
# ----------------------------------------------------------------------
fig, (axB, axC, axD) = plt.subplots(
    3, 1, figsize=(4.2, 11.4), constrained_layout=True)

# ================= PANEL B: degree by class (log) =================
data = [np.log10(core_deg + 1), np.log10(acc_deg + 1)]
parts = axB.violinplot(data, showmedians=True, showextrema=True,
                       widths=0.7, bw_method=0.3)
for pc, col in zip(parts['bodies'], [CORE_C, ACC_C]):
    pc.set_facecolor(col); pc.set_alpha(1.0); pc.set_linewidth(0)
for key in ['cbars', 'cmins', 'cmaxes', 'cmedians']:
    parts[key].set_color(LINE_C); parts[key].set_linewidth(1.5)
axB.set_axisbelow(True); axB.yaxis.grid(True, color="0.8", linewidth=0.8)
axB.set_xticks([1, 2]); axB.set_xticklabels(["Core", "Accessory"], fontsize=FS_TICK)
yt = [1, 3, 10, 30, 100, 300]
axB.set_yticks([np.log10(v + 1) for v in yt]); axB.set_yticklabels(yt)
axB.tick_params(axis='y', labelsize=FS_TICK)
axB.set_ylabel("Number of interactions\nper gene", fontsize=FS_LABEL)
#axB.yaxis.set_label_coords(-0.18, 0.5)

# ================= PANEL C: interactions by type =================
labelsC = ["Core\u2013\ncore", "Accessory\u2013\naccessory", "Core\u2013\naccessory"]
valsC = [cc, aa, mx]; colsC = [CORE_C, ACC_C, MIX_C]
axC.set_axisbelow(True); axC.yaxis.grid(True, color="0.8", linewidth=0.8)
barsC = axC.bar(labelsC, valsC, color=colsC, linewidth=0, width=0.65)
for b, v in zip(barsC, valsC):
    axC.text(b.get_x()+b.get_width()/2, v + max(valsC)*0.015, f"{v:,}",
             ha="center", va="bottom", fontsize=FS_VAL)
axC.set_ylabel("Total number of\ninteractions", fontsize=FS_LABEL)
#axB.yaxis.set_label_coords(-0.18, 0.5)
axC.set_ylim(0, max(valsC)*1.12)
axC.tick_params(axis='x', labelsize=FS_TICK)
axC.tick_params(axis='y', labelsize=FS_TICK)


valsD = [n_conn, n_iso]
labelsD = ["Interacting\nwith core\ngenes", "Interacting with\nother\n accessory\ngenes"]
axD.set_axisbelow(True); axD.yaxis.grid(True, color="0.8", linewidth=0.8)
barsD = axD.bar([0, 1], valsD, color=[MIX_C, ACC_C], linewidth=0, width=0.6)
for b, v in zip(barsD, valsD):
    axD.text(b.get_x()+b.get_width()/2, v + max(valsD)*0.015, f"{v}",
             ha="center", va="bottom", fontsize=FS_VAL)
axD.set_xticks([0, 1]); axD.set_xticklabels(labelsD, fontsize=FS_TICK)
axD.set_ylabel("Number of\naccessory genes", fontsize=FS_LABEL)
#axB.yaxis.set_label_coords(-0.16, 0.5)
axD.set_ylim(0, max(valsD)*1.15)
axD.tick_params(axis='y', labelsize=FS_TICK)

for a in (axB, axC, axD):
    a.set_box_aspect(1)

plt.savefig("/path/combined_BCD_v8.svg", dpi=300, bbox_inches="tight")
print("saved combined_BCD.png")
