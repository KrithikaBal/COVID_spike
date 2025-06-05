#!/usr/bin/env python3
"""
Tree_by_Spike_GIDs_AboveBranches_HQ_Wide.py

1) Reads nextclade_output_Mar_May/nextclade.csv → qc=good rows.
2) Extracts Nextclade_pango lineage and Spike-profile (frozenset) for each sequence.
3) Counts “good” genomes per lineage; uses MIN_COUNT = 135 to identify “common” lineages.
4) Collapses “rare” lineages (<135) up to their nearest “common” ancestor.
5) Picks exactly one representative per collapsed lineage (mode spikeProfile).
6) Builds a grouped pseudo-tree: 
   - Internal node labels = G1, G2, … (one per unique spikeProfile).
   - Leaf labels = collapsed lineage name (e.g. “XEC”, “XEC.4”).
7) Plots the tree with:
   - G-ID labels **above** their branch.
   - Darker branches (linewidth=1.5).
   - Reduced vertical spacing (0.4″ per group).
   - Extra-wide canvas (40″ wide).
   - Larger colored dots (size=80) with thin black border.
8) Saves:
   - **Vector PDF** (`*_HQ_Wide.pdf`)
   - **SVG** (`*_HQ_Wide.svg`) — infinitely scalable.
   - **High‐dpi PNG** (600 dpi).
9) Also builds a legend mapping each G-ID → full spike‐mutation list, saved as PDF/SVG/PNG.
"""

import os
import sys
import pandas as pd
from collections import Counter, defaultdict

from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ──────────────────────────────────────────────────────────────────────────────
# 1) FORCE a color-capable backend + reset Matplotlib style
# ──────────────────────────────────────────────────────────────────────────────
matplotlib.use("Agg")
plt.style.use("default")
plt.rcParams["axes.prop_cycle"] = plt.rcParamsDefault["axes.prop_cycle"]
plt.rcParams["image.cmap"]     = "viridis"

# Increase default font size (tip labels, G-ID labels, legend)
plt.rcParams["font.size"] = 10

# ──────────────────────────────────────────────────────────────────────────────
# 2) File paths & sanity checks
# ──────────────────────────────────────────────────────────────────────────────
OUTDIR        = "nextclade_output_Mar_May"
NEXTCLADE_CSV = os.path.join(OUTDIR, "nextclade.csv")

# Tree outputs:
TREE_PDF = os.path.join(OUTDIR, "spike_grouped_tree_GIDs_Above_HQ_Wide.pdf")
TREE_SVG = os.path.join(OUTDIR, "spike_grouped_tree_GIDs_Above_HQ_Wide.svg")
TREE_PNG = os.path.join(OUTDIR, "spike_grouped_tree_GIDs_Above_HQ_Wide.png")

# Legend outputs:
LEG_PDF = os.path.join(OUTDIR, "spike_grouped_legend_GIDs_Above_HQ_Wide.pdf")
LEG_SVG = os.path.join(OUTDIR, "spike_grouped_legend_GIDs_Above_HQ_Wide.svg")
LEG_PNG = os.path.join(OUTDIR, "spike_grouped_legend_GIDs_Above_HQ_Wide.png")

if not os.path.isfile(NEXTCLADE_CSV):
    print(f"ERROR: {NEXTCLADE_CSV} not found. Run Nextclade first.", file=sys.stderr)
    sys.exit(1)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Read Nextclade CSV → filter qc=good → extract pango lineage & spikeProfile
# ──────────────────────────────────────────────────────────────────────────────
print("→ Reading Nextclade CSV and filtering qc=good …")
df = pd.read_csv(NEXTCLADE_CSV, sep=";")
df_good = df[df["qc.overallStatus"] == "good"].copy()
print(f"   • {len(df_good)} genomes passed QC.")

if "Nextclade_pango" not in df_good.columns:
    print("ERROR: Nextclade_pango column missing.", file=sys.stderr)
    sys.exit(1)

def extract_spike_set(aa_subs):
    """Return a frozenset of the ‘S:…’ entries from aaSubstitutions (or empty)."""
    if pd.isna(aa_subs) or aa_subs.strip() == "":
        return frozenset()
    return frozenset(x for x in aa_subs.split(",") if x.startswith("S:"))

df_good["spikeProfile"] = df_good["aaSubstitutions"].apply(extract_spike_set)
seq_to_lineage = dict(zip(df_good["seqName"], df_good["Nextclade_pango"]))
seq_to_spike   = dict(zip(df_good["seqName"], df_good["spikeProfile"]))

# ──────────────────────────────────────────────────────────────────────────────
# 4) Count “good” genomes per lineage & define MIN_COUNT = 135
# ──────────────────────────────────────────────────────────────────────────────
MIN_COUNT = 135
print(f"→ MIN_COUNT = {MIN_COUNT} …")

lineage_counts  = Counter(df_good["Nextclade_pango"])
common_lineages = {lin for lin, cnt in lineage_counts.items() if cnt >= MIN_COUNT}
rare_lineages   = set(lineage_counts.keys()) - common_lineages
print(f"   • {len(common_lineages)} lineages have ≥ {MIN_COUNT} genomes (common).")
print(f"   • {len(rare_lineages)} lineages have < {MIN_COUNT} (rare).")

# ──────────────────────────────────────────────────────────────────────────────
# 5) Collapse rare lineages → nearest common ancestor (strip “.<suffix>”)
# ──────────────────────────────────────────────────────────────────────────────
def find_common_ancestor(lin, common_set):
    """
    If ‘lin’ is in common_set, return it.
    Otherwise strip rightmost “.<suffix>” until in common_set or no dot remains.
    """
    if lin in common_set:
        return lin
    while "." in lin:
        lin = lin.rsplit(".", 1)[0]
        if lin in common_set:
            return lin
    return lin  # if nothing left

collapsed_for = {}
for seq, orig_lin in seq_to_lineage.items():
    if orig_lin in common_lineages:
        collapsed_for[seq] = orig_lin
    else:
        collapsed_for[seq] = find_common_ancestor(orig_lin, common_lineages)

collapsed_counts = Counter(collapsed_for.values())
print(f"   • After collapsing rare lineages, {len(collapsed_counts)} collapsed lineages remain.")

# ──────────────────────────────────────────────────────────────────────────────
# 6) Pick exactly ONE representative per collapsed lineage (mode of spikeProfile)
# ──────────────────────────────────────────────────────────────────────────────
print("→ Picking one representative (mode Spike-profile) per collapsed lineage …")

lineage_to_profiles = defaultdict(list)
lineage_to_seqs     = defaultdict(list)

for seq, col_lin in collapsed_for.items():
    lineage_to_profiles[col_lin].append(seq_to_spike[seq])
    lineage_to_seqs[col_lin].append(seq)

lin_to_mode_spike = {}
reps = set()

for lin, profiles in lineage_to_profiles.items():
    mode_profile = Counter(profiles).most_common(1)[0][0]
    lin_to_mode_spike[lin] = mode_profile
    # pick one seq in that lineage whose spikeProfile == mode_profile
    for seq in lineage_to_seqs[lin]:
        if seq_to_spike[seq] == mode_profile:
            reps.add(seq)
            break

print(f"   • {len(reps)} representatives chosen (one per collapsed lineage).")

# ──────────────────────────────────────────────────────────────────────────────
# 7) Build a grouped pseudo-tree with short SpikeGroup IDs (G1, G2, …)
# ──────────────────────────────────────────────────────────────────────────────
print("→ Building grouped-by-spike pseudo-tree with G-IDs …")

# 7a) Gather reps under each mode_spike_profile
spike_to_reps = defaultdict(list)
for seq in reps:
    col_lin = collapsed_for[seq]
    mode_profile = lin_to_mode_spike[col_lin]
    spike_to_reps[mode_profile].append(seq)

mode_profiles = sorted(spike_to_reps.keys(), key=lambda s: ",".join(sorted(s)))
n_profiles     = len(mode_profiles)

# 7b) Assign short IDs: G1, G2, …, G<n_profiles>
profile_to_gid = { prof: f"G{i+1}" for i, prof in enumerate(mode_profiles) }
gid_to_profile = { f"G{i+1}": prof for i, prof in enumerate(mode_profiles) }

# 7c) Generate a distinct color for each profile
cmap = plt.get_cmap("gist_ncar", n_profiles)
profile_to_color = { prof: cmap(i) for i, prof in enumerate(mode_profiles) }

# 7d) Build the Newick: root → (G1 → leaves), (G2 → leaves), … (Gk → leaves)
root = Clade(name="", branch_length=0.0)

for prof in mode_profiles:
    gid   = profile_to_gid[prof]      # e.g. "G1"
    color = profile_to_color[prof]

    # Internal node = Clade(name="G<ID>", branch_length=0)
    profile_node = Clade(name=gid, branch_length=0.0)

    # Under that node, add one leaf per representative, labeled by collapsed lineage
    for seq in sorted(spike_to_reps[prof]):
        collapsed_lin = collapsed_for[seq]  # e.g. "XEC" or "XEC.4"
        tip = Clade(name=collapsed_lin, branch_length=0.0)
        tip.tip_color = color
        profile_node.clades.append(tip)

    root.clades.append(profile_node)

grouped_tree = Tree(root, rooted=True)

# ──────────────────────────────────────────────────────────────────────────────
# 8) Plot the grouped tree with:
#    - G-ID labels above the branch
#    - darker, thicker branch lines
#    - reduced vertical spacing (0.4" per group)
#    - extra-wide canvas (40" wide)
#    - larger colored dots
# ──────────────────────────────────────────────────────────────────────────────
print("→ Plotting grouped-by-spike tree (G-IDs above branches; high-quality; wide) …")

# 8a) Very wide figure: 40″ × (~0.4" × n_profiles, min 8″)
fig_width  = 40
fig_height = max(8, 0.4 * n_profiles)
fig = plt.figure(figsize=(fig_width, fig_height))
ax  = fig.add_subplot(1, 1, 1)

# Increase default line width for branches:
plt.rcParams["lines.linewidth"] = 1.5

# Draw branches and labels (internal = G<ID>, leaf = collapsed lineage) in black
# show_confidence=False hides bootstrap numbers
Phylo.draw(
    grouped_tree,
    do_show=False,
    axes=ax,
    label_colors=lambda c: "k",
    show_confidence=False
)

# 8b) DARKEN all branch lines by iterating over Line2D objects:
for line in ax.get_lines():
    line.set_color("black")
    line.set_linewidth(1.5)

# 8c) Gather exact (x, y) for each terminal (leaf) if Biopython set them
exact_xy = {}
for clade in grouped_tree.get_terminals():
    if hasattr(clade, "x") and hasattr(clade, "y"):
        exact_xy[clade.name] = (clade.x, clade.y)

# 8d) Fallback shift if some terminal lacks (x, y)
x_min, x_max = ax.get_xlim()
delta_label  = (x_max - x_min) * 0.0025  # ~0.25% of horizontal span

# 8e) Overlay larger colored dots at each leaf coordinate
dot_size = 80
for clade in grouped_tree.get_terminals():
    leaf_label = clade.name  # e.g. "XEC" or "XEC.4"
    color_rgba = clade.tip_color
    if leaf_label in exact_xy:
        x_i, y_i = exact_xy[leaf_label]
        ax.scatter(
            [x_i], [y_i],
            s=dot_size,
            facecolor=[color_rgba],
            edgecolors="black",
            linewidths=0.3,
            transform=ax.transData
        )
    else:
        # fallback: shift left from label text
        for txt in ax.texts:
            if txt.get_text().strip() == leaf_label:
                x_text, y_text = txt.get_position()
                ax.scatter(
                    [x_text - delta_label], [y_text],
                    s=dot_size,
                    facecolor=[color_rgba],
                    edgecolors="black",
                    linewidths=0.3,
                    transform=ax.transData
                )
                break

# 8f) Move each G-ID label **above** its branch
#     We find text objects whose .get_text() == "G1", "G2", etc. and shift them upward.
#     Compute a small y-offset as fraction of y-axis range:
y_min, y_max = ax.get_ylim()
delta_y = (y_max - y_min) * 0.02  # ~2% of vertical span

for txt in ax.texts:
    txt_str = txt.get_text().strip()
    if txt_str in gid_to_profile:  # it's an internal node label "G1", "G2", ...
        x0, y0 = txt.get_position()
        txt.set_position((x0, y0 + delta_y))
        txt.set_fontweight("bold")
        txt.set_fontsize(11)

    else:
        # Leaf labels (collapsed lineage names): keep them horizontally just to the right of their dot
        txt.set_fontsize(10)
        txt.set_rotation(0)
        txt.set_ha("left")

ax.set_title(
    "Spike-Profile Grouped Tree (G-ID labels above branches)",
    fontsize=14,
    fontweight="bold"
)
ax.axis("off")

plt.tight_layout()

# 8g) Save as vector PDF, SVG, and high-dpi PNG
print(f"   • Saving vector PDF: {TREE_PDF}")
fig.savefig(TREE_PDF, bbox_inches="tight")

print(f"   • Saving SVG (vector): {TREE_SVG}")
fig.savefig(TREE_SVG, bbox_inches="tight")

print(f"   • Saving high-dpi PNG (600 dpi): {TREE_PNG}")
fig.savefig(TREE_PNG, bbox_inches="tight", dpi=600)

plt.close(fig)
print(f"→ Saved grouped tree:\n   {TREE_PDF}\n   {TREE_SVG}\n   {TREE_PNG}")

# ──────────────────────────────────────────────────────────────────────────────
# 9) Build & save legend (G-ID → comma-joined S:…) in PDF/SVG/PNG
# ──────────────────────────────────────────────────────────────────────────────
print("→ Building separate legend (G<ID> → spike mutations) …")

leg_fig = plt.figure(figsize=(10, max(3, 0.4 * n_profiles)))
ax2     = leg_fig.add_subplot(1, 1, 1)
ax2.axis("off")

handles, labels = [], []
for prof in mode_profiles:
    gid       = profile_to_gid[prof]
    label_str = ", ".join(sorted(prof)) if prof else "(no S: mutation)"
    color_rgba = profile_to_color[prof]

    handles.append(Line2D(
        [0], [0],
        marker="o",
        color="black",
        markerfacecolor=color_rgba,
        markersize=10,
        linewidth=0.3
    ))
    labels.append(f"{gid}: {label_str}")

legend = ax2.legend(
    handles=handles,
    labels=labels,
    title="SpikeGroup ID  →  Mode Spike-profile (S:…)",
    loc="upper left",
    bbox_to_anchor=(0, 1),
    fontsize=10,
    title_fontsize=11,
    ncol=1
)

for text in legend.get_texts():
    text.set_fontsize(10)
legend.get_title().set_fontsize(11)

plt.tight_layout()

print(f"   • Saving vector PDF: {LEG_PDF}")
leg_fig.savefig(LEG_PDF, bbox_inches="tight")

print(f"   • Saving SVG (vector): {LEG_SVG}")
leg_fig.savefig(LEG_SVG, bbox_inches="tight")

print(f"   • Saving high-dpi PNG (600 dpi): {LEG_PNG}")
leg_fig.savefig(LEG_PNG, bbox_inches="tight", dpi=600)

plt.close(leg_fig)
print(f"→ Saved legend:\n   {LEG_PDF}\n   {LEG_SVG}\n   {LEG_PNG}")

print("All done (G-IDs above branches; wide; high-quality).")
