#!/usr/bin/env python3
"""
Tree.py

1) Read nextclade_output_april_may/nextclade.csv → keep only qc=good rows.
2) Extract Nextclade_pango lineage and Spike‐profile (frozenset) for each sequence.
3) Count “good” genomes per lineage; define MIN_COUNT = 50 as “common.”
4) Collapse any “rare” lineage (<50) by stripping trailing “.<suffix>” until you
   reach a “common” one.
5) For each collapsed lineage, pick exactly one representative (the genome whose
   Spike‐profile is most frequent—i.e., the mode—in that lineage).
6) Prune nextclade_output_april_may/nextclade.nwk to keep only those one‐per‐lineage
   representatives.
7) Draw the pruned tree (branches + black tip‐labels) so Biopython assigns clade.x/clade.y.
8) For each tip, place a colored dot exactly at the node (clade.x, clade.y).
   If clade.x/clade.y aren’t available, fall back to a small leftward shift from the tip label.
9) Save the resulting figure as PDF + 300 dpi PNG, plus a separate legend mapping
   each Spike‐profile → dot color.
"""

import os
import sys
import pandas as pd
from collections import Counter, defaultdict
from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ──────────────────────────────────────────────────────────────────────────────
# 1) Force a color‐capable backend + reset Matplotlib to default style
# ──────────────────────────────────────────────────────────────────────────────
matplotlib.use("Agg")
plt = matplotlib.pyplot
plt.style.use("default")
plt.rcParams["axes.prop_cycle"] = plt.rcParamsDefault["axes.prop_cycle"]
plt.rcParams["image.cmap"]     = "viridis"

# ──────────────────────────────────────────────────────────────────────────────
# 2) File paths & sanity checks
# ──────────────────────────────────────────────────────────────────────────────
OUTDIR        = "nextclade_output_april_may"
NEXTCLADE_CSV = os.path.join(OUTDIR, "nextclade.csv")
FULL_TREE     = os.path.join(OUTDIR, "nextclade.nwk")

TREE_PDF = os.path.join(OUTDIR, "lineage_tree_dots_fallback.pdf")
TREE_PNG = os.path.join(OUTDIR, "lineage_tree_dots_fallback.png")
LEG_PDF  = os.path.join(OUTDIR, "spike_legend.pdf")
LEG_PNG  = os.path.join(OUTDIR, "spike_legend.png")

for fp in (NEXTCLADE_CSV, FULL_TREE):
    if not os.path.isfile(fp):
        print(f"ERROR: {fp} not found. Run Nextclade first.", file=sys.stderr)
        sys.exit(1)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Read Nextclade CSV → filter qc=good → extract Nextclade_pango & spikeProfile
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
# 4) Count “good” genomes per lineage & define MIN_COUNT = 50
# ──────────────────────────────────────────────────────────────────────────────
MIN_COUNT = 25
print(f"→ MIN_COUNT = {MIN_COUNT} …")

lineage_counts = Counter(df_good["Nextclade_pango"])
common_lineages = {lin for lin, cnt in lineage_counts.items() if cnt >= MIN_COUNT}
print(f"   • {len(common_lineages)} lineages have ≥ {MIN_COUNT} genomes (common).")
print(f"   • {len(lineage_counts) - len(common_lineages)} lineages have < {MIN_COUNT} (rare).")

# ──────────────────────────────────────────────────────────────────────────────
# 5) Collapse rare lineages → nearest common ancestor (strip “.<suffix>”)
# ──────────────────────────────────────────────────────────────────────────────
def find_common_ancestor(lin, common_set):
    """
    If ‘lin’ is in common_set, return it.
    Otherwise strip off the rightmost “.<suffix>” until you land in common_set
    or no dot remains.
    """
    if lin in common_set:
        return lin
    while "." in lin:
        lin = lin.rsplit(".", 1)[0]
        if lin in common_set:
            return lin
    return lin

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
print("→ Picking one representative (mode Spike‐profile) per collapsed lineage …")

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
    # pick any one seq in that lineage with mode_profile
    for seq in lineage_to_seqs[lin]:
        if seq_to_spike[seq] == mode_profile:
            reps.add(seq)
            break

print(f"   • {len(reps)} representatives chosen (one per collapsed lineage).")

# ──────────────────────────────────────────────────────────────────────────────
# 7) Prune the full Nextclade tree to keep only these reps
# ──────────────────────────────────────────────────────────────────────────────
print("→ Loading & pruning full Nextclade tree …")
tree = Phylo.read(FULL_TREE, "newick")
all_tips = [c.name for c in tree.get_terminals()]

pruned = 0
for tip in all_tips:
    if tip not in reps:
        try:
            tree.prune(tip)
            pruned += 1
        except Exception:
            pass

print(f"   • Pruned {pruned} tips; now {len(tree.get_terminals())} tips remain.")

# ──────────────────────────────────────────────────────────────────────────────
# 8) Build a color palette for each lineage’s mode spikeProfile
# ──────────────────────────────────────────────────────────────────────────────
mode_profiles = set(lin_to_mode_spike.values())
n_profiles    = len(mode_profiles)
print(f"→ {n_profiles} unique mode Spike‐profiles across {len(reps)} lineages.")

if n_profiles <= 20:
    cmap = plt.cm.tab20
else:
    cmap = plt.cm.get_cmap("tab20b", n_profiles)

profile_to_color = {
    prof: cmap(i / (n_profiles - 1 if n_profiles > 1 else 1))
    for i, prof in enumerate(sorted(mode_profiles, key=lambda s: ",".join(sorted(s))))
}

# ──────────────────────────────────────────────────────────────────────────────
# 9) Annotate each tip (one per lineage) with collapsed lineage label & RGBA color
# ──────────────────────────────────────────────────────────────────────────────
print("→ Annotating tips with collapsed lineage label and RGBA color …")
for clade in tree.get_terminals():
    # Force the tip’s “name” to be just the collapsed lineage (so Biopython draws exactly that)
    collapsed_lin    = collapsed_for[clade.name]
    clade.name       = collapsed_lin
    clade.label_name = collapsed_lin
    clade.tip_color  = profile_to_color.get(
        lin_to_mode_spike[collapsed_lin], (0.6, 0.6, 0.6, 1.0)
    )

# ──────────────────────────────────────────────────────────────────────────────
# 10) Draw pruned tree (branches + black tip‐labels) and overlay smaller colored dots
# ──────────────────────────────────────────────────────────────────────────────
print("→ Plotting tree with colored‐dot fallback … (dots placed exactly at nodes)")

n_tips     = len(tree.get_terminals())
fig_height = max(6, 0.01 * n_tips)   # ~0.01" per tip

# Make figure wider to accommodate long branches:
fig = plt.figure(figsize=(15, fig_height))
ax  = fig.add_subplot(1, 1, 1)
plt.rcParams["font.size"] = 5    # slightly smaller font for many tips

# 10a) Draw branches + tip labels in black.  This will populate clade.x/clade.y if Biopython supports it
Phylo.draw(
    tree,
    label_func=lambda c: c.label_name if c.is_terminal() else "",
    label_colors=lambda c: "k",
    do_show=False,
    axes=ax
)

# ──────────────────────────────────────────────────────────────────────────────
# 10b) Gather exact (x,y) for any clade that got coordinates
# ──────────────────────────────────────────────────────────────────────────────
exact_xy = {}
for clade in tree.get_terminals():
    if hasattr(clade, "x") and hasattr(clade, "y"):
        exact_xy[clade.label_name] = (clade.x, clade.y)

# ──────────────────────────────────────────────────────────────────────────────
# 10c) Compute delta_label = small left‐shift for fallback (in case exact coords missing)
# ──────────────────────────────────────────────────────────────────────────────
x_min, x_max = ax.get_xlim()
delta_label = (x_max - x_min) * 0.0017  # ~1% of horizontal span

# ──────────────────────────────────────────────────────────────────────────────
# 10d) Overlay a smaller circle at each tip
# ──────────────────────────────────────────────────────────────────────────────
for clade in tree.get_terminals():
    lbl  = clade.label_name
    rgba = clade.tip_color

    # If Biopython provided exact (clade.x, clade.y), place dot exactly there
    if lbl in exact_xy:
        x_i, y_i = exact_xy[lbl]
        ax.scatter(
            [x_i], [y_i],
            s=20,                # smaller circle
            facecolor=[rgba],
            edgecolors="none",
            transform=ax.transData
        )
        continue

    # Otherwise, fallback: find the text position, then shift left by delta_label
    for txt in ax.texts:
        drawn = txt.get_text().strip()
        if drawn == lbl:
            x_text, y_text = txt.get_position()
            ax.scatter(
                [x_text - delta_label], [y_text],
                s=20,                # smaller circle
                facecolor=[rgba],
                edgecolors="none",
                transform=ax.transData
            )
            break

ax.set_title(
    f"One‐Tip‐Per‐Lineage Tree (MIN_COUNT={MIN_COUNT})\n"
    f"{n_tips} tips, each with a smaller colored dot at its node"
)
ax.axis("off")

plt.tight_layout()
fig.savefig(TREE_PDF, bbox_inches="tight", dpi=300)
fig.savefig(TREE_PNG, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"→ Saved:\n   {TREE_PDF}\n   {TREE_PNG}")

# ──────────────────────────────────────────────────────────────────────────────
# 11) Build & save a separate legend (Spike‐profile → dot color)
# ──────────────────────────────────────────────────────────────────────────────
print("→ Building separate Spike‐profile legend …")
profile_labels = [
    (", ".join(sorted(prof)) if prof else "(no S: mutation)", prof)
    for prof in sorted(mode_profiles, key=lambda s: ",".join(sorted(s)))
]

legend_height = max(2, 0.2 * len(profile_labels))
fig2 = plt.figure(figsize=(6, legend_height))
ax2  = fig2.add_subplot(1, 1, 1)
ax2.axis("off")

handles, labels = [], []
for label_str, prof in profile_labels:
    rgba = profile_to_color.get(prof, (0.6, 0.6, 0.6, 1.0))
    handles.append(Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor=rgba,
        markersize=6
    ))
    labels.append(label_str)

ax2.legend(
    handles,
    labels,
    title="Mode Spike‐profile (S:…)",
    loc="upper left",
    bbox_to_anchor=(0, 1),
    fontsize="small",
    title_fontsize="small",
    ncol=1
)

plt.tight_layout()
fig2.savefig(LEG_PDF, bbox_inches="tight")
fig2.savefig(LEG_PNG, dpi=300, bbox_inches="tight")
plt.close(fig2)
print(f"→ Saved legend:\n   {LEG_PDF}\n   {LEG_PNG}\nAll done.")
