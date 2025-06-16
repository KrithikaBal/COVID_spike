#!/usr/bin/env python3
"""
Tree_by_Spike.py
=================================
Optimized version ensuring each lineage has a visible tip with colored labels,
shorter branches, and polished layout.

1) Load qc=good Nextclade CSV → extract Spike mutations.
2) Collapse rare lineages (<135) → common ancestor.
3) One representative per collapsed lineage (mode profile).
4) Encode profiles in binary matrix; cluster (Jaccard + UPGMA).
5) Scale branch lengths (×0.1) for compact layout.
6) Convert to Newick; read into Biopython.
7) Draw polished tree; color tip labels by Spike-profile.
8) Save PDF/SVG/PNG (600 dpi), Newick, and legend.

Requires: pandas, numpy, scipy, BioPython, matplotlib
"""
import os, sys
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Settings
matplotlib.use('Agg')
plt.rcParams.update({
    'font.size': 20,
    'lines.linewidth': 2.5,
    'image.cmap': 'viridis'
})

# Paths
out_dir = 'nextclade_output_Mar_June'
csv_in  = os.path.join(out_dir, 'nextclade_combined.csv')
tree_nwk= os.path.join(out_dir, 'spike_clustered_tree.nwk')
tree_pdf= os.path.join(out_dir, 'spike_clustered_tree.pdf')
tree_svg= os.path.join(out_dir, 'spike_clustered_tree.svg')
tree_png= os.path.join(out_dir, 'spike_clustered_tree.png')
leg_pdf = os.path.join(out_dir, 'spike_clustered_legend.pdf')
leg_svg = os.path.join(out_dir, 'spike_clustered_legend.svg')
leg_png = os.path.join(out_dir, 'spike_clustered_legend.png')
if not os.path.exists(csv_in):
    sys.exit(f"Missing {csv_in}")

# 1) Load and filter
print("→ Loading Nextclade data…")
df = pd.read_csv(csv_in, sep=';')
df = df[df['qc.overallStatus']=='good']
print(f"   • {len(df):,} QC-good genomes")

# Extract Spike profiles
def extract_spike(x):
    return frozenset(s for s in str(x).split(',') if s.startswith('S:'))

df['profile'] = df['aaSubstitutions'].apply(extract_spike)
seq2lin = dict(zip(df.seqName, df.Nextclade_pango))
seq2prof= dict(zip(df.seqName, df.profile))

# 2) Collapse rare lineages
MIN_COUNT = 100
counts = Counter(df.Nextclade_pango)
common = {lin for lin,c in counts.items() if c >= MIN_COUNT}
def collapse_lin(lin):
    while lin not in common and '.' in lin:
        lin = lin.rsplit('.',1)[0]
    return lin
collapsed = {s:(lin if lin in common else collapse_lin(lin)) for s,lin in seq2lin.items()}
print(f"   • {len(set(collapsed.values()))} collapsed lineages")

# 3) Representatives
by_lin = defaultdict(list)
for s,lin in collapsed.items():
    by_lin[lin].append(seq2prof[s])
reps, lineages = [], []
for lin,profs in by_lin.items():
    mode = Counter(profs).most_common(1)[0][0]
    for s in seq2lin:
        if collapsed[s]==lin and seq2prof[s]==mode:
            reps.append(s)
            lineages.append(lin)
            break
print(f"   • {len(reps)} representatives selected")

# Color mapping
profiles = [seq2prof[s] for s in reps]
unique_profiles = sorted(set(profiles), key=lambda p: ','.join(sorted(p)))
# Color mapping: use a qualitative palette visible on white background
unique_profiles = sorted(set(profiles), key=lambda p: ','.join(sorted(p)))
# switch to a distinct qualitative palette suitable for white background
cmap = plt.get_cmap('tab20', len(unique_profiles))
prof2col = {p: cmap(i) for i, p in enumerate(unique_profiles)}
lin2col  = {lin: prof2col[seq2prof[s]] for lin, s in zip(lineages, reps)}
prof2col = {p:cmap(i) for i,p in enumerate(unique_profiles)}
lin2col  = {lin:prof2col[seq2prof[s]] for lin,s in zip(lineages,reps)}

# 4) Binary matrix
all_muts = sorted({m for p in profiles for m in p})
bin_mat  = np.array([[1 if m in p else 0 for m in all_muts] for p in profiles])

# 5) Clustering + scaling branch lengths
print("→ Clustering Spike profiles…")
dists = pdist(bin_mat, lambda u,v:1-np.logical_and(u,v).sum()/np.logical_or(u,v).sum())
link = sch.linkage(dists, method='average')
link[:,2] *= 0.1  # shorten branches for compact layout

# 6) Convert to Newick
print("→ Converting to Newick…")
root = sch.to_tree(link)
def to_newick(node):
    if node.is_leaf(): return lineages[node.id]
    left = to_newick(node.left)
    right= to_newick(node.right)
    return f"({left}:{node.dist:.2f},{right}:{node.dist:.2f})"
nw = to_newick(root) + ';'
with open(tree_nwk,'w') as fh:
    fh.write(nw)
print(f"   • Newick saved: {tree_nwk}")

# 7) Draw polished tree with colored tip labels
print("→ Drawing polished tree with colored labels…")
fig_height = max(12, 0.3 * len(lineages))  # ~0.3" per lineage for better spacing
fig, ax = plt.subplots(figsize=(35, fig_height))
# ensure white background for figure and axes
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# Draw branches and colored labels
Phylo.draw(
    Phylo.read(tree_nwk, 'newick'),
    label_func=lambda c: c.name if c.is_terminal() else "",
    label_colors=lambda label: lin2col.get(label, 'black'),
    do_show=False,
    axes=ax,
    show_confidence=False
)
# Make tip labels bold
for txt in ax.texts:
    txt.set_fontweight('bold')


# Clean up and save
ax.axis('off')
plt.tight_layout()
for path in (tree_pdf, tree_svg):
    fig.savefig(path, bbox_inches='tight')
fig.savefig(tree_png, dpi=600, bbox_inches='tight')
print(f"   • Polished tree saved: {tree_pdf}, {tree_svg}, {tree_png}")

# 8) Legend
print("→ Building legend…")
fig2, ax2 = plt.subplots(figsize=(6, 0.4 * len(unique_profiles) + 2))
ax2.axis('off')
handles, labels = [], []
for p in unique_profiles:
    handles.append(Line2D([0],[0], marker='o', color='black', markerfacecolor=prof2col[p], markersize=10, linewidth=0))
    labels.append(', '.join(sorted(p)) if p else '(no S)')
ax2.legend(handles, labels, title='Spike profile', loc='upper left', bbox_to_anchor=(0,1), fontsize=12, title_fontsize=14)
plt.tight_layout()
fig2.savefig(leg_pdf, bbox_inches='tight')
fig2.savefig(leg_svg, bbox_inches='tight')
fig2.savefig(leg_png, dpi=600, bbox_inches='tight')
print('Done.')
