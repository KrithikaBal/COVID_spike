#!/usr/bin/env python3
import os
import sys
import pandas as pd
from collections import defaultdict

# ──────────────────────────────────────────────────────────────────────────────
# 1) Paths & sanity checks
# ──────────────────────────────────────────────────────────────────────────────
OUTDIR        = "nextclade_output_april_may"
NEXTCLADE_CSV = os.path.join(OUTDIR, "nextclade.csv")
OUT_CSV       = "muttable.csv"

if not os.path.isfile(NEXTCLADE_CSV):
    print(f"ERROR: {NEXTCLADE_CSV} not found. Run Nextclade first.", file=sys.stderr)
    sys.exit(1)

# ──────────────────────────────────────────────────────────────────────────────
# 2) Load Nextclade CSV & filter for qc=good
# ──────────────────────────────────────────────────────────────────────────────
df = pd.read_csv(NEXTCLADE_CSV, sep=";")
df_good = df[df["qc.overallStatus"] == "good"].copy()
print(f"→ Loaded {len(df_good)} ‘good’ genomes from Nextclade CSV.")

if "Nextclade_pango" not in df_good.columns:
    print("ERROR: Nextclade_pango column missing.", file=sys.stderr)
    sys.exit(1)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Extract each genome’s spike‐mutation set
# ──────────────────────────────────────────────────────────────────────────────
def extract_spike_set(aa_subs):
    """
    From Nextclade’s ‘aaSubstitutions’ string (comma‐separated),
    keep only those that start with "S:" and return as a Python set.
    """
    if pd.isna(aa_subs) or aa_subs.strip() == "":
        return set()
    return set(x.strip() for x in aa_subs.split(",") if x.startswith("S:"))

df_good["spikeProfileSet"] = df_good["aaSubstitutions"].apply(extract_spike_set)

# ──────────────────────────────────────────────────────────────────────────────
# 4) For each lineage, union all spike sets
# ──────────────────────────────────────────────────────────────────────────────
lineage_to_mutset = defaultdict(set)

for _, row in df_good.iterrows():
    lineage = row["Nextclade_pango"]
    muts    = row["spikeProfileSet"]
    lineage_to_mutset[lineage].update(muts)

print(f"→ Found {len(lineage_to_mutset)} distinct lineages among ‘good’ genomes.")

# ──────────────────────────────────────────────────────────────────────────────
# 5) Build a sorted list of ALL spike‐mutations across all lineages
# ──────────────────────────────────────────────────────────────────────────────
all_mutations = sorted({m for s in lineage_to_mutset.values() for m in s})
print(f"→ Total distinct spike mutations observed: {len(all_mutations)}")

# ──────────────────────────────────────────────────────────────────────────────
# 6) Construct a DataFrame: rows=lineages, cols=mutations, 'X' if present
# ──────────────────────────────────────────────────────────────────────────────
# Initialize an empty DataFrame with index=lineages and columns=all_mutations
df_table = pd.DataFrame(
    "", 
    index=sorted(lineage_to_mutset.keys()),
    columns=all_mutations,
)

# Fill in "X" wherever a lineage’s union‐set contains that mutation
for lineage, muts in lineage_to_mutset.items():
    for m in muts:
        df_table.at[lineage, m] = "X"

# Optionally, if you want “Lineage” as its own column instead of index:
df_table.insert(0, "Lineage", df_table.index)

# ──────────────────────────────────────────────────────────────────────────────
# 7) Write to CSV
# ──────────────────────────────────────────────────────────────────────────────
df_table.to_csv(OUT_CSV, index=False)
print(f"→ Wrote mutation‐table to {OUT_CSV}")
