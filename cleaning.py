import pandas as pd
from Bio import SeqIO

# 1) Load Nextclade’s CSV.
#    (Note: this file uses semicolons as separators, not commas.)
df = pd.read_csv("nextclade_output_april_may/nextclade.csv", sep=";")

# 2) Filter for qc.overallStatus == "good"
good_df = df[df["qc.overallStatus"] == "good"].copy()

print(f"Total sequences: {len(df)}")
print(f"Sequences flagged ‘good’: {len(good_df)}")

# 3) Grab the list of sequence names to keep
good_names = set(good_df["seqName"].tolist())

# 4) Now filter the aligned FASTA so only those “good” tips remain.
#    Nextclade’s aligned FASTA file is typically named “nextclade.aligned.fasta”.
in_fasta  = "nextclade_output_april_may/nextclade.aligned.fasta"
out_fasta = "nextclade_output_april_may/good_aligned.fasta"

with open(in_fasta) as infile, open(out_fasta, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        # record.id exactly matches the seqName in the CSV
        if record.id in good_names:
            SeqIO.write(record, outfile, "fasta")

print(f"Wrote {sum(1 for _ in SeqIO.parse(out_fasta,'fasta'))} good sequences to {out_fasta}")
