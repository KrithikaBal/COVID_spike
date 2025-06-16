#!/usr/bin/env python3
"""
combine_nextclade.py

Script to concatenate two Nextclade CSV outputs (semicolon-delimited) into one combined CSV.

Usage:
    python combine_nextclade.py

This will read `nextclade_1.csv` and `nextclade.csv` from the
`nextclade_output_Mar_May` directory and write `nextclade_combined.csv`.
"""
import pandas as pd
import os
import sys

# Directory containing the CSV files
dir_path = os.path.join(os.path.dirname(__file__), 'nextclade_output_Mar_May')
# CSV filenames to combine
input_files = ['nextclade_1.csv', 'nextclade.csv']

def main():
    dataframes = []
    for fname in input_files:
        path = os.path.join(dir_path, fname)
        if not os.path.isfile(path):
            print(f"Error: file not found: {path}", file=sys.stderr)
            sys.exit(1)
        print(f"Loading {path}...")
        df = pd.read_csv(path, sep=';')
        dataframes.append(df)
    # Concatenate, preserving header from first file
    combined = pd.concat(dataframes, ignore_index=True)
    # Output path
    out_path = os.path.join(dir_path, 'nextclade_combined.csv')
    print(f"Writing combined CSV to {out_path}...")
    combined.to_csv(out_path, sep=';', index=False)
    print(f"Done: combined {len(input_files)} files into {len(combined)} rows.")

if __name__ == '__main__':
    main()
