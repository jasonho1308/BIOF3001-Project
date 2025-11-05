#!/usr/bin/env python3
"""Stream-extract beta-values from the large gzipped methylation matrix.

Writes a tab-separated gzipped file with the probe ID in the first column and
all sample columns preserved. This script processes the entire file, filtering out rows with NA values.

Output: data/geo/GSE213478_betas_only.tsv.gz
"""
from pathlib import Path
import gzip
import csv
import sys

PROJECT = Path.cwd()
INP = PROJECT / 'data' / 'geo' / 'GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz'
OUTP = PROJECT / 'data' / 'geo' / 'GSE213478_betas_only.tsv.gz'

def main():
    if not INP.exists():
        print('Input file not found:', INP, file=sys.stderr)
        sys.exit(2)

    # Stream read and write
    with gzip.open(INP, 'rt', newline='') as inf, gzip.open(OUTP, 'wt', newline='') as outf:
        reader = csv.reader(inf, delimiter='\t')  # Use tab as delimiter
        writer = csv.writer(outf, delimiter='\t', lineterminator='\n')
        
        header = next(reader)
        # Normalize first column header to 'probe' if empty
        if header[0].strip(' \"') == '':
            header[0] = 'probe'
        writer.writerow(header)

        # Stream lines and write, filtering out rows with NA values
        n = 0
        for row in reader:
            # Check for NA values (assuming NA or empty strings are to be filtered)
            if any(value.strip() == '' or value.lower() == 'na' for value in row):
                continue  # Skip this row if it contains any NA values

            writer.writerow(row)
            n += 1
            if n % 500000 == 0:
                print(f'Processed {n:,} rows...', file=sys.stderr)

    print('Wrote:', OUTP, '- rows written (approx):', n)

if __name__ == '__main__':
    main()