#!/usr/bin/env python3
"""Stream-extract beta-values from the large gzipped methylation matrix.

Writes a tab-separated gzipped file with the probe ID in the first column and
all sample columns preserved. The script detects the input delimiter (comma
or tab) and streams the file to avoid using large memory.

Output: data/geo/GSE213478_betas_only.tsv.gz
"""
from pathlib import Path
import gzip
import csv
import sys

PROJECT = Path(__file__).resolve().parents[1]
INP = PROJECT / 'data' / 'geo' / 'GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz'
OUTP = PROJECT / 'data' / 'geo' / 'GSE213478_betas_only.tsv.gz'

def detect_delim(sample_line: str) -> str:
    # crude heuristic: comma vs tab
    if sample_line.count(',') > sample_line.count('\t'):
        return ','
    return '\t'

def main():
    if not INP.exists():
        print('Input file not found:', INP, file=sys.stderr)
        sys.exit(2)

    # open once to read the first line to detect delimiter
    with gzip.open(INP, 'rt', newline='') as fh:
        first = fh.readline()
        if not first:
            print('Empty input', file=sys.stderr)
            sys.exit(3)
        delim = detect_delim(first)

    print('Detected delimiter:', repr(delim))

    # stream read and write
    with gzip.open(INP, 'rt', newline='') as inf, gzip.open(OUTP, 'wt', newline='') as outf:
        reader = csv.reader(inf, delimiter=delim)
        writer = csv.writer(outf, delimiter='\t', lineterminator='\n')
        header = next(reader)
        # normalize first column header to 'probe' if empty
        if header[0].strip(' \"') == '':
            header[0] = 'probe'
        writer.writerow(header)
        # stream lines and write
        n = 0
        for row in reader:
            writer.writerow(row)
            n += 1
            if n % 500000 == 0:
                print(f'Processed {n:,} rows...', file=sys.stderr)

    print('Wrote:', OUTP, '- rows written (approx):', n)

if __name__ == '__main__':
    main()
