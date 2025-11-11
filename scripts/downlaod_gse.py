#!/usr/bin/env python3
"""
Download GSE74193 beta matrix + metadata
into:  data/geo/raw/

- Downloads .gz files
- Decompresses to .csv and .txt
- Keeps .gz originals (for DVC)
"""

import os
import sys
import gzip
import shutil
import requests
from pathlib import Path

# ----------------------------------------------------------------------
# Output directory (relative to script location)
OUT_DIR = Path("data/geo/raw")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# 1. Beta matrix
BETA_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74193"
    "&format=file&file=GSE74193_GEO_procData.csv.gz"
)
BETA_GZ = OUT_DIR / "GSE74193_GEO_procData.csv.gz"
BETA_CSV = OUT_DIR / "GSE74193_GEO_procData.csv"

# 2. Series matrix (metadata)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/matrix/"
    "GSE74193_series_matrix.txt.gz"
)
MATRIX_GZ = OUT_DIR / "GSE74193_series_matrix.txt.gz"
MATRIX_TXT = OUT_DIR / "GSE74193_series_matrix.txt"

# ----------------------------------------------------------------------
def download_file(url: str, dest: Path, desc: str) -> None:
    dest = Path(dest)
    if dest.exists():
        print(f"[OK] {desc} already exists: {dest.name}")
        return

    print(f"[Download] {desc} → {dest.name}")
    try:
        headers = {}
        if dest.exists():
            headers["Range"] = f"bytes={dest.stat().st_size}-"

        r = requests.get(url, stream=True, headers=headers, timeout=60)
        r.raise_for_status()
    except Exception as e:
        print(f"[ERROR] Download failed: {e}")
        sys.exit(1)

    total = int(r.headers.get("content-length", 0)) or 0
    mode = "ab" if "Range" in headers else "wb"
    downloaded = dest.stat().st_size if mode == "ab" else 0

    with open(dest, mode) as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if total > 0:
                    done = int(50 * downloaded / total)
                    mb = downloaded / 1e6
                    tot_mb = total / 1e6
                    print(f"\r[{'█' * done:<50}] {mb:.1f}/{tot_mb:.1f} MB", end="", flush=True)
    print(f"\n[OK] Saved: {dest.name}")

# ----------------------------------------------------------------------
def decompress_gz(gz_path: Path, out_path: Path, desc: str) -> None:
    if out_path.exists():
        print(f"[OK] Decompressed {desc} already exists: {out_path.name}")
        return

    print(f"[Decompress] {gz_path.name} → {out_path.name}")
    with gzip.open(gz_path, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"[OK] Created: {out_path.name}")

# ----------------------------------------------------------------------
def main() -> None:
    print(f"Target folder: {OUT_DIR.resolve()}\n")

    # Step 1: Download .gz files
    download_file(BETA_URL, BETA_GZ, "Beta matrix (compressed)")
    download_file(MATRIX_URL, MATRIX_GZ, "Series matrix (compressed)")

    # Step 2: Decompress
    decompress_gz(BETA_GZ, BETA_CSV, "beta matrix")
    decompress_gz(MATRIX_GZ, MATRIX_TXT, "series matrix")

    print("\nAll done! Files created:")
    print(f"   • {BETA_GZ.name}  (original)")
    print(f"   • {BETA_CSV.name} (decompressed CSV)")
    print(f"   • {MATRIX_GZ.name}  (original)")
    print(f"   • {MATRIX_TXT.name} (decompressed TXT)")
    print(f"\nFull path: {OUT_DIR.resolve()}")

    # Optional: Show first few lines of metadata
    print("\nFirst 10 lines of series matrix:")
    try:
        with open(MATRIX_TXT, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if i >= 10:
                    break
                print(f"   {line.rstrip()}")
    except Exception as e:
        print(f"   [Could not preview: {e}]")

if __name__ == "__main__":
    main()