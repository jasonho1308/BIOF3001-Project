#!/usr/bin/env python3
"""
Download GSE74193 beta matrix + metadata
into:  data/geo/raw/

- Downloads .gz files
- Decompresses to .csv and .txt
- Keeps .gz originals (for DVC)
"""

import sys
import gzip
import shutil
import requests
from pathlib import Path
import argparse
import logging
from typing import Optional
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ----------------------------------------------------------------------
LOG = logging.getLogger(__name__)

DEFAULT_OUT = Path("data/geo/raw")

# ----------------------------------------------------------------------
# 1. Beta matrix
BETA_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74193"
    "&format=file&file=GSE74193_GEO_procData.csv.gz"
)
BETA_NAME = "GSE74193_GEO_procData.csv.gz"
BETA_CSV_NAME = "GSE74193_GEO_procData.csv"
BETA_GZ = DEFAULT_OUT / BETA_NAME
BETA_CSV = DEFAULT_OUT / BETA_CSV_NAME

# 2. Series matrix (metadata)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/matrix/"
    "GSE74193_series_matrix.txt.gz"
)
MATRIX_NAME = "GSE74193_series_matrix.txt.gz"
MATRIX_TXT_NAME = "GSE74193_series_matrix.txt"
MATRIX_GZ = DEFAULT_OUT / MATRIX_NAME
MATRIX_TXT = DEFAULT_OUT / MATRIX_TXT_NAME

# ----------------------------------------------------------------------
def _make_session(retries: int = 3, backoff: float = 0.3) -> requests.Session:
    """Create a requests.Session with retry/backoff configured.

    Returns
    -------
    requests.Session
    """
    s = requests.Session()
    retry = Retry(total=retries, backoff_factor=backoff, status_forcelist=(500, 502, 503, 504))
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "biof3001-downloader/1.0"})
    return s


def download_file(url: str, dest: Path, desc: str, session: Optional[requests.Session] = None, overwrite: bool = False) -> None:
    """Download a remote URL to `dest` safely.

    The function downloads into a temporary .part file and renames into place on success.

    Parameters
    ----------
    url : str
        Remote URL to download.
    dest : Path
        Destination file path.
    desc : str
        Human-readable description for printing.
    session : requests.Session | None
        Optional requests session to use (recommended for retries).
    overwrite : bool
        If True, always re-download even when `dest` exists.
    """
    dest = Path(dest)
    if dest.exists() and not overwrite:
        LOG.info("%s already exists: %s", desc, dest.name)
        return

    LOG.info("Downloading %s → %s", desc, dest.name)
    s = session or _make_session()

    tmp = dest.with_suffix(dest.suffix + ".part")
    headers = {}
    mode = "wb"
    downloaded = 0
    # resume support: if a .part file exists, request remaining bytes and append
    if tmp.exists() and not overwrite:
        start = tmp.stat().st_size
        if start > 0:
            headers["Range"] = f"bytes={start}-"
            mode = "ab"
            downloaded = start

    try:
        # Stream to temporary file
        with s.get(url, stream=True, timeout=60, headers=headers) as r:
            r.raise_for_status()
            # content-length is remaining bytes in case of Range request
            remaining = int(r.headers.get("content-length", 0)) or None
            total = (downloaded + remaining) if remaining is not None else None
            with open(tmp, mode) as fh:
                for chunk in r.iter_content(chunk_size=8192):
                    if not chunk:
                        continue
                    fh.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        done = int(50 * downloaded / total)
                        mb = downloaded / 1e6
                        tot_mb = total / 1e6
                        print(f"\r[{'█' * done:<50}] {mb:.1f}/{tot_mb:.1f} MB", end="", flush=True)
        # rename into final destination
        tmp.replace(dest)
        LOG.info("Saved: %s", dest.name)
    except requests.RequestException as e:
        LOG.error("Download failed (%s): %s", desc, e)
        # keep partial file for potential resume; do not delete
        raise
    except Exception as e:
        LOG.error("Download failed (%s): %s", desc, e)
        if tmp.exists():
            try:
                tmp.unlink()
            except Exception:
                pass
        raise

# ----------------------------------------------------------------------
def decompress_gz(gz_path: Path, out_path: Path, desc: str, overwrite: bool = False) -> None:
    """Decompress a .gz file to `out_path`.

    Parameters
    ----------
    gz_path : Path
        Source .gz file.
    out_path : Path
        Destination decompressed file.
    desc : str
        Description printed to logs.
    overwrite : bool
        If True, overwrite existing decompressed file.
    """
    if out_path.exists() and not overwrite:
        LOG.info("Decompressed %s already exists: %s", desc, out_path.name)
        return

    LOG.info("Decompressing %s → %s", gz_path.name, out_path.name)
    with gzip.open(gz_path, "rb") as f_in:
        with open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    LOG.info("Created: %s", out_path.name)

# ----------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Download GSE74193 beta matrix and series matrix into data/geo/raw")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT, help="Output directory")
    parser.add_argument("--no-decompress", action="store_true", help="Skip decompression step")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
    parser.add_argument("--show-preview", action="store_true", help="Show first lines of series matrix after download")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # adjust global paths to chosen outdir
    beta_gz = out_dir / BETA_GZ.name
    beta_csv = out_dir / BETA_CSV.name
    matrix_gz = out_dir / MATRIX_GZ.name
    matrix_txt = out_dir / MATRIX_TXT.name

    LOG.info("Target folder: %s", out_dir.resolve())

    session = _make_session()

    # Step 1: Download .gz files
    try:
        download_file(BETA_URL, beta_gz, "Beta matrix (compressed)", session=session, overwrite=args.overwrite)
        download_file(MATRIX_URL, matrix_gz, "Series matrix (compressed)", session=session, overwrite=args.overwrite)
    except Exception as exc:
        LOG.error("Download(s) failed: %s", exc)
        sys.exit(2)

    # Step 2: Decompress (optional)
    if not args.no_decompress:
        try:
            decompress_gz(beta_gz, beta_csv, "beta matrix", overwrite=args.overwrite)
            decompress_gz(matrix_gz, matrix_txt, "series matrix", overwrite=args.overwrite)
        except Exception as exc:
            LOG.error("Decompression failed: %s", exc)
            sys.exit(3)

    LOG.info("All done! Files created:")
    LOG.info("  • %s  (original)", beta_gz.name)
    LOG.info("  • %s  (decompressed CSV)", beta_csv.name)
    LOG.info("  • %s  (original)", matrix_gz.name)
    LOG.info("  • %s  (decompressed TXT)", matrix_txt.name)
    LOG.info("Full path: %s", out_dir.resolve())

    # Optional: Show first few lines of metadata
    if args.show_preview:
        LOG.info("First 10 lines of series matrix:")
        try:
            with open(matrix_txt, "r", encoding="utf-8", errors="ignore") as f:
                for i, line in enumerate(f):
                    if i >= 10:
                        break
                    print(f"   {line.rstrip()}")
        except Exception as e:
            LOG.warning("Could not preview file: %s", e)

if __name__ == "__main__":
    main()