#!/usr/bin/env python3
"""DEPRECATED wrapper for the renamed script.

This project renamed `scripts/analysis_GEO.py` to
`scripts/extract_meta_GEO.py`. This small shim prints an informational message
and exits so callers see the new script name. If you want to run the
functionality, use `python scripts/extract_meta_GEO.py` instead.
"""
import sys


def main():
    print("The script has been renamed to 'scripts/extract_meta_GEO.py'.")
    print("Please run: python scripts/extract_meta_GEO.py --series <series_path> --outdir <outdir>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
