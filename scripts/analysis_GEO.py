#!/usr/bin/env python3
"""Reproduce the extraction/enrichment from the notebook `scripts/analysis_GEO.ipynb`.

This script reads a GEO series_matrix file, parses the !Sample_ header lines,
extracts sample metadata, expands free-text `characteristics` into key:value
columns, normalizes common fields (age, sex, tissue, donor id), and writes two
CSV outputs into the specified output directory:

- gse213478_metadata_parsed.csv (basic parsed fields)
- gse213478_metadata_parsed_enriched.csv (characteristics expanded)

Usage:
    python scripts/analysis_GEO.py --series data/geo/raw/GSE213478_series_matrix.txt \
        --outdir data/geo/processed
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict

import pandas as pd


def read_series_matrix(series_path: Path) -> str:
    if not series_path.exists():
        raise FileNotFoundError(f"Series matrix not found: {series_path}")
    return series_path.read_text(encoding="utf-8")


def parse_sample_headers(text: str) -> Dict[str, list]:
    lines = text.splitlines()
    sample_lines = [L for L in lines if L.startswith("!Sample_")]
    sample_fields: Dict[str, list] = {}
    for L in sample_lines:
        parts = L.split("\t")
        key = parts[0].lstrip("!")
        vals = [v.strip().strip('"') for v in parts[1:]]
        sample_fields.setdefault(key, []).append(vals)
    return sample_fields


def safe_col(k: str) -> str:
    s = k.replace("Sample_", "")
    return re.sub(r"[^0-9a-zA-Z]+", "_", s).strip("_")


def build_records(sample_fields: Dict[str, list]):
    # Determine number of samples
    if "Sample_geo_accession" in sample_fields:
        n = len(sample_fields["Sample_geo_accession"][0])
    elif "Sample_title" in sample_fields:
        n = len(sample_fields["Sample_title"][0])
    else:
        raise RuntimeError("No sample id lines found in the series_matrix header")

    records = []
    for i in range(n):
        rec = {}
        rec["GSM"] = sample_fields.get("Sample_geo_accession", [[""] * n])[0][i] if "Sample_geo_accession" in sample_fields else ""
        rec["title"] = sample_fields.get("Sample_title", [[""] * n])[0][i] if "Sample_title" in sample_fields else ""
        # Aggregate characteristics lines
        char_parts = []
        for k, vlist in sample_fields.items():
            if k.startswith("Sample_characteristics"):
                for line_vals in vlist:
                    if i < len(line_vals) and line_vals[i] != "":
                        char_parts.append(line_vals[i])
        rec["characteristics"] = " | ".join(char_parts)

        # capture other Sample_* fields
        for k, vlist in sample_fields.items():
            if k in ("Sample_geo_accession", "Sample_title") or k.startswith("Sample_characteristics"):
                continue
            try:
                val = vlist[0][i] if i < len(vlist[0]) else ""
            except Exception:
                val = ""
            if val != "":
                rec[safe_col(k)] = val
        records.append(rec)
    return records


def parse_age_raw(ar: str):
    if not ar:
        return "", None
    if "-" in ar:
        try:
            a, b = [float(x) for x in ar.split("-")[:2]]
            return ar, (a + b) / 2
        except Exception:
            return ar, None
    try:
        v = float(ar)
        return ar, v
    except Exception:
        return ar, None


def regex_extract_basic(records):
    re_tissue = re.compile(r"(?i)tissue: *([^|;]+)")
    re_age = re.compile(r"(?i)age: *([0-9]+(?:-[0-9]+)?)")
    re_sex = re.compile(r"(?i)sex: *([12]|[MFmf]|male|female)")
    re_subject_collab = re.compile(r"(?i)collaborator_participant_id: *([A-Za-z0-9_-]+)")
    re_subject_part = re.compile(r"(?i)participant_id: *([A-Za-z0-9_-]+)")
    re_subject_generic = re.compile(r"(?i)(?:subjectid|subject_id|patient|patient_id|donor): *([A-Za-z0-9_-]+)")

    out = []
    for r in records:
        ch = r.get("characteristics", "") or ""
        tissue = ""
        age_raw = ""
        age_num = None
        sex = ""
        donor = ""
        m = re_tissue.search(ch)
        if m:
            tissue = m.group(1).strip()
        m = re_age.search(ch)
        if m:
            age_raw, age_num = parse_age_raw(m.group(1).strip())
        m = re_sex.search(ch)
        if m:
            s = m.group(1).strip()
            if s in ("1", "2"):
                sex = "M" if s == "1" else "F"
            else:
                s2 = s[0].upper() if s else ""
                if s2 in ("M", "F"):
                    sex = s2
        m = re_subject_collab.search(ch)
        if m:
            donor = m.group(1).strip()
        else:
            m = re_subject_part.search(ch)
            if m:
                donor = m.group(1).strip()
            else:
                m = re_subject_generic.search(ch)
                if m:
                    donor = m.group(1).strip()

        rec = {
            "GSM": r.get("GSM", ""),
            "title": r.get("title", ""),
            "tissue": tissue,
            "age_raw": age_raw,
            "age_num": age_num,
            "sex": sex,
            "donor_id": donor,
            "characteristics": ch,
        }
        for k, v in r.items():
            if k in ("GSM", "title", "characteristics"):
                continue
            if k not in rec:
                rec[k] = v
        out.append(rec)
    df = pd.DataFrame(out)
    df["age_num"] = pd.to_numeric(df.get("age_num", pd.Series([None] * len(df))), errors="coerce")
    return df


def normalize_colname(k: str) -> str:
    k = k.strip()
    k = re.sub(r"(?i)collaborator[_ ]?participant[_ ]?id", "collaborator_participant_id", k)
    k = re.sub(r"(?i)participant[_ ]?id", "participant_id", k)
    k = re.sub(r"[^0-9a-zA-Z]+", "_", k).strip("_").lower()
    return k


def parse_characteristics_text(s: str) -> Dict[str, str]:
    if not s or pd.isna(s):
        return {}
    parts = re.split(r"\s*\|\s*|\s*;\s*", s)
    out: Dict[str, str] = {}
    for p in parts:
        if not p:
            continue
        if ":" in p:
            k, v = p.split(":", 1)
        elif "=" in p:
            k, v = p.split("=", 1)
        else:
            k, v = "note", p
        k = normalize_colname(k)
        v = v.strip().strip('"').strip()
        if k in out:
            out[k] = out[k] + " | " + v
        else:
            out[k] = v
    return out


def coerce_age_from_fields(row):
    for candidate in ("age", "age_raw", "age_y", "age_years"):
        if candidate in row and pd.notna(row[candidate]) and str(row[candidate]).strip() != "":
            ar = str(row[candidate]).strip()
            if "-" in ar:
                try:
                    a, b = [float(x) for x in ar.split("-")[:2]]
                    return (ar, (a + b) / 2)
                except Exception:
                    return (ar, None)
            try:
                v = float(re.sub(r"[^0-9.]", "", ar))
                return (ar, v)
            except Exception:
                return (ar, None)
    return (row.get("age_raw", ""), row.get("age_num", None))


def coerce_sex_from_fields(row):
    for candidate in ("sex", "gender"):
        if candidate in row and pd.notna(row[candidate]) and str(row[candidate]).strip() != "":
            s = str(row[candidate]).strip()
            if s in ("1", "2"):
                return "M" if s == "1" else "F"
            s2 = s[0].upper()
            if s2 in ("M", "F"):
                return s2
    return row.get("sex", "")


def normalize_tissue(t):
    if pd.isna(t) or str(t).strip() == "":
        return ""
    s = str(t).strip()
    s = re.sub(r"\s+", " ", s)
    return s


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--series", default="data/geo/raw/GSE213478_series_matrix.txt", help="Path to series_matrix.txt")
    p.add_argument("--outdir", default="data/geo/processed", help="Directory to write outputs")
    args = p.parse_args()

    series_path = Path(args.series)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    text = read_series_matrix(series_path)
    sample_fields = parse_sample_headers(text)
    records = build_records(sample_fields)

    # Basic parsed dataframe (regex extracts)
    df_basic = regex_extract_basic(records)

    # Enriched: expand characteristics into columns
    parsed_dicts = [parse_characteristics_text(ch) for ch in df_basic["characteristics"].fillna("")]
    parsed_df = pd.DataFrame(parsed_dicts)

    for col in parsed_df.columns:
        if col in df_basic.columns:
            df_basic[col] = df_basic[col].astype("object").fillna(parsed_df[col])
        else:
            df_basic[col] = parsed_df[col]

    age_results = df_basic.apply(coerce_age_from_fields, axis=1)
    df_basic["age_raw_parsed"] = [x[0] for x in age_results]
    df_basic["age_num_parsed"] = [x[1] for x in age_results]
    df_basic["sex_parsed"] = df_basic.apply(coerce_sex_from_fields, axis=1)
    df_basic["tissue"] = df_basic["tissue"].fillna("").apply(normalize_tissue)

    out_basic = outdir / "gse213478_metadata_parsed.csv"
    out_enriched = outdir / "gse213478_metadata_parsed_enriched.csv"
    df_basic.to_csv(out_basic, index=False)
    df_basic.to_csv(out_enriched, index=False)

    print(f"Wrote parsed CSV: {out_basic}")
    print(f"Wrote enriched CSV: {out_enriched}")


if __name__ == "__main__":
    main()
