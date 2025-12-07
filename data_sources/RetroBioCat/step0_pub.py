#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Select YAML files that contain any PMID from a given (built-in) target list.
- Recursively scans .yml/.yaml under --input-dir
- Robustly extracts PMIDs from:
  * pubmed_id / pubmedid / pmid
  * uri/link fields that contain https://pubmed.ncbi.nlm.nih.gov/<PMID>/
  * any string leaf that looks like a PMID / PubMed URL / "PMID:xxxxx"
- Normalizes PMIDs: strip URLs, remove '.0', keep digits only
- Copies (or moves with --move) matched YAMLs to out/by_target_pmids/
- Writes reports:
  * pmid_matches.csv: file, matched_pmids, all_pmids_in_file
  * target_pmids_not_found.csv: pmids from target list not found in any YAML
"""

import argparse
import os
import re
import csv
import sys
import shutil
from typing import Any, Dict, Iterable, List, Set

try:
    import yaml
except ImportError:
    print("Missing dependency: PyYAML. Install with `pip install pyyaml`.", file=sys.stderr)
    sys.exit(1)

# === Your target PMID list (from your earlier message) ===
TARGET_PMIDS: Set[str] = {
    "23075382", "26404833", "29067382", "30601002", "28450969",
    "20455228", "30333895", "28452041", "28010060", "29393988",
    "27547270", "29024400", "28937665", "25372591", "26689856"
}

DIGITS_RE = re.compile(r"\d+")
PUBMED_URL_RE = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", re.IGNORECASE)

def normalize_pmid(value: Any) -> str:
    """Return normalized PMID (digits only) or '' if invalid."""
    if value is None:
        return ""
    s = str(value).strip()

    # Prefer extracting from PubMed URL if present
    m = PUBMED_URL_RE.search(s)
    if m:
        s = m.group(1)

    s = s.strip().rstrip("/")
    if s.lower().startswith("pmid:"):
        s = s[5:].strip()

    if s.endswith(".0"):
        s = s[:-2]

    digits = "".join(DIGITS_RE.findall(s))
    if 6 <= len(digits) <= 10:  # PMID length heuristic
        return digits
    return ""

def iter_yaml_files(root: str) -> Iterable[str]:
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            if os.path.splitext(fn)[1].lower() in {".yml", ".yaml"}:
                yield os.path.join(dirpath, fn)

def safe_load_yaml(path: str) -> Any:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        return yaml.safe_load(f)

def extract_pmids(obj: Any, parent_key: str = "") -> Set[str]:
    """Recursively collect PMIDs from a loaded YAML object."""
    pmids: Set[str] = set()

    if isinstance(obj, dict):
        for k, v in obj.items():
            k_l = str(k).lower()

            # direct fields
            if any(tag in k_l for tag in ("pubmed_id", "pubmedid", "pmid")):
                pm = normalize_pmid(v)
                if pm:
                    pmids.add(pm)

            # uri/link-like fields
            if "uri" in k_l or "link" in k_l:
                pm = normalize_pmid(v)
                if pm:
                    pmids.add(pm)

            # citation-ish blocks
            if any(tag in k_l for tag in ("pub", "citation", "reference", "ref", "bibliograph")):
                pmids |= extract_pmids(v, parent_key=k_l)
                continue

            # generic recursion
            pmids |= extract_pmids(v, parent_key=k_l)

    elif isinstance(obj, list):
        for it in obj:
            pmids |= extract_pmids(it, parent_key=parent_key)

    else:
        # leaf
        if isinstance(obj, (str, int, float)):
            s = str(obj)
            m = PUBMED_URL_RE.search(s)
            if m:
                pmids.add(m.group(1))
            else:
                pm = normalize_pmid(s)
                if pm:
                    pmids.add(pm)

    return pmids

def copy_or_move(src: str, dst_dir: str, move: bool = False) -> str:
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, os.path.basename(src))
    if move:
        return shutil.move(src, dst)
    shutil.copy2(src, dst)
    return dst

def write_csv(path: str, header: List[str], rows: Iterable[Iterable[str]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(list(r))

def main():
    ap = argparse.ArgumentParser(description="Select YAMLs that contain any target PMIDs.")
    ap.add_argument("--input-dir", required=True, help="Folder containing .yml/.yaml")
    ap.add_argument("--out-dir", required=True, help="Output folder")
    ap.add_argument("--move", action="store_true", help="Move files instead of copying")
    args = ap.parse_args()

    in_dir = os.path.abspath(args.input_dir)
    out_dir = os.path.abspath(args.out_dir)
    out_pick = os.path.join(out_dir, "by_target_pmids")
    os.makedirs(out_pick, exist_ok=True)

    if not os.path.isdir(in_dir):
        print(f"[ERROR] input dir not found: {in_dir}", file=sys.stderr)
        sys.exit(2)

    yaml_files = list(iter_yaml_files(in_dir))
    print(f"[INFO] Found {len(yaml_files)} YAML files under {in_dir}")

    file_to_all_pmids: Dict[str, Set[str]] = {}
    file_to_matched_pmids: Dict[str, Set[str]] = {}
    found_targets: Set[str] = set()
    picked_files: List[str] = []

    for path in yaml_files:
        try:
            data = safe_load_yaml(path)
        except Exception as e:
            print(f"[WARN] Failed to load YAML: {path} ({e})")
            continue

        all_pmids = extract_pmids(data)
        file_to_all_pmids[path] = all_pmids
        matched = all_pmids & TARGET_PMIDS
        if matched:
            file_to_matched_pmids[path] = matched
            found_targets |= matched
            dst = copy_or_move(path, out_pick, move=args.move)
            picked_files.append(dst)

    # reports
    matches_csv = os.path.join(out_dir, "pmid_matches.csv")
    rows = []
    for f in sorted(file_to_all_pmids.keys()):
        matched = sorted(file_to_matched_pmids.get(f, set()))
        all_ = sorted(file_to_all_pmids.get(f, set()))
        rows.append((f, ";".join(matched), ";".join(all_)))
    write_csv(matches_csv, header=["file", "matched_pmids", "all_pmids_in_file"], rows=rows)
    print(f"[INFO] Wrote {matches_csv}")

    missing_targets = sorted(list(TARGET_PMIDS - found_targets))
    missing_csv = os.path.join(out_dir, "target_pmids_not_found.csv")
    write_csv(missing_csv, header=["pmid"], rows=((pm,) for pm in missing_targets))
    print(f"[INFO] Wrote {missing_csv}")

    print(f"[INFO] {'Moved' if args.move else 'Copied'} {len(picked_files)} YAML files to {out_pick}")
    print("[DONE]")

if __name__ == "__main__":
    main()
