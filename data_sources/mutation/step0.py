#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch-fill final sequences into EnzymeML/YAML files and export categorized copies.

What you get
------------
1) In-place update: writes proteins[*].sequence to the original YAMLs.
2) Automatic variant_type setting and sequence cleanup:
   - wildtype: sets variant_type="wildtype", keeps wildtype_sequence, sets mutant_sequence=null
   - mutant: sets variant_type="mutant", sets wildtype_sequence=null, keeps mutant_sequence
3) Reports:
   - report_ok.csv       (WT or successfully mutated)
   - report_issues.csv   (mutation failed -> reference sequence used, or fetch errors)
4) YAML exports:
   - OK YAMLs   -> --ok-dir (default: mutation_exports/ok)
   - Issue YAMLs (reference_only or fetch errors) -> --issues-dir (default: mutation_exports/issues)
   Folder structure is preserved relative to the --input root.

Usage
-----
# Process a folder recursively
python test.py --input ./output/sabio_rk

# Process a single file
python test.py --input "SABIO_..._P52020_F228A_... .yaml"

# Custom export/report dirs
python test.py --input ./output/sabio_rk --ok-dir ./out_ok --issues-dir ./out_issues --report-dir ./reports
"""

import argparse
import csv
import glob
import os
import re
import sys
import time
from typing import Dict, Tuple, List, Optional, Any

import requests
import yaml

# -------- Amino acid maps (HGVS support) --------
THREE_TO_ONE = {
    "Ala": "A","Cys": "C","Asp": "D","Glu": "E","Phe": "F","Gly": "G","His": "H",
    "Ile": "I","Lys": "K","Leu": "L","Met": "M","Asn": "N","Pro": "P","Gln": "Q",
    "Arg": "R","Ser": "S","Thr": "T","Val": "V","Trp": "W","Tyr": "Y","Ter": "*",
    "Sec": "U","Pyl": "O","Xaa": "X"
}

def three_to_one(aa3: str) -> str:
    key = aa3[:1].upper() + aa3[1:].lower()
    if key not in THREE_TO_ONE:
        raise ValueError(f"Unknown 3-letter AA code: {aa3}")
    return THREE_TO_ONE[key]

# -------- Parse variant_description --------
def _tokenize_variants(desc: str) -> List[str]:
    if not desc:
        return []
    s = desc.strip()
    s = re.sub(r'^\s*p\.\[?', '', s, flags=re.IGNORECASE)  # strip "p." and opening '['
    s = re.sub(r'\]$', '', s)                              # strip trailing ']'
    parts = re.split(r'[;,\s]+', s)
    return [p for p in parts if p]

def parse_variant_description(desc: Optional[str]) -> List[Tuple[str, int, str]]:
    """
    Returns a list of (fromAA, pos, toAA).
    Supports:
      - One-letter: 'F228A', 'H338N;A123G'
      - HGVS three-letter: 'p.Cys257Leu', 'p.[Ala123Gly;Cys257Leu]'
      - Ignores '=' tokens
    """
    if not desc:
        return []
    tokens = _tokenize_variants(desc)
    muts: List[Tuple[str, int, str]] = []
    for t in tokens:
        if t == "=":
            continue
        m1 = re.fullmatch(r'([A-Z\*])(\d+)([A-Z\*])', t)
        if m1:
            muts.append((m1.group(1), int(m1.group(2)), m1.group(3)))
            continue
        m2 = re.fullmatch(r'([A-Za-z]{3})(\d+)([A-Za-z]{3}|=)', t)
        if m2:
            if m2.group(3) == "=":
                continue
            aa0 = three_to_one(m2.group(1))
            aa1 = three_to_one(m2.group(3))
            muts.append((aa0, int(m2.group(2)), aa1))
            continue
        raise ValueError(f"Unrecognized mutation token: '{t}' (from '{desc}')")
    # de-duplicate while preserving order
    seen = set()
    out = []
    for a, p, b in muts:
        key = (a, p, b)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out

# -------- UniProt fetch with small cache --------
_seq_cache: Dict[Tuple[str, Optional[int]], str] = {}

def fetch_uniprot_seq(accession: str, isoform: Optional[int] = None, timeout: int = 30) -> str:
    key = (accession, isoform)
    if key in _seq_cache:
        return _seq_cache[key]
    query = f"accession:{accession}" if isoform is None else f"accession:{accession}-{isoform}"
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={query}"
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    text = r.text.strip()
    if not text.startswith(">"):
        raise ValueError(f"Invalid FASTA for {accession} (isoform={isoform})")
    seq = "".join(line.strip() for line in text.splitlines()[1:] if line and not line.startswith(">"))
    if not seq:
        raise ValueError(f"No sequence parsed for {accession} (isoform={isoform})")
    _seq_cache[key] = seq
    return seq

# -------- Helpers --------
def apply_mutations(ref_seq: str, muts: List[Tuple[str, int, str]], strict: bool = True) -> str:
    s = list(ref_seq)
    for aa0, pos, aa1 in muts:
        if not (1 <= pos <= len(s)):
            raise ValueError(f"Position {pos} out of range (len={len(s)})")
        ref = s[pos - 1]
        if strict and ref != aa0:
            raise ValueError(f"Ref mismatch at {pos}: seq has '{ref}', mutation expects '{aa0}'")
        s[pos - 1] = aa1
    return "".join(s)

def looks_like_wt(desc: Optional[str]) -> bool:
    if desc is None:
        return True
    s = desc.strip().lower()
    return s in {"wt", "wildtype", "-", "none", ""}

def list_yaml_paths(root: str) -> List[str]:
    if os.path.isdir(root):
        return sorted(glob.glob(os.path.join(root, "**", "*.yaml"), recursive=True))
    return [root] if root.lower().endswith(".yaml") else []

def ensure_parent_dir(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)

def rel_export_path(src_path: str, input_root: str, out_root: str) -> str:
    # preserve folder structure relative to input_root
    rel = os.path.relpath(src_path, start=input_root)
    dest = os.path.join(out_root, rel)
    ensure_parent_dir(dest)
    return dest

# -------- Core processing (writes YAML in-place + exports copies) --------
def process_yaml_file(path: str,
                      input_root: str,
                      ok_rows: List[List[Any]],
                      issue_rows: List[List[Any]],
                      ok_dir: str,
                      issues_dir: str) -> Tuple[int, int]:
    """
    Returns (num_written, num_attempted).
    Also writes categorized copies to ok_dir / issues_dir.
    """
    with open(path, "r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)


    # Handle empty or invalid YAML files
    if doc is None or not isinstance(doc, dict):
        print(f"[SKIP] {os.path.basename(path)} - Empty or invalid YAML file", file=sys.stderr)
        return (0, 0)

    proteins = doc.get("proteins")
    if not isinstance(proteins, list) or not proteins:
        return (0, 0)

    changed = False
    attempted = 0
    written = 0

    # Track the strongest status for export decision:
    # If any protein is "reference_only" or "no_sequence" -> export to issues_dir,
    # otherwise (all WT/mutated) -> export to ok_dir
    file_has_issue = False

    for i, prot in enumerate(proteins):
        if not isinstance(prot, dict):
            continue

        prot_id = prot.get("id")
        uniprotid = prot.get("uniprotid") or prot.get("uniprot_id") or prot.get("uniprot")
        vd = prot.get("variant_description") or prot.get("variant") or prot.get("mutation")
        iso = prot.get("isoform")
        try:
            iso = int(iso) if iso is not None else None
        except Exception:
            iso = None

        if not uniprotid:
            continue

        attempted += 1
        status = "WT"
        msg = ""
        final_seq = None
        muts_str = ""

        try:
            ref_seq = fetch_uniprot_seq(uniprotid, isoform=iso)
            if looks_like_wt(vd) or not vd:
                # WT: wildtype_sequence has sequence, mutant_sequence is null
                prot["variant_type"] = "wildtype"
                prot["wildtype_sequence"] = ref_seq
                prot["mutant_sequence"] = None
                status = "WT"
                ok_rows.append([path, i, prot_id, uniprotid, iso, vd or "", status, "", len(ref_seq), ""])
            else:
                mutations = parse_variant_description(vd)
                muts_str = ";".join([f"{a}{p}{b}" for a, p, b in mutations])
                try:
                    mutant_seq = apply_mutations(ref_seq, mutations, strict=True)
                    # Successful mutation: variant_type=mutant, keep only mutant_sequence
                    prot["variant_type"] = "mutant"
                    prot["wildtype_sequence"] = None
                    prot["mutant_sequence"] = mutant_seq
                    status = "mutated"
                    ok_rows.append([path, i, prot_id, uniprotid, iso, vd, status, muts_str, len(mutant_seq), ""])
                except Exception as e_mut:
                    # fall back to reference sequence; mutant_sequence is null
                    prot["variant_type"] = "wildtype"
                    prot["wildtype_sequence"] = ref_seq
                    prot["mutant_sequence"] = None
                    status = "reference_only"
                    msg = f"mutation_failed: {e_mut}"
                    file_has_issue = True
                    issue_rows.append([path, i, prot_id, uniprotid, iso, vd, status, muts_str, len(ref_seq), msg])
        except Exception as e_fetch:
            # fetch failed -> cannot write sequence; record issue
            status = "no_sequence"
            msg = f"fetch_failed: {e_fetch}"
            file_has_issue = True
            issue_rows.append([path, i, prot_id, uniprotid, iso, vd or "", status, "", "", msg])
            print(f"[ISSUE] {os.path.basename(path)} :: proteins[{i}] {uniprotid} - {msg}", file=sys.stderr)
            # continue to next protein without writing sequence
            continue

        # Remove old "sequence" field if it exists
        if "sequence" in prot:
            del prot["sequence"]
        
        proteins[i] = prot
        changed = True
        written += 1

        tag = "OK" if status in {"WT", "mutated"} else "ISSUE"
        if status == "WT" or status == "reference_only":
            seq_info = f"variant_type=wildtype, WT len={len(ref_seq)}, MT=null"
        else:  # mutated
            seq_info = f"variant_type=mutant, WT=null, MT len={len(prot['mutant_sequence'])}"
        print(f"[{tag}] {os.path.basename(path)} :: proteins[{i}] {uniprotid} -> "
              f"{seq_info} ({status}{'' if not muts_str else ' | ' + muts_str})")

        time.sleep(0.05)  # polite rate limiting

    # save the in-place YAML if anything changed
    if changed:
        with open(path, "w", encoding="utf-8") as fw:
            yaml.safe_dump(doc, fw, allow_unicode=True, sort_keys=False)

    # export a copy of the (possibly modified) YAML into the appropriate folder
    export_root = issues_dir if file_has_issue else ok_dir
    dest_path = rel_export_path(path, input_root=input_root, out_root=export_root)
    with open(dest_path, "w", encoding="utf-8") as ef:
        yaml.safe_dump(doc, ef, allow_unicode=True, sort_keys=False)

    return (written, attempted)

# -------- Main --------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", "-i", required=True, help="YAML file or folder (recursive)")
    ap.add_argument("--report-dir", default="./mutation_reports", help="Directory for CSV reports")
    ap.add_argument("--ok-dir", default="../../output/mutation_exports/ok", help="Directory for OK YAML exports")
    ap.add_argument("--issues-dir", default="../../output/mutation_exports/issues", help="Directory for Issue YAML exports")
    args = ap.parse_args()
    # Normalize all directory paths to absolute paths to avoid Windows path issues
    args.report_dir = os.path.abspath(args.report_dir)
    args.ok_dir = os.path.abspath(args.ok_dir)
    args.issues_dir = os.path.abspath(args.issues_dir)


    # list files and determine a base root for relative export paths
    paths = list_yaml_paths(args.input)
    if not paths:
        print("No YAML files found.", file=sys.stderr)
        sys.exit(1)
    if os.path.isdir(args.input):
        input_root = os.path.abspath(args.input)
    else:
        input_root = os.path.dirname(os.path.abspath(args.input)) or "."

    os.makedirs(args.report_dir, exist_ok=True)
    os.makedirs(args.ok_dir, exist_ok=True)
    os.makedirs(args.issues_dir, exist_ok=True)

    ok_rows: List[List[Any]] = []
    issue_rows: List[List[Any]] = []

    total_written = 0
    total_attempted = 0

    for p in paths:
        w, a = process_yaml_file(
            path=p,
            input_root=input_root,
            ok_rows=ok_rows,
            issue_rows=issue_rows,
            ok_dir=args.ok_dir,
            issues_dir=args.issues_dir,
        )
        total_written += w
        total_attempted += a

    # write reports
    ok_path = os.path.join(args.report_dir, "report_ok.csv")
    issues_path = os.path.join(args.report_dir, "report_issues.csv")
    headers = [
        "file", "protein_index", "protein_id", "uniprotid", "isoform", "variant_description",
        "status", "mutations", "sequence_length", "message"
    ]

    with open(ok_path, "w", newline="", encoding="utf-8") as f_ok:
        csv.writer(f_ok).writerow(headers)
        csv.writer(f_ok).writerows(ok_rows)

    with open(issues_path, "w", newline="", encoding="utf-8") as f_is:
        csv.writer(f_is).writerow(headers)
        csv.writer(f_is).writerows(issue_rows)

    print(f"\nDone: wrote sequences for {total_written} / attempted {total_attempted} proteins.")
    print(f"Reports:\n  OK     -> {ok_path}\n  Issues -> {issues_path}")
    print(f"Exports:\n  OK     -> {os.path.abspath(args.ok_dir)}\n  Issues -> {os.path.abspath(args.issues_dir)}")

if __name__ == "__main__":
    main()
