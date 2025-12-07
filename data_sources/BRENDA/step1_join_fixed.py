#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BRENDA -> EnzymeML exporter (TN-only references; NSP-only equations)
Changes in this version:
- TN (kcat) values are ALWAYS exported in 1/s. If '1/min' is detected, value is divided by 60.
  Parenthetical buffers like "50 mM" no longer affect kcat units.
- AC (activators): small-molecule id is just "activator_<Name>", and name is the bare compound
  (e.g., "Isopropanol"); parenthetical info goes to notes only (never in id/name).
- protein.name comes only from PR (organism), otherwise fallback "protein#<protein_idx>".
- Equations are collected ONLY from NSP, deduplicated, with parentheses + brace tokens merged into notes.
- Top-level references include ONLY those actually used by TN (kcat) entries.
- ST (tissues) now correctly filtered by protein_idx with segment-aware note extraction.
- IDs: parameter.id 和 measurement.id 包含 uniprotid；若缺失，则省略该片段（不再使用 NA）。
"""

import argparse
import io
import json
import re
import sys
import zipfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False
    print("Warning: PyYAML not installed. Install: pip install pyyaml")


# ---------------------------- Regex ----------------------------
PROTIDX_RE = re.compile(r"#\s*([0-9]+(?:,[0-9]+)*)\s*#")          # #10,39#
ANGLE_GROUP_RE = re.compile(r"<([^<>]+)>")                         # <49,95>
ALL_PARENS_RE = re.compile(r"\([^()]*\)")                          # ( ... )
CURLY_RE = re.compile(r"\{([^{}]+)\}")                             # {substrate}
PH_RE = re.compile(r"pH\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)
TEMP_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*°?\s*[Cc]")
MW_MULT_RE = re.compile(r"(\d+)\s*\*\s*(\d+)")
UNIPROT_IN_BRACES = re.compile(r"\{([A-Z][0-9][A-Z0-9]{3,}[0-9](?:[A-Z0-9]{0,2})?)")
PUBMED_RE = re.compile(r"Pubmed:(\d+)", re.IGNORECASE)


# ---------------------------- Helpers ----------------------------
def sanitize_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def extract_leading_pidx_list(value: str) -> List[str]:
    m = PROTIDX_RE.search(value)
    if not m:
        return []
    return [x.strip() for x in m.group(1).split(",") if x.strip()]


def extract_all_pidxs(text: str) -> List[str]:
    pidxs: List[str] = []
    for g in PROTIDX_RE.findall(text):
        pidxs.extend([x.strip() for x in g.split(",") if x.strip()])
    return list(dict.fromkeys(pidxs))


def extract_parenthetical_notes(raw: str) -> Optional[str]:
    """Concatenate all (...) contents; remove #..# and <..> tokens; compact whitespace."""
    parts = [m.group(1).strip() for m in re.finditer(r"\(([^()]*)\)", raw)]
    if not parts:
        return None
    text = "; ".join(parts)
    text = PROTIDX_RE.sub("", text)
    text = ANGLE_GROUP_RE.sub("", text)
    return re.sub(r"\s+", " ", text).strip(" ;,") or None


def extract_substrate(value: str) -> Optional[str]:
    m = CURLY_RE.search(value)
    return m.group(1).strip() if m else None


def extract_value_and_unit(text: str) -> Tuple[Optional[float], Optional[str]]:
    """Generic parser for concentration-like or rate-like values."""
    s = PROTIDX_RE.sub("", text, count=1)
    s = CURLY_RE.sub("", s, count=1)
    m = re.search(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", s)
    value = None
    if m:
        try:
            value = float(m.group(1))
        except Exception:
            value = None
    unit = None
    if re.search(r"(s-1|s⁻¹|1/s)\b", s):
        unit = "1/s"
    elif re.search(r"(min-1|min⁻¹|1/min)\b", s):
        unit = "1/min"
    elif re.search(r"\b(µM|uM)\b", s):
        unit = "µM"
    elif re.search(r"\bmM\b", s):
        unit = "mM"
    elif re.search(r"\bnM\b", s):
        unit = "nM"
    elif re.search(r"\bpM\b", s):
        unit = "pM"
    elif re.search(r"(^|\s)M(\s|$)", s):
        unit = "M"
    return value, unit


def parse_kcat_value_unit_brenda(text: str) -> Tuple[Optional[float], str]:
    """
    TN-specific (kcat) parser for BRENDA:
    - Remove parentheses to avoid buffer concentrations (e.g., '50 mM') polluting the unit.
    - Extract numeric value; default unit is ALWAYS '1/s'.
    - If an explicit '1/min' is present, convert value to per-second (divide by 60).
    """
    s = PROTIDX_RE.sub("", text, count=1)
    s = CURLY_RE.sub("", s, count=1)
    s_no_paren = ALL_PARENS_RE.sub("", s)

    # numeric value
    m = re.search(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", s_no_paren)
    value: Optional[float] = None
    if m:
        try:
            value = float(m.group(1))
        except Exception:
            value = None

    # detect per-minute and convert
    if re.search(r"(min-1|min⁻¹|1/min)\b", s_no_paren) and value is not None:
        value = value / 60.0

    # unify to 1/s
    return value, "1/s"


def parse_ph_temp(text: str) -> Tuple[Optional[float], Optional[float]]:
    ph = None
    tm = None
    pm = PH_RE.search(text)
    if pm:
        try:
            ph = float(pm.group(1))
        except Exception:
            ph = None
    tm_ = TEMP_RE.search(text)
    if tm_:
        try:
            tm = float(tm_.group(1))
        except Exception:
            tm = None
    return ph, tm


def parse_range(text: str) -> Tuple[Optional[float], Optional[float]]:
    """
    Parse a range string like '7-9.6' or '7 - 9.6' into lower and upper bounds.
    Returns (lower_bound, upper_bound) or (None, None) if parsing fails.
    """
    if not text:
        return None, None
    
    # Match patterns like: 7-9.6, 7 - 9.6, 7–9.6 (em dash), 7—9.6 (em dash)
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*[-\u2013\u2014]\s*([0-9]+(?:\.[0-9]+)?)", text)
    if m:
        try:
            lower = float(m.group(1))
            upper = float(m.group(2))
            return lower, upper
        except Exception:
            pass
    
    return None, None


def parse_per_ref_segments(raw: str) -> List[Dict[str, Any]]:
    """
    Parse ';'-separated segments inside parentheses.
    Each segment may look like: "#39# text <95>" or "#10,12# text <49>".
    Returns list of dicts with fields: pidxs, ref_idx, note, pH, temperature.
    """
    segs: List[Dict[str, Any]] = []
    for m in re.finditer(r"\(([^()]*)\)", raw):
        inner = m.group(1)
        for part in [p.strip() for p in inner.split(";") if p.strip()]:
            ids = PROTIDX_RE.findall(part)
            pidxs: List[str] = []
            for g in ids:
                pidxs.extend([x.strip() for x in g.split(",") if x.strip()])
            if not pidxs:
                continue
            idxs = ANGLE_GROUP_RE.findall(part)
            ref_idx = None
            if idxs:
                last = idxs[-1]
                ref_idx = last.split(",")[0].strip()
            note_txt = PROTIDX_RE.sub("", part)
            note_txt = ANGLE_GROUP_RE.sub("", note_txt)
            note_txt = note_txt.strip(" ,.;")
            pH, tC = parse_ph_temp(note_txt)
            segs.append({
                "pidxs": list(dict.fromkeys(pidxs)),
                "ref_idx": ref_idx,
                "note": note_txt if note_txt else None,
                "pH": pH,
                "temperature": tC
            })
    return segs


def extract_trailing_refs_outside_parens(raw: str) -> List[str]:
    no_paren = ALL_PARENS_RE.sub("", raw)
    groups = ANGLE_GROUP_RE.findall(no_paren)
    refs: List[str] = []
    for g in groups:
        refs.extend([x.strip() for x in g.split(",") if x.strip()])
    return list(dict.fromkeys(refs))


# ------------------------- Data classes -------------------------
@dataclass
class AuxData:
    organism: Optional[str] = None
    uniprot_id: Optional[str] = None
    synonyms: List[str] = field(default_factory=list)
    systematic_name: Optional[str] = None
    recommended_name: Optional[str] = None

    km_list: List[Dict] = field(default_factory=list)
    kkm_list: List[Dict] = field(default_factory=list)
    ki_list: List[Dict] = field(default_factory=list)

    activators: List[Dict] = field(default_factory=list)   # AC
    inhibitors: List[Dict] = field(default_factory=list)   # IN
    cofactors: List[Dict] = field(default_factory=list)    # CF
    metals: List[Dict] = field(default_factory=list)       # ME

    equations: List[Dict] = field(default_factory=list)    # NSP only

    tissues: List[str] = field(default_factory=list)       # ST

    ph_optimum: Optional[float] = None
    ph_range: Optional[str] = None
    ph_stability: Optional[str] = None

    temp_optimum: Optional[float] = None
    temp_range: Optional[str] = None
    temp_stability: Optional[str] = None
    storage_stability: Optional[str] = None
    oxygen_stability: Optional[str] = None

    protein_molecular_weight: Optional[str] = None
    protein_notes: Optional[str] = None

    post_trans_mods: List[Dict] = field(default_factory=list)  # PM
    mutations: List[Dict] = field(default_factory=list)        # EN

    specific_activities: List[Dict] = field(default_factory=list)  # SA

    used_refs: set = field(default_factory=set)  # local bookkeeping only


# --------------------- Load BRENDA data ---------------------
def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    if len(df.columns) >= 3:
        df.columns = ['REC_ID', 'TAG', 'VALUE_RAW'][:len(df.columns)]
    elif len(df.columns) == 2:
        df.columns = ['REC_ID', 'TAG']
    return df


def load_data_from_zip(zip_path: Path) -> Dict[str, pd.DataFrame]:
    data = {}
    with zipfile.ZipFile(zip_path, 'r') as zf:
        for name in zf.namelist():
            if not name.endswith(".csv"):
                continue
            with zf.open(name) as f:
                try:
                    df = pd.read_csv(io.TextIOWrapper(f, 'utf-8'), sep=None, engine='python',
                                     dtype=str, on_bad_lines='skip')
                    df = standardize_columns(df)
                    if not df.empty:
                        key = df['TAG'].iloc[0] if 'TAG' in df.columns else Path(name).stem
                        key = str(key).strip().upper()
                        data.setdefault(key, pd.DataFrame(columns=df.columns))
                        data[key] = pd.concat([data[key], df], ignore_index=True)
                except Exception as e:
                    print(f"  Warning: cannot load {name}: {e}")
    return data


def load_data_from_dir(dir_path: Path) -> Dict[str, pd.DataFrame]:
    data = {}
    for p in dir_path.glob("*.csv"):
        try:
            df = pd.read_csv(p, sep=None, engine='python', dtype=str, on_bad_lines='skip')
            df = standardize_columns(df)
            if not df.empty:
                key = df['TAG'].iloc[0] if 'TAG' in df.columns else p.stem
                key = str(key).strip().upper()
                data.setdefault(key, pd.DataFrame(columns=df.columns))
                data[key] = pd.concat([data[key], df], ignore_index=True)
        except Exception as e:
            print(f"  Warning: cannot load {p}: {e}")
    return data


# --------------------- Schema completion utils -----------------
def load_schema(schema_path: Path) -> Dict[str, Any]:
    if not HAVE_YAML:
        raise RuntimeError("PyYAML is required. Install: pip install pyyaml")
    with open(schema_path, "r", encoding="utf-8") as f:
        sch = yaml.safe_load(f)
    return sch.get("classes", {})


def build_class_template(class_name: str, classes: Dict[str, Any], cache: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    if class_name in cache:
        return cache[class_name].copy()
    cls = classes.get(class_name)
    if not cls:
        cache[class_name] = {}
        return {}
    tmpl = {}
    for attr in (cls.get("attributes") or {}):
        tmpl[attr] = None
    cache[class_name] = tmpl.copy()
    return tmpl


def _range_is_class(meta: Dict[str, Any], classes: Dict[str, Any]) -> Optional[str]:
    rng = meta.get("range")
    if rng and rng in classes:
        return rng
    return None


def fill_instance_to_schema(class_name: str, instance: Any, classes: Dict[str, Any], cache: Dict[str, Dict[str, Any]]) -> Any:
    cls = classes.get(class_name)
    if not cls:
        return instance
    attrs = cls.get("attributes") or {}
    if not isinstance(instance, dict):
        return build_class_template(class_name, classes, cache)
    result = build_class_template(class_name, classes, cache)
    for k, v in instance.items():
        result[k] = v
    for attr_name, meta in attrs.items():
        rng_class = _range_is_class(meta, classes)
        multival = bool(meta.get("multivalued"))
        val = result.get(attr_name)
        if rng_class:
            if multival:
                if isinstance(val, list):
                    result[attr_name] = [
                        fill_instance_to_schema(rng_class, el if isinstance(el, dict) else {}, classes, cache)
                        for el in val
                    ]
                else:
                    result[attr_name] = None
            else:
                if isinstance(val, dict):
                    result[attr_name] = fill_instance_to_schema(rng_class, val, classes, cache)
                else:
                    result[attr_name] = None
        else:
            if attr_name not in result:
                result[attr_name] = None
    return result


def complete_doc_to_schema(doc: Dict[str, Any], classes: Dict[str, Any]) -> Dict[str, Any]:
    root = "EnzymeMLDocument"
    cache: Dict[str, Dict[str, Any]] = {}
    completed = fill_instance_to_schema(root, doc, classes, cache)
    for attr_name, meta in (classes[root].get("attributes") or {}).items():
        rng_class = _range_is_class(meta, classes)
        multival = bool(meta.get("multivalued"))
        val = completed.get(attr_name)
        if rng_class and multival and isinstance(val, list):
            completed[attr_name] = [
                fill_instance_to_schema(rng_class, el if isinstance(el, dict) else {}, classes, cache)
                for el in val
            ]
        elif rng_class and not multival and isinstance(val, dict):
            completed[attr_name] = fill_instance_to_schema(rng_class, val, classes, cache)
        elif attr_name not in completed:
            completed[attr_name] = None
    return completed


def clean_unwanted_fields(doc: Dict[str, Any]) -> Dict[str, Any]:
    """
    Remove unwanted fields from the completed document.
    Specifically removes: fit, stderr, constant from parameters, and equations from document
    """
    # Clean parameters
    if doc.get("parameters") and isinstance(doc["parameters"], list):
        for param in doc["parameters"]:
            if isinstance(param, dict):
                # Remove unwanted fields
                param.pop("fit", None)
                param.pop("stderr", None)
                param.pop("constant", None)
    
    # Remove equations field entirely
    doc.pop("equations", None)
    doc.pop("tissues", None)


    return doc


# -------------------- RF helpers (reference mapping) -----------------
def build_rf_index_for_ec(data: Dict[str, pd.DataFrame], ec: str) -> Dict[str, Dict[str, Any]]:
    idx_map: Dict[str, Dict[str, Any]] = {}
    if 'RF' not in data:
        return idx_map
    rf_df = data['RF'][data['RF']['REC_ID'] == ec]
    for _, row in rf_df.iterrows():
        raw = str(row['VALUE_RAW'])
        idxs = ANGLE_GROUP_RE.findall(raw)
        ref_idxs: List[str] = []
        for g in idxs:
            ref_idxs.extend([x.strip() for x in g.split(",") if x.strip()])
        if not ref_idxs:
            continue
        pubmed_id = None
        m = PUBMED_RE.search(raw)
        if m:
            pubmed_id = m.group(1)
        clean = re.sub(r"<[^>]*>", "", raw)
        clean = re.sub(r"\{[^}]+\}", "", clean)
        parts = clean.split(':')
        title = None; authors = None; journal = "BRENDA Database"; year = None
        if len(parts) >= 2:
            authors = parts[0].strip()
            title = parts[1].strip()
            if len(parts) > 2:
                journal_part = ':'.join(parts[2:])
                year_m = re.search(r"\((\d{4})\)", journal_part)
                if year_m:
                    try:
                        year = int(year_m.group(1))
                    except Exception:
                        year = None
                journal = re.sub(r"\(\d{4}\)", "", journal_part).strip() or journal
        meta = {"title": title, "authors": authors, "journal": journal,
                "year": year, "pubmed_id": pubmed_id, "uri": None, "doi": None}
        for idx in ref_idxs:
            idx_map[idx] = meta
    return idx_map


# -------------------- Collect pair-level data ------------------
@dataclass
class CollectCtx:
    data: Dict[str, pd.DataFrame]
    ec: str
    protein_idx: str
    substrate: str


def collect_aux_for_pair(data: Dict[str, pd.DataFrame], ec: str, protein_idx: str, substrate: str) -> AuxData:
    aux = AuxData()

    def has_leading_pidx(raw: str) -> bool:
        return protein_idx in extract_leading_pidx_list(raw)

    def keep_seg_for_this(seg: Dict[str, Any]) -> bool:
        return protein_idx in (seg.get("pidxs") or [])

    # PR -> organism, uniprot, protein_notes
    if 'PR' in data:
        pr_df = data['PR'][data['PR']['REC_ID'] == ec]
        for _, row in pr_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            
            # Extract UniProt ID from braces
            m = UNIPROT_IN_BRACES.search(raw)
            if m:
                aux.uniprot_id = m.group(1)
            
            # Extract organism name (text before first parenthesis, after protein_idx)
            clean = PROTIDX_RE.sub("", raw, count=1)
            clean = re.sub(r"\s*\{[^}]*\}", "", clean)  # Remove UniProt braces
            clean = ANGLE_GROUP_RE.sub("", clean)  # Remove reference markers
            
            # Split at first parenthesis to get just the organism name
            organism_name = clean.split('(')[0].strip()
            if organism_name:
                aux.organism = organism_name
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            pr_notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        pr_notes.append(seg["note"])
            
            # If no segments with pidx markers but there are parentheses, extract general notes
            if not pr_notes:
                paren_content = extract_parenthetical_notes(raw)
                if paren_content:
                    pr_notes.append(paren_content)
            
            # Add PR notes to protein_notes (will be combined with MW notes later)
            if pr_notes:
                note_text = "; ".join(pr_notes)
                if aux.protein_notes:
                    aux.protein_notes = note_text + "; " + aux.protein_notes
                else:
                    aux.protein_notes = note_text
            
            break

    # SY -> synonyms
    if 'SY' in data:
        sy_df = data['SY'][data['SY']['REC_ID'] == ec]
        for _, row in sy_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            clean = PROTIDX_RE.sub("", raw, count=1)
            clean = ANGLE_GROUP_RE.sub("", clean).strip()
            if clean:
                aux.synonyms.append(clean)

    # RN / SN -> descriptive names (not used as protein.name)
    if 'RN' in data:
        rn_df = data['RN'][data['RN']['REC_ID'] == ec]
        if not rn_df.empty:
            txt = re.sub(r"^[\d.]+\s+RN\s+", "", str(rn_df.iloc[0]['VALUE_RAW'])).strip()
            if txt:
                aux.recommended_name = txt
    if 'SN' in data:
        sn_df = data['SN'][data['SN']['REC_ID'] == ec]
        if not sn_df.empty:
            txt = re.sub(r"^[\d.]+\s+SN\s+", "", str(sn_df.iloc[0]['VALUE_RAW'])).strip()
            if txt:
                aux.systematic_name = txt

    # KM
    if 'KM' in data:
        km_df = data['KM'][data['KM']['REC_ID'] == ec]
        for _, row in km_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            sub = extract_substrate(raw)
            if not sub or sub != substrate:
                continue
            val, unit = extract_value_and_unit(raw)
            segs = parse_per_ref_segments(raw)
            if segs:
                for seg in segs:
                    if not keep_seg_for_this(seg):
                        continue
                    aux.km_list.append({
                        "value": val, "unit": unit or "mM",
                        "substrate": sub, "pH": seg.get("pH"), "temperature": seg.get("temperature"),
                        "conditions": seg.get("note"), "ref_idx": seg.get("ref_idx")
                    })
                    if seg.get("ref_idx"):
                        aux.used_refs.add(seg["ref_idx"])
            else:
                refs = extract_trailing_refs_outside_parens(raw) or [None]
                pH, tC = parse_ph_temp(ALL_PARENS_RE.sub("", raw))
                note_txt = extract_parenthetical_notes(raw)
                for ref in refs:
                    aux.km_list.append({
                        "value": val, "unit": unit or "mM",
                        "substrate": sub, "pH": pH, "temperature": tC,
                        "conditions": note_txt, "ref_idx": ref
                    })
                    if ref:
                        aux.used_refs.add(ref)

    # KKM
    if 'KKM' in data:
        kkm_df = data['KKM'][data['KKM']['REC_ID'] == ec]
        for _, row in kkm_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            sub = extract_substrate(raw)
            if not sub or sub != substrate:
                continue
            val, unit = extract_value_and_unit(raw)
            segs = parse_per_ref_segments(raw)
            if segs:
                for seg in segs:
                    if not keep_seg_for_this(seg):
                        continue
                    note = seg.get("note")
                    extra = []
                    if seg.get("pH") is not None:
                        extra.append(f"pH {seg['pH']}")
                    if seg.get("temperature") is not None:
                        extra.append(f"{seg['temperature']}°C")
                    if extra:
                        note = (note + "; " if note else "") + ", ".join(extra)
                    aux.kkm_list.append({
                        "value": val, "unit": unit,
                        "substrate": sub, "pH": seg.get("pH"), "temperature": seg.get("temperature"),
                        "conditions": note, "ref_idx": seg.get("ref_idx")
                    })
                    if seg.get("ref_idx"):
                        aux.used_refs.add(seg["ref_idx"])
            else:
                refs = extract_trailing_refs_outside_parens(raw) or [None]
                pH, tC = parse_ph_temp(ALL_PARENS_RE.sub("", raw))
                note_txt = extract_parenthetical_notes(raw)
                extra = []
                if pH is not None:
                    extra.append(f"pH {pH}")
                if tC is not None:
                    extra.append(f"{tC}°C")
                if extra:
                    note_txt = (note_txt + "; " if note_txt else "") + ", ".join(extra)
                for ref in refs:
                    aux.kkm_list.append({
                        "value": val, "unit": unit,
                        "substrate": sub, "pH": pH, "temperature": tC,
                        "conditions": note_txt, "ref_idx": ref
                    })
                    if ref:
                        aux.used_refs.add(ref)

    # KI
    if 'KI' in data:
        ki_df = data['KI'][data['KI']['REC_ID'] == ec]
        for _, row in ki_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            sub = extract_substrate(raw)
            val, unit = extract_value_and_unit(raw)
            segs = parse_per_ref_segments(raw)
            if segs:
                for seg in segs:
                    if not keep_seg_for_this(seg):
                        continue
                    aux.ki_list.append({
                        "value": val, "unit": unit or "mM",
                        "substrate": sub, "pH": seg.get("pH"), "temperature": seg.get("temperature"),
                        "conditions": seg.get("note"), "ref_idx": seg.get("ref_idx")
                    })
                    if seg.get("ref_idx"):
                        aux.used_refs.add(seg["ref_idx"])
            else:
                refs = extract_trailing_refs_outside_parens(raw) or [None]
                pH, tC = parse_ph_temp(ALL_PARENS_RE.sub("", raw))
                note_txt = extract_parenthetical_notes(raw)
                for ref in refs:
                    aux.ki_list.append({
                        "value": val, "unit": unit or "mM",
                        "substrate": sub, "pH": pH, "temperature": tC,
                        "conditions": note_txt, "ref_idx": ref
                    })
                    if ref:
                        aux.used_refs.add(ref)

    # SA
    if 'SA' in data:
        sa_df = data['SA'][data['SA']['REC_ID'] == ec]
        for _, row in sa_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            val, unit = extract_value_and_unit(raw)
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            note = "; ".join(notes) if notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            refs = extract_trailing_refs_outside_parens(raw) or [None]
            for ref in refs:
                aux.specific_activities.append({"value": val, "unit": unit, "notes": note, "ref_idx": ref})
                if ref:
                    aux.used_refs.add(ref)

    # pH/Temperature related tags
    def one_note_ref(tag: str, attr_target: str):
        if tag not in data:
            return
        df = data[tag][data[tag]['REC_ID'] == ec]
        for _, row in df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            note = "; ".join(notes) if notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            refs = extract_trailing_refs_outside_parens(raw)
            val = None
            if tag in ('PHO', 'PI', 'TO'):
                m = re.search(r"([0-9]+(?:\.[0-9]+)?)", PROTIDX_RE.sub("", raw, count=1))
                if m:
                    try:
                        val = float(m.group(1))
                    except Exception:
                        val = None
            elif tag in ('PHR', 'TR'):
                m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*[-\u2013\u2014]\s*([0-9]+(?:\.[0-9]+)?)", raw)
                if m:
                    val = f"{m.group(1)}-{m.group(2)}"
            if attr_target == "ph_optimum" and aux.ph_optimum is None and val is not None:
                aux.ph_optimum = val
            if attr_target == "ph_range" and aux.ph_range is None and val:
                aux.ph_range = val
            if attr_target == "temp_optimum" and aux.temp_optimum is None and val is not None:
                aux.temp_optimum = val
            if attr_target == "temp_range" and aux.temp_range is None and val:
                aux.temp_range = val
            if attr_target == "ph_stability" and aux.ph_stability is None and note:
                aux.ph_stability = note
            if attr_target == "temp_stability" and aux.temp_stability is None and note:
                aux.temp_stability = note
            if attr_target == "storage_stability" and aux.storage_stability is None and note:
                aux.storage_stability = note
            if attr_target == "oxygen_stability" and aux.oxygen_stability is None and note:
                aux.oxygen_stability = note
            for r in refs or []:
                aux.used_refs.add(r)

    one_note_ref('PHO', 'ph_optimum')
    one_note_ref('PHR', 'ph_range')
    one_note_ref('PHS', 'ph_stability')
    one_note_ref('TO', 'temp_optimum')
    one_note_ref('TR', 'temp_range')
    one_note_ref('TS', 'temp_stability')
    one_note_ref('SS', 'storage_stability')
    one_note_ref('OS', 'oxygen_stability')

    # CF
    if 'CF' in data:
        cf_df = data['CF'][data['CF']['REC_ID'] == ec]
        for _, row in cf_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if protein_idx not in extract_all_pidxs(raw):
                continue
            
            # Find the leading protein index pattern
            m_pidx = PROTIDX_RE.search(raw)
            if not m_pidx:
                continue
            
            # Get text after the protein index, before any (, <, or ;
            start_pos = m_pidx.end()
            remainder = raw[start_pos:].strip()
            
            # Extract cofactor name: first token before (, <, ;, or {
            # Should be alphanumeric starting with a letter, no # or ,
            m_name = re.match(r"([A-Za-z][A-Za-z0-9\-+/]*)", remainder)
            if not m_name:
                continue
            name = m_name.group(1).strip()
            
            # Skip if name looks invalid
            if not name or '#' in name or ',' in name:
                continue
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
                    if seg.get("ref_idx"):
                        aux.used_refs.add(seg["ref_idx"])
            
            # Only use segment-specific notes, no generic fallback
            note_txt = "; ".join(notes) if notes else None
            
            if name:
                aux.cofactors.append({"name": name, "notes": note_txt})

    # IN
    if 'IN' in data:
        in_df = data['IN'][data['IN']['REC_ID'] == ec]
        for _, row in in_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if protein_idx not in extract_all_pidxs(raw):
                continue
            
            # Find leading protein index and extract name after it
            m_pidx = PROTIDX_RE.search(raw)
            if m_pidx:
                start_pos = m_pidx.end()
                remainder = raw[start_pos:].strip()
                m_name = re.match(r"([A-Za-z][A-Za-z0-9\-+/,.\s]*?)(?:\s*[(<;]|$)", remainder)
                nm = m_name.group(1).strip() if m_name else None
            else:
                base = PROTIDX_RE.sub("", raw, count=1).strip()
                m = re.match(r"([^\s(;<]+)", base)
                nm = m.group(1) if m else None
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
                    if seg.get("ref_idx"):
                        aux.used_refs.add(seg["ref_idx"])
            
            # Only use segment-specific notes
            note_txt = "; ".join(notes) if notes else None
            
            if nm:
                aux.inhibitors.append({"name": nm, "notes": note_txt})

    # AC (FIXED): only bare compound name in 'name' and id; parentheses go to notes
    if 'AC' in data:
        ac_df = data['AC'][data['AC']['REC_ID'] == ec]
        for _, row in ac_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            base = PROTIDX_RE.sub("", raw, count=1)
            base = ANGLE_GROUP_RE.sub("", base).strip()
            m = re.match(r"([^\s(;<]+)", base)
            nm = m.group(1) if m else None  # e.g., Isopropanol
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
            
            # Only use segment-specific notes
            note = "; ".join(notes) if notes else None
            
            if nm:
                aux.activators.append({"name": nm, "notes": note})

    # ME
    if 'ME' in data:
        me_df = data['ME'][data['ME']['REC_ID'] == ec]
        for _, row in me_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            nm = PROTIDX_RE.sub("", raw, count=1)
            nm = ANGLE_GROUP_RE.sub("", nm)
            nm = nm.strip().split(" ")[0]
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            # (this means the note applies to all leading pidxs)
            note = "; ".join(notes) if notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            refs = extract_trailing_refs_outside_parens(raw)
            for r in refs or []:
                aux.used_refs.add(r)
            if nm:
                aux.metals.append({"name": nm, "notes": note})

    # MW
    if 'MW' in data:
        mw_df = data['MW'][data['MW']['REC_ID'] == ec]
        notes = []
        for _, row in mw_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            clean = PROTIDX_RE.sub("", raw, count=1)
            mw_mult = MW_MULT_RE.search(clean)
            mws = None
            if mw_mult:
                mws = f"{mw_mult.group(1)} × {mw_mult.group(2)} Da"
            else:
                simple = re.search(r"(\d+)\s*(?:Da|kDa|)\b", clean)
                if simple:
                    try:
                        v = int(simple.group(1))
                        mws = f"{v} Da" if v > 1000 else f"{v*1000} Da"
                    except Exception:
                        pass
            if mws and not aux.protein_molecular_weight:
                aux.protein_molecular_weight = mws
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            seg_notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        seg_notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            note = "; ".join(seg_notes) if seg_notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            if note:
                notes.append(note)
            for r in extract_trailing_refs_outside_parens(raw) or []:
                aux.used_refs.add(r)
        if notes:
            mw_notes = "; ".join(notes)
            if aux.protein_notes:
                aux.protein_notes = aux.protein_notes + "; " + mw_notes
            else:
                aux.protein_notes = mw_notes

    # NSP ONLY - collect equations, deduplicate, merge notes (parentheses + brace tokens)
    if 'NSP' in data:
        df = data['NSP'][data['NSP']['REC_ID'] == ec]
        eq_map: Dict[str, Dict[str, Any]] = {}  # normalized_key -> {"equation": str, "notes_set": set}
        for _, row in df.iterrows():
            raw = str(row['VALUE_RAW'])
            if protein_idx not in extract_all_pidxs(raw):
                continue
            brace_tokens = [t.strip() for t in re.findall(r"\{([^{}]+)\}", raw)]
            eq_txt = re.sub(r"^[\d.]+\s+NSP\s+", "", raw).strip()
            eq_txt = PROTIDX_RE.sub("", eq_txt).strip()
            eq_txt = re.split(r"\s*\(|<|\{", eq_txt)[0].strip()
            if not eq_txt or re.fullmatch(r"#\s*\d+(?:,\d+)*\s*#", eq_txt):
                continue
            if not re.search(r"(=|⇌|→|←)", eq_txt):
                continue
            key = re.sub(r"\s+", " ", eq_txt).strip().lower()
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            paren_notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        paren_notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            paren_note = "; ".join(paren_notes) if paren_notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            note_parts = []
            if paren_note:
                note_parts.append(paren_note)
            note_parts.extend([t for t in brace_tokens if t])
            if key not in eq_map:
                eq_map[key] = {"equation": re.sub(r"\s+", " ", eq_txt).strip(),
                               "notes_set": set(note_parts)}
            else:
                eq_map[key]["notes_set"].update(note_parts)
        for v in eq_map.values():
            note = "; ".join(sorted(x for x in v["notes_set"] if x)) if v["notes_set"] else None
            aux.equations.append({"equation": v["equation"], "notes": note})

    # EN mutation
    if 'EN' in data:
        en_df = data['EN'][data['EN']['REC_ID'] == ec]
        for _, row in en_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            clean = PROTIDX_RE.sub("", raw, count=1).strip()
            m = re.match(r"([A-Z]\d+[A-Z](?:/[A-Z]\d+[A-Z])*)", clean)
            mutation = m.group(1) if m else None
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            note = "; ".join(notes) if notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            refs = extract_trailing_refs_outside_parens(raw)
            for r in refs or []:
                aux.used_refs.add(r)
            if mutation:
                aux.mutations.append({"mutation": mutation, "notes": note,
                                      "brenda_reference_index": (refs or [None])[0]})

    # PM post-translational
    if 'PM' in data:
        pm_df = data['PM'][data['PM']['REC_ID'] == ec]
        for _, row in pm_df.iterrows():
            raw = str(row['VALUE_RAW'])
            if not has_leading_pidx(raw):
                continue
            txt = PROTIDX_RE.sub("", raw, count=1).strip()
            txt = re.split(r"\s*<", txt)[0].strip()
            
            # Extract notes ONLY for segments matching this protein_idx
            segs = parse_per_ref_segments(raw)
            notes = []
            for seg in segs:
                if protein_idx in (seg.get("pidxs") or []):
                    if seg.get("note"):
                        notes.append(seg["note"])
            
            # If no segments with pidx markers, use all parenthetical content
            note = "; ".join(notes) if notes else (extract_parenthetical_notes(raw) if not segs else None)
            
            refs = extract_trailing_refs_outside_parens(raw)
            for r in refs or []:
                aux.used_refs.add(r)
            aux.post_trans_mods.append({"type": txt, "notes": note,
                                        "brenda_reference_index": (refs or [None])[0]})

    # ST tissues (FIXED: now properly filters by protein_idx and extracts segment-specific notes)
    if 'ST' in data:
        st_df = data['ST'][data['ST']['REC_ID'] == ec]
        for _, row in st_df.iterrows():
            raw = str(row['VALUE_RAW'])
            
            # Check if this protein_idx appears anywhere in the raw data
            if protein_idx not in extract_all_pidxs(raw):
                continue
            
            # Find the leading protein index and extract tissue name after it
            m_pidx = PROTIDX_RE.search(raw)
            if m_pidx:
                start_pos = m_pidx.end()
                remainder = raw[start_pos:].strip()
                # Extract tissue name: text before (, <, or ;
                m_name = re.match(r"([A-Za-z][A-Za-z0-9\s,.-]*?)(?:\s*[(<;]|$)", remainder)
                tissue_name = m_name.group(1).strip() if m_name else None
            else:
                # Fallback if no leading pidx marker
                base = PROTIDX_RE.sub("", raw, count=1).strip()
                base = ANGLE_GROUP_RE.sub("", base).strip()
                tissue_name = base.split('(')[0].strip() if base else None
            
            if tissue_name:
                aux.tissues.append(tissue_name)

    return aux


# --------------------- Build EnzymeML document -----------------
def build_doc_for_group(ec: str, protein_idx: str, substrate: str,
                        aux: AuxData,
                        tn_params: List[Dict],
                        data: Dict[str, pd.DataFrame]) -> Tuple[Dict[str, Any], set]:
    # protein.name strictly from PR (organism); otherwise fallback
    protein_name = aux.organism or f"protein#{protein_idx}"
    title = f"{ec} {protein_name} {substrate}"

    doc: Dict[str, Any] = {
        "name": title,
        "version": "2.0",
        "description": f"EnzymeML for EC {ec}, protein #{protein_idx}, substrate {substrate}",
        "created": datetime.now().isoformat(),
        "modified": None,
        "creators": [{"given_name": "BRENDA", "family_name": "Database", "mail": "info@brenda-enzymes.org"}]
    }

    prot_id = f"protein_{ec.replace('.', '_')}_p{protein_idx}"
    
    # Determine variant_type and variant_description
    if aux.mutations:
        # If there are mutations, this is a mutant
        variant_type = "mutant"
        # Combine all mutations into description
        mutation_strs = [m.get("mutation") for m in aux.mutations if m.get("mutation")]
        if mutation_strs:
            variant_description = "; ".join(mutation_strs)
        else:
            variant_description = "mutant"
    else:
        # No mutations means wildtype
        variant_type = "wildtype"
        variant_description = "wildtype"
    
    protein = {
        "id": prot_id,
        "name": protein_name,
        "ecnumber": ec,
        "uniprotid": aux.uniprot_id,
        "organism": aux.organism,
        "variant_type": variant_type,
        "variant_description": variant_description,
        "sequence": None,
        "organism_tax_id": None,
        "synonymous_names": aux.synonyms or None
    }
    # Don't add mutations list - it's now in variant_description
    if aux.post_translational_modifications if hasattr(aux, 'post_translational_modifications') else False:
        protein["post_translational_modifications"] = aux.post_translational_modifications
    elif aux.post_trans_mods:
        protein["post_translational_modifications"] = aux.post_trans_mods
    doc["proteins"] = [protein]

    # Small molecules
    small_molecules = []
    species_map = {}
    sid = f"substrate_{sanitize_id(substrate)}"
    small_molecules.append({
        "id": sid, "name": substrate,
        "synonymous_names": None,
        "inchi": None, "inchikey": None, "canonical_smiles": None,
        "CHEbIid": None, "KEGGid": None,
        "notes": None
    })
    species_map[substrate] = sid

    def ensure_mol(name: Optional[str], prefix: str, notes: Optional[str] = None):
        if not name:
            return None
        if name not in species_map:
            mid = f"{prefix}{sanitize_id(name)}"
            small_molecules.append({
                "id": mid, "name": name,
                "synonymous_names": None,
                "inchi": None, "inchikey": None, "canonical_smiles": None,
                "CHEbIid": None, "KEGGid": None,
                "notes": notes
            })
            species_map[name] = mid
        return species_map[name]

    # Note: ensure_mol will produce ids like 'activator_Isopropanol'
    for it in aux.inhibitors:
        ensure_mol(it.get("name"), "inhibitor_", it.get("notes"))
    for it in aux.activators:  # AC fixed: 'name' is bare compound name
        ensure_mol(it.get("name"), "activator_", it.get("notes"))
    for it in aux.cofactors:
        ensure_mol(it.get("name"), "cofactor_", it.get("notes"))
    for it in aux.metals:
        ensure_mol(it.get("name"), "metal_", it.get("notes"))

    doc["small_molecules"] = small_molecules if small_molecules else None

    # Reaction (with equation from NSP merged in)
    equation_text = None
    equation_notes = None
    smiles_equation = None
    if aux.equations:
        equation_text = aux.equations[0].get("equation")
        equation_notes = aux.equations[0].get("notes")
    
    if not equation_text:
        equation_text = f"substrate → product"
    
    organism_part = sanitize_id(aux.organism) if aux.organism else f"p{protein_idx}"
    reaction_id = f"reaction_{ec.replace('.', '_')}_{organism_part}_{sanitize_id(substrate)}"
    
    reaction = {
        "id": reaction_id,
        "name": f"Reaction catalyzed by {protein_name} on {substrate}",
        "equation": equation_text,
        "smiles_equation": smiles_equation,
        "reversible": True,
        "reactants": None, 
        "products": None,
        "modifiers": [{"species_id": prot_id, "role": "BIOCATALYST"}],
        "selectivity": None, 
        "Rheaid": None, 
        "KEGGid": None
    }
    
    # Add inhibitors/activators to reaction modifiers
    for inh in aux.inhibitors[:8]:
        nm = inh.get("name")
        if nm and nm in species_map:
            reaction["modifiers"].append({"species_id": species_map[nm], "role": "INHIBITOR"})
    for act in aux.activators[:8]:
        nm = act.get("name")
        if nm and nm in species_map:
            reaction["modifiers"].append({"species_id": species_map[nm], "role": "ACTIVATOR"})
    doc["reactions"] = [reaction]

    # ---------------- IDs with conditional uniprotid ----------------
    def _compose_id_parts(base_parts: List[str]) -> str:
        """Join non-empty parts with '_' after sanitizing."""
        sanitized = [sanitize_id(p) for p in base_parts if p is not None and str(p) != ""]
        return "_".join(sanitized)

    # Parameters
    parameters: List[Dict[str, Any]] = []
    param_counters: Dict[str, int] = {}  # Counter for each parameter type
    rf_index = build_rf_index_for_ec(data, ec)  # Get reference index mapping

    def add_param(name: str, value: Optional[float], unit: Optional[str],
                  assoc: Optional[str], pH: Optional[float], tC: Optional[float],
                  notes: Optional[str], ref_idx: Optional[str], ptype: str,
                  lower_bound: Optional[float] = None, upper_bound: Optional[float] = None):
        # Generate unique parameter ID
        param_key = f"{name}"
        if param_key not in param_counters:
            param_counters[param_key] = 0
        param_counters[param_key] += 1
        
        organism_part_local = sanitize_id(aux.organism) if aux.organism else f"p{protein_idx}"
        # Build parts conditionally including uniprotid if exists
        parts = [
            ec.replace('.', '_'),
            organism_part_local,
        ]
        if aux.uniprot_id:
            parts.append(aux.uniprot_id)
        parts.extend([
            substrate,
            name,
            str(param_counters[param_key])
        ])
        param_id = _compose_id_parts(parts)
        
        # Convert ref_idx to pubmed_id using RF index
        pubmed_id = None
        if ref_idx and ref_idx in rf_index:
            pubmed_id = rf_index[ref_idx].get("pubmed_id")
        
        parameters.append({
            "id": param_id,
            "name": name,
            "symbol": None,
            "value": None,
            "unit": unit,
            "initial_value": value,
            "upper_bound": upper_bound,
            "lower_bound": lower_bound,
            "associated_species": assoc,
            "original_name": None,
            "parameter.type": ptype,
            "pH": pH,
            "temperature": tC,
            "temperature_unit": "°C" if tC is not None else None,
            "notes": notes,
            "pubmed_id": pubmed_id
        })

    # TN -> kcat (ALWAYS in 1/s)
    for rec in tn_params:
        add_param("kcat", rec.get("value"),
                  "1/s",  # enforce 1/s
                  species_map.get(substrate), rec.get("pH"), rec.get("temperature"),
                  rec.get("note"), rec.get("ref_idx"), "KINETIC_CONSTANT")

    # KM
    for km in aux.km_list:
        add_param("Km", km.get("value"), km.get("unit") or "mM",
                  species_map.get(substrate), km.get("pH"), km.get("temperature"),
                  km.get("conditions"), km.get("ref_idx"), "KINETIC_CONSTANT")

    # kcat/Km
    for kkm in aux.kkm_list:
        add_param("kcat/Km", kkm.get("value"), kkm.get("unit"),
                  species_map.get(substrate), kkm.get("pH"), kkm.get("temperature"),
                  kkm.get("conditions"), kkm.get("ref_idx"), "SPECIFICITY_CONSTANT")

    # KI
    for ki in aux.ki_list:
        add_param("Ki", ki.get("value"), ki.get("unit") or "mM",
                  species_map.get(substrate) if ki.get("substrate") else None,
                  ki.get("pH"), ki.get("temperature"),
                  ki.get("conditions"), ki.get("ref_idx"), "INHIBITION_CONSTANT")

    # SA
    for sa in aux.specific_activities:
        add_param("SpecificActivity", sa.get("value"), sa.get("unit"),
                  None, None, None, sa.get("notes"), sa.get("ref_idx"), "EXPERIMENTAL")

    # Other condition parameters
    if aux.ph_range:
        lower, upper = parse_range(aux.ph_range)
        notes = None if (lower is not None and upper is not None) else aux.ph_range
        add_param("pH_range", None, None, None, None, None, notes, None, "CONDITION",
                  lower_bound=lower, upper_bound=upper)
    
    if aux.ph_stability:
        add_param("pH_stability", None, None, None, None, None, aux.ph_stability, None, "CONDITION")
    
    if aux.temp_range:
        lower, upper = parse_range(aux.temp_range)
        notes = None if (lower is not None and upper is not None) else aux.temp_range
        add_param("Temperature_range", None, None, None, None, None, notes, None, "CONDITION",
                  lower_bound=lower, upper_bound=upper)
    
    if aux.temp_stability:
        add_param("Temperature_stability", None, None, None, None, None, aux.temp_stability, None, "CONDITION")
    
    if aux.storage_stability:
        add_param("Storage_stability", None, None, None, None, None, aux.storage_stability, None, "CONDITION")
    
    if aux.oxygen_stability:
        add_param("Oxygen_stability", None, None, None, None, None, aux.oxygen_stability, None, "CONDITION")

    doc["parameters"] = parameters if parameters else None

    # Measurements - use pH_optimum and temp_optimum if available, fallback to TN
    measurement_ph = aux.ph_optimum if aux.ph_optimum is not None else next((r.get("pH") for r in tn_params if r.get("pH") is not None), None)
    measurement_temp = aux.temp_optimum if aux.temp_optimum is not None else next((r.get("temperature") for r in tn_params if r.get("temperature") is not None), None)
    
    if measurement_ph is not None or measurement_temp is not None:
        organism_part_local = sanitize_id(aux.organism) if aux.organism else f"p{protein_idx}"
        # measurement id parts: include uniprotid only if present
        m_parts = [
            "measurement",
            ec.replace('.', '_'),
            organism_part_local,
        ]
        if aux.uniprot_id:
            m_parts.append(aux.uniprot_id)
        m_parts.append(substrate)
        measurement_id = _compose_id_parts(m_parts)
        
        doc["measurements"] = [{
            "id": measurement_id,
            "name": "Experimental conditions (TN)",
            "species_data": None,
            "ph_opt": measurement_ph, 
            "temperature_opt": measurement_temp,
            "temperature_unit": "°C" if measurement_temp is not None else None
        }]
    else:
        doc["measurements"] = None

    # Additional lists
    if aux.cofactors:
        doc["cofactors"] = aux.cofactors
    if aux.inhibitors:
        doc["inhibitors"] = aux.inhibitors
    if aux.activators:
        doc["activators"] = aux.activators
    if aux.metals:
        doc["metals"] = aux.metals
    if aux.tissues:
        doc["tissues"] = aux.tissues

    # Top-level references: collect from ALL parameters that were actually added
    all_pubmed_ids = set()
    for param in parameters:
        pm_id = param.get("pubmed_id")
        if pm_id:
            all_pubmed_ids.add(pm_id)
    
    refs = []
    rf_index = build_rf_index_for_ec(data, ec)
    for pm_id in sorted(list(all_pubmed_ids), key=lambda x: int(x) if x and str(x).isdigit() else 10**9):
        # Find the ref_idx that maps to this pubmed_id
        ref_meta = None
        for idx, meta in rf_index.items():
            if meta.get("pubmed_id") == pm_id:
                ref_meta = meta
                break
        
        if ref_meta:
            refs.append({
                "title": ref_meta["title"],
                "authors": ref_meta["authors"],
                "journal": ref_meta["journal"],
                "year": ref_meta["year"],
                "pubmed_id": pm_id,
                "uri": f"https://pubmed.ncbi.nlm.nih.gov/{pm_id}/",
                "doi": ref_meta.get("doi")
            })
        else:
            refs.append({
                "title": f"Reference for EC {ec}",
                "authors": None,
                "journal": "BRENDA Database",
                "year": None,
                "pubmed_id": pm_id,
                "uri": f"https://pubmed.ncbi.nlm.nih.gov/{pm_id}/",
                "doi": None
            })
    doc["references"] = refs or None

    return doc, all_pubmed_ids


# ------------------------------ main ---------------------------
def main():
    parser = argparse.ArgumentParser(description="BRENDA -> EnzymeML exporter (TN-only refs; NSP-only equations; kcat=1/s; AC id/name clean; ST protein-filtered)")
    parser.add_argument("--brenda-dir", type=Path, required=True, help="Path to BRENDA ZIP or directory with CSVs")
    parser.add_argument("--schema", type=Path, required=True, help="Path to enzymeml-v2-extended.yaml")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--format", choices=["yaml", "json", "both"], default="json")
    parser.add_argument("--group-by", choices=["pair", "pair-ref", "tn-row"], default="pair",
                        help="Output granularity")
    parser.add_argument("--limit", type=int, default=0, help="Limit number of groups (0=no limit)")
    args = parser.parse_args()

    if args.format in ["yaml", "both"] and not HAVE_YAML:
        print("Error: PyYAML not installed. Install: pip install pyyaml")
        args.format = "json"

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load schema
    try:
        classes = load_schema(args.schema)
    except Exception as e:
        print(f"Error: cannot load schema: {e}")
        return 1

    # Load BRENDA
    print("Loading BRENDA data...")
    if args.brenda_dir.is_file() and args.brenda_dir.suffix == ".zip":
        data = load_data_from_zip(args.brenda_dir)
    elif args.brenda_dir.is_dir():
        data = load_data_from_dir(args.brenda_dir)
    else:
        print(f"Error: {args.brenda_dir} is not a valid ZIP or directory")
        return 1
    if not data:
        print("Error: no data loaded")
        return 1

    # Build index from TN
    print("\nBuilding index from TN...")
    if 'TN' not in data:
        print("Error: TN table not found.")
        return 1

    pair_index: Dict[Tuple, List[Dict[str, Any]]] = {}
    tn_df = data['TN']
    total_records = 0
    row_no = 0
    for _, row in tn_df.iterrows():
        row_no += 1
        ec = str(row['REC_ID'])
        raw = str(row['VALUE_RAW'])
        pidxs = extract_leading_pidx_list(raw)
        if not pidxs:
            continue
        substrate = extract_substrate(raw)
        if not substrate:
            continue

        # TN (kcat) parse: always 1/s (convert from 1/min if present), ignore parentheses for unit detection
        val, unit_kcat = parse_kcat_value_unit_brenda(raw)

        segments = parse_per_ref_segments(raw)
        if segments:
            for seg in segments:
                for pidx in seg["pidxs"]:
                    if pidx not in pidxs:
                        continue
                    rec = {"ec": ec, "pidx": pidx, "substrate": substrate,
                           "value": val, "unit": unit_kcat,
                           "pH": seg.get("pH"), "temperature": seg.get("temperature"),
                           "note": seg.get("note"), "ref_idx": seg.get("ref_idx"),
                           "tn_row": row_no}
                    key = (ec, pidx, substrate) if args.group_by == "pair" else \
                          (ec, pidx, substrate, seg.get("ref_idx") if args.group_by == "pair-ref" else row_no)
                    pair_index.setdefault(key, []).append(rec)
                    total_records += 1
        else:
            trailing_refs = extract_trailing_refs_outside_parens(raw) or [None]
            pH, tC = parse_ph_temp(ALL_PARENS_RE.sub("", raw))
            note_txt = extract_parenthetical_notes(raw)
            for pidx in pidxs:
                for ref_idx in trailing_refs:
                    rec = {"ec": ec, "pidx": pidx, "substrate": substrate,
                           "value": val, "unit": unit_kcat,
                           "pH": pH, "temperature": tC,
                           "note": note_txt, "ref_idx": ref_idx,
                           "tn_row": row_no}
                    key = (ec, pidx, substrate) if args.group_by == "pair" else \
                          (ec, pidx, substrate, ref_idx if args.group_by == "pair-ref" else row_no)
                    pair_index.setdefault(key, []).append(rec)
                    total_records += 1

    print(f"Built {total_records} TN-derived parameter records into {len(pair_index)} groups ({args.group_by}).")

    # Optional limit
    items = list(pair_index.items())
    if args.limit > 0:
        items = items[:args.limit]
        print(f"Limiting to {len(items)} groups.")

    # Export
    print("\nConverting groups to EnzymeML...")
    ok = 0
    err = 0
    used_names = set()

    for key, recs in items:
        try:
            ec, pidx, substrate = key[0], key[1], key[2]
            aux = collect_aux_for_pair(data, ec, pidx, substrate)

            tn_params = [{"value": r["value"], "unit": r["unit"], "pH": r["pH"],
                          "temperature": r["temperature"], "note": r["note"], "ref_idx": r["ref_idx"]}
                         for r in recs]

            raw_doc, _ = build_doc_for_group(ec, pidx, substrate, aux, tn_params, data)

            # Schema completion (every attribute present; missing -> null)
            completed = complete_doc_to_schema(raw_doc, classes)
            
            # Remove unwanted fields
            completed = clean_unwanted_fields(completed)

            # File name (no parentheses in names; clean tokens)
            prot_name = aux.organism or "unknown_protein"
            parts = [f"brenda_ec{ec}", sanitize_id(prot_name)]
            if aux.organism and aux.organism.lower() not in prot_name.lower():
                parts.append("org-" + sanitize_id(aux.organism))
            elif aux.uniprot_id:
                parts.append("uniprot-" + sanitize_id(aux.uniprot_id))
            parts.append("s_" + sanitize_id(substrate))
            if args.group_by == "pair-ref":
                this_ref = next((r.get("ref_idx") for r in recs if r.get("ref_idx")), None)
                if this_ref:
                    parts.append(f"ref{this_ref}")
            elif args.group_by == "tn-row":
                parts.append(f"row{recs[0]['tn_row']}")
            fname_base = "_".join(parts)
            candidate = fname_base
            suffix = 1
            while candidate in used_names or \
                  (args.format in ["yaml", "both"] and (args.output_dir / f"{candidate}.yaml").exists()) or \
                  (args.format in ["json", "both"] and (args.output_dir / f"{candidate}.json").exists()):
                suffix += 1
                candidate = f"{fname_base}__{suffix}"
            used_names.add(candidate)
            fname = candidate

            if args.format in ["json", "both"]:
                with open(args.output_dir / f"{fname}.json", "w", encoding="utf-8") as f:
                    json.dump(completed, f, ensure_ascii=False, indent=2)
            if args.format in ["yaml", "both"]:
                with open(args.output_dir / f"{fname}.yaml", "w", encoding="utf-8") as f:
                    yaml.safe_dump(completed, f, allow_unicode=True, sort_keys=False, default_flow_style=False)

            ok += 1
        except Exception as e:
            print(f"  Error on {key}: {e}")
            err += 1

    print("\nDone.")
    print(f"  Success: {ok}")
    print(f"  Errors : {err}")
    print(f"  Outdir : {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
