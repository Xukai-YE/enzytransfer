#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhea Step1: Convert to EnzymeML Format
Modified to output EXACTLY the same format as SABIO-RK (BRENDA-compatible output)
Usage:
    python step1_join.py  --input ../../data/Rhea/   --schema ../../schemas/enzymeml-v2-extended.yaml  --output-dir ../../output/Rhea/   --format yaml   --limit 100
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple

import json
import yaml
import pandas as pd

# ---------------------------------------------------------------------
# Project path
# ---------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from enzymeml_utils import (
    SchemaLoader,
    TypeConverter,
    CommonBuilders
)


# ============================================================
# SCHEMA COMPLETION (EXACT COPY FROM SABIO-RK/BRENDA)
# ============================================================

def build_class_template(class_name: str, classes: dict, cache: dict) -> dict:
    """Build a template with all attributes set to None"""
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


def _range_is_class(meta: dict, classes: dict):
    """Check if range is a class"""
    rng = meta.get("range")
    if rng and rng in classes:
        return rng
    return None


def fill_instance_to_schema(class_name: str, instance, classes: dict, cache: dict):
    """Fill instance with all schema attributes"""
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


def complete_doc_to_schema(doc: dict, classes: dict) -> dict:
    """Complete document to match full schema (EXACT COPY FROM SABIO-RK/BRENDA)"""
    root = "EnzymeMLDocument"
    cache = {}
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


def clean_unwanted_fields(doc: dict) -> dict:
    """
    Remove unwanted fields to match SABIO-RK/BRENDA output exactly
    Specifically removes: fit, stderr, constant from parameters, constant from all objects,
    stoichiometry: 1.0 from reactions, and equations from document
    """
    # Clean parameters
    if doc.get("parameters") and isinstance(doc["parameters"], list):
        for param in doc["parameters"]:
            if isinstance(param, dict):
                # Remove unwanted fields
                param.pop("fit", None)
                param.pop("stderr", None)
                param.pop("constant", None)
    
    # Remove constant from proteins
    if doc.get("proteins") and isinstance(doc["proteins"], list):
        for protein in doc["proteins"]:
            if isinstance(protein, dict):
                protein.pop("constant", None)
    
    # Remove constant from small_molecules
    if doc.get("small_molecules") and isinstance(doc["small_molecules"], list):
        for sm in doc["small_molecules"]:
            if isinstance(sm, dict):
                sm.pop("constant", None)
    
    # Remove constant from complexes
    if doc.get("complexes") and isinstance(doc["complexes"], list):
        for cpx in doc["complexes"]:
            if isinstance(cpx, dict):
                cpx.pop("constant", None)
    
    # Clean reactions: remove stoichiometry if it's 1.0
    if doc.get("reactions") and isinstance(doc["reactions"], list):
        for reaction in doc["reactions"]:
            if isinstance(reaction, dict):
                # Clean reactants
                if reaction.get("reactants") and isinstance(reaction["reactants"], list):
                    for reactant in reaction["reactants"]:
                        if isinstance(reactant, dict) and reactant.get("stoichiometry") == 1.0:
                            reactant.pop("stoichiometry", None)
                
                # Clean products
                if reaction.get("products") and isinstance(reaction["products"], list):
                    for product in reaction["products"]:
                        if isinstance(product, dict) and product.get("stoichiometry") == 1.0:
                            product.pop("stoichiometry", None)
    
    # Remove equations field entirely
    doc.pop("equations", None)
    
    return doc


# =====================================================================
# Utilities
# =====================================================================

def extract_chebi_id(chebi_url: str) -> str:
    """Return ChEBI ID in 'CHEBI:12345' form."""
    if not chebi_url or pd.isna(chebi_url):
        return ""
    s = str(chebi_url).strip()
    if 'CHEBI_' in s:
        return f'CHEBI:{s.split("CHEBI_")[-1].strip()}'
    if s.startswith('CHEBI:'):
        return s
    if s.isdigit():
        return f'CHEBI:{s}'
    return ""


def extract_chebi_number(chebi_id: str) -> str:
    """Extract only the numeric part from ChEBI ID (e.g., 'CHEBI:31011' -> '31011')"""
    if not chebi_id or pd.isna(chebi_id):
        return None
    s = str(chebi_id).strip()
    # Remove 'CHEBI:' prefix
    if s.startswith('CHEBI:'):
        s = s[6:]  # Remove 'CHEBI:'
    # Extract only digits
    if s.isdigit():
        return s
    return None


def load_tsv(filepath: Path) -> pd.DataFrame:
    """Load TSV safely; return empty DataFrame on failure."""
    try:
        return pd.read_csv(filepath, sep='\t', comment='#', header=None, low_memory=False)
    except Exception:
        return pd.DataFrame()

# =====================================================================
# Data loader
# =====================================================================

def load_rhea_data(rhea_dir: Path, schema_path: Path = None) -> Dict:
    """Load Rhea files into dict."""
    print(f"\n[Loading Rhea data from: {rhea_dir}]")
    data = {
        'ec': load_tsv(rhea_dir / 'rhea2ec.tsv'),
        'uniprot': load_tsv(rhea_dir / 'rhea2uniprot_sprot.tsv'),
        'smiles': load_tsv(rhea_dir / 'rhea-reaction-smiles.tsv'),
        'go': load_tsv(rhea_dir / 'rhea2go.tsv'),
        'kegg': load_tsv(rhea_dir / 'rhea2kegg_reaction.tsv'),
        'chebi_names': {},
        'chebi_smiles': {},
        'participants': pd.DataFrame(),
        'schema_templates': {},
        'schema_classes': {}
    }

    # schema templates and classes
    if schema_path and schema_path.exists():
        print(f"  ✓ Using schema: {schema_path}")
        schema = SchemaLoader(schema_path)
        data['schema_templates'] = schema.templates
        # Load schema classes for completion
        with open(schema_path, "r", encoding="utf-8") as f:
            sch = yaml.safe_load(f)
        data['schema_classes'] = sch.get("classes", {})
        print(f"  ✓ Loaded schema with {len(data['schema_classes'])} classes")
    else:
        print(f"  ⚠️  No schema file provided")

    # SMILES lookup
    if not data['smiles'].empty:
        data['smiles_lookup'] = {}
        for _, row in data['smiles'].iterrows():
            rhea_id = str(row[0]).strip()
            smiles_eq = str(row[1]).strip()
            if rhea_id and smiles_eq and smiles_eq != 'nan':
                data['smiles_lookup'][rhea_id] = smiles_eq
        print(f"  ✓ SMILES: {len(data['smiles_lookup'])} reactions")
    else:
        data['smiles_lookup'] = {}

    # ChEBI names
    for fname in ['chebiId_name.tsv', 'chebild_name.tsv']:
        df = load_tsv(rhea_dir / fname)
        if not df.empty:
            for _, row in df.iterrows():
                cid = extract_chebi_id(row[0])
                name = TypeConverter.safe_str(row[1])
                if cid and name:
                    data['chebi_names'][cid] = name
            print(f"  ✓ ChEBI names: {len(data['chebi_names'])}")
            break

    # ChEBI smiles
    df = load_tsv(rhea_dir / 'rhea-chebi-smiles.tsv')
    if not df.empty:
        for _, row in df.iterrows():
            cid = extract_chebi_id(row[0])
            smiles = TypeConverter.safe_str(row[1])
            if cid and smiles:
                data['chebi_smiles'][cid] = smiles
        print(f"  ✓ ChEBI SMILES: {len(data['chebi_smiles'])}")

    # participants file
    for fname in ['test.csv', 'rhea_participants.csv', 'test.xlsx']:
        fpath = rhea_dir / fname
        if fpath.exists():
            try:
                df = pd.read_excel(fpath) if fname.endswith('xlsx') else pd.read_csv(fpath)
                df.columns = df.columns.str.strip()
                data['participants'] = df
                print(f"  ✓ Participants: {len(df)} rows")
                break
            except Exception:
                pass

    return data


def get_mapping(rhea_id: str, df: pd.DataFrame, col_idx: int = 1, multi: bool = False):
    """Get mapped value(s) from a mapping table."""
    if df.empty:
        return [] if multi else ""
    mask = df[0].astype(str).str.strip() == str(rhea_id)
    matches = df[mask]
    if multi:
        return [TypeConverter.safe_str(row[col_idx]) for _, row in matches.iterrows() if col_idx < len(row)]
    return TypeConverter.safe_str(matches.iloc[0][col_idx]) if not matches.empty and col_idx < len(df.columns) else ""


def extract_participants(rhea_id: str, data: Dict) -> Tuple[List[Dict], List[Dict], str]:
    """From participants table, get substrates, products, and equation."""
    df = data['participants']
    if df.empty:
        return [], [], ""

    cols = {name: next((c for c in df.columns if name in c.lower()), None)
            for name in ['accession', 'equation', 'chebi', 'name', 'side']}
    if not all(cols.values()):
        return [], [], ""

    mask = df[cols['accession']].astype(str).str.contains(str(rhea_id), regex=False, na=False)
    rows = df[mask]
    if rows.empty:
        return [], [], ""

    equation = TypeConverter.safe_str(rows.iloc[0][cols['equation']])
    substrates, products = [], []

    for _, row in rows.iterrows():
        cid = extract_chebi_id(TypeConverter.safe_str(row[cols['chebi']]))
        if not cid:
            continue
        name = data['chebi_names'].get(cid, TypeConverter.safe_str(row[cols['name']]))
        smiles = data['chebi_smiles'].get(cid, "")
        particle = {'chebi_id': cid, 'name': name, 'smiles': smiles}
        side = TypeConverter.safe_str(row[cols['side']])
        if '_L' in side:
            substrates.append(particle)
        elif '_R' in side:
            products.append(particle)

    return substrates, products, equation

# =====================================================================
# EnzymeML document builder
# =====================================================================

def build_enzymeml_doc(rhea_id: str, data: Dict) -> Dict:
    """Build EnzymeML document (schema-aligned)."""
    ec = get_mapping(rhea_id, data['ec'])
    uniprot_ids = [uid for uid in get_mapping(rhea_id, data['uniprot'], col_idx=3, multi=True) if uid]
    substrates, products, equation = extract_participants(rhea_id, data)
    smiles = data['smiles_lookup'].get(str(rhea_id), "")
    kegg = get_mapping(rhea_id, data['kegg'])
    go_terms = [g for g in get_mapping(rhea_id, data['go'], multi=True) if g]  # reserved

    templates = data.get('schema_templates', {})

    # document
    doc_template = templates.get('EnzymeMLDocument', {})
    doc = {k: None for k in doc_template.keys()}
    doc.update({
        "name": f"Rhea Reaction {rhea_id}",
        "version": "2.0",
        "description": f"Enzymatic reaction from Rhea database, ID {rhea_id}",
        "created": datetime.now().isoformat(),
        "modified": None,
        "creators": None,
        "references": None,
        "vessels": None,
        "proteins": [],
        "complexes": None,
        "small_molecules": [],
        "reactions": None,
        "measurements": None,
        "parameters": None,
        "equations": None,
        "buffer": None,
        "notes": f"Rhea reaction {rhea_id}. No LLM validation applied."
    })

    # proteins
    if uniprot_ids:
        protein_template = templates.get('Protein', {})
        for i, uid in enumerate(uniprot_ids, 1):
            protein = {k: None for k in protein_template.keys()}
            protein.update({
                "id": f"protein{i}_rhea{rhea_id}",
                "name": f"Enzyme for rhea {rhea_id}",
                "uniprotid": uid,
                "ecnumber": ec or None,
                "organism": None,
                "organism_tax_id": None,
                "sequence": None
            })
            doc["proteins"].append(protein)

    # small molecules and reactions
    if substrates or products:
        sm_template = templates.get('SmallMolecule', {})
        species_map = {}

        for i, s in enumerate(substrates, 1):
            sid = f"substrate{i}_rhea{rhea_id}"
            species_map[f"s{i}"] = sid
            sm = {k: None for k in sm_template.keys()}
            sm.update({
                "id": sid,
                "name": s['name'],
                "canonical_smiles": s['smiles'] or None,
                "CHEbIid": extract_chebi_number(s['chebi_id'])
            })
            doc["small_molecules"].append(sm)

        for i, p in enumerate(products, 1):
            pid = f"product{i}_rhea{rhea_id}"
            species_map[f"p{i}"] = pid
            sm = {k: None for k in sm_template.keys()}
            sm.update({
                "id": pid,
                "name": p['name'],
                "canonical_smiles": p['smiles'] or None,
                "CHEbIid": extract_chebi_number(p['chebi_id'])
            })
            doc["small_molecules"].append(sm)

        # reaction
        reaction_template = templates.get('Reaction', {})
        reaction = {k: None for k in reaction_template.keys()}
        try:
            rid = int(rhea_id)
            is_reversible = (rid % 2 == 0)
        except Exception:
            is_reversible = True

        reaction.update({
            "id": f"reaction_rhea{rhea_id}",
            "name": f"rhea_reaction_{rhea_id}",
            "reversible": is_reversible,
            "reactants": [CommonBuilders.build_reaction_element(species_map[f"s{i}"], 1.0)
                          for i in range(1, len(substrates) + 1)],
            "products": [CommonBuilders.build_reaction_element(species_map[f"p{i}"], 1.0)
                         for i in range(1, len(products) + 1)],
            "modifiers": [CommonBuilders.build_modifier_element(p["id"], "BIOCATALYST")
                          for p in doc["proteins"]],
            "Rheaid": str(rhea_id),
            "KEGGid": kegg or None
        })

        # kinetic law/equation
        equation_template = templates.get('Equation', {})
        if equation_template:
            kinetic_law = {k: None for k in equation_template.keys()}
            kinetic_law.update({
                "species_id": f"reaction_rhea{rhea_id}",
                "equation": equation or None,
                "equation_type": "RATE_LAW",
                "smiles_equation": smiles or None,
                "variables": None
            })
            reaction["kinetic_law"] = kinetic_law

        doc["reactions"] = [reaction]

        # complex
        complex_template = templates.get('Complex', {})
        complex_obj = {k: None for k in complex_template.keys()}
        complex_obj.update({"id": f"complex1_rhea{rhea_id}", "name": "reaction_complex"})
        doc["complexes"] = [complex_obj]

    return doc

# =====================================================================
# Local writer to force .yaml/.json outputs
# =====================================================================

def write_document_local(doc: Dict, out_dir: Path, filename_stem: str, fmt: str) -> None:
    """
    Explicitly write a single file per reaction with extension.
    fmt in {'yaml', 'json', 'both'}
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if fmt in ('yaml', 'both'):
        ypath = out_dir / f"{filename_stem}.yaml"
        with open(ypath, "w", encoding="utf-8") as f:
            yaml.safe_dump(doc, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
        print(f"  ✓ {ypath.name}")

    if fmt in ('json', 'both'):
        jpath = out_dir / f"{filename_stem}.json"
        with open(jpath, "w", encoding="utf-8") as f:
            json.dump(doc, f, indent=2, ensure_ascii=False)
        print(f"  ✓ {jpath.name}")

# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="Convert Rhea data to EnzymeML V2 format (SABIO-RK/BRENDA-compatible output)")
    parser.add_argument("--input", required=True, help="Input directory containing Rhea TSV files (or file path, will use parent dir)")
    parser.add_argument("--schema", required=True, help="EnzymeML v2 extended YAML schema file")
    parser.add_argument("--output-dir", required=True, help="Output directory for EnzymeML files")
    parser.add_argument("--format", choices=['yaml', 'json', 'both'], default='yaml', help='Output format: yaml, json, or both (default: yaml)')
    parser.add_argument("--rhea-ids", "-r", nargs='+', help="Specific Rhea IDs to process")
    parser.add_argument("--limit", "-l", type=int, default=20000, help="Maximum number of reactions to process (default: 10)")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"❌ Error: Input path not found: {input_path}")
        return 1

    rhea_dir = input_path.parent if input_path.is_file() else input_path

    schema_path = Path(args.schema)
    if not schema_path.exists():
        print(f"❌ Error: Schema file not found: {schema_path}")
        return 1

    data = load_rhea_data(rhea_dir, schema_path)

    if args.rhea_ids:
        rhea_ids = args.rhea_ids
        print(f"Processing specified Rhea IDs: {rhea_ids}")
    elif not data['participants'].empty:
        acc_col = [c for c in data['participants'].columns if 'accession' in c.lower()]
        if acc_col:
            rhea_ids = data['participants'][acc_col[0]].astype(str).str.extract(r'(\d+)')[0].unique().tolist()
            rhea_ids = [r for r in rhea_ids if r and r != 'nan'][:args.limit]
            print(f"Found {len(rhea_ids)} unique Rhea IDs from participants file")
        else:
            print("❌ Error: Could not find 'accession' column in participants file")
            return 1
    else:
        print("❌ Error: No Rhea IDs found. Provide --rhea-ids or ensure participants file exists")
        return 1

    print(f"\n{'='*60}")
    print(f"Processing {len(rhea_ids)} Rhea reactions")
    print(f"Output format: {args.format} (SABIO-RK/BRENDA-compatible)")
    print(f"{'='*60}\n")

    out_dir = Path(args.output_dir)
    success_count = 0
    error_count = 0

    schema_classes = data.get('schema_classes', {})

    for rid in rhea_ids:
        try:
            # Build basic document
            doc = build_enzymeml_doc(rid, data)
            
            # Complete to full schema (fills all missing fields with null)
            if schema_classes:
                doc = complete_doc_to_schema(doc, schema_classes)
            
            # Clean unwanted fields to match SABIO-RK/BRENDA exactly
            doc = clean_unwanted_fields(doc)
            
            filename = f"rhea_{rid}"
            write_document_local(doc, out_dir, filename, args.format)
            success_count += 1
        except Exception as e:
            print(f"  ⚠️  Error processing Rhea {rid}: {e}")
            import traceback
            traceback.print_exc()
            error_count += 1
            continue

    print(f"\n{'='*60}")
    print(f"✅ Done!")
    print(f"   Success: {success_count}/{len(rhea_ids)} files")
    if error_count > 0:
        print(f"   Errors:  {error_count}")
    print(f"   Output:  {out_dir}/")
    print(f"   Format:  EXACT SABIO-RK/BRENDA match (UNIFIED)")
    print(f"{'='*60}\n")
    return 0


if __name__ == "__main__":
    exit(main())
