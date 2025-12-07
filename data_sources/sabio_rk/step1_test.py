#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SABIO-RK Step1: Convert to EnzymeML Format (BRENDA-compatible output)
======================================================================
Modified to output EXACTLY the same format as BRENDA step1_join.py
UNIFIED VERSION - Output format identical to BRENDA

python step1_extract.py --input ../../data/Sabio_rk/final_standardized.csv --schema ../../schemas/enzymeml-v2-extended.yaml --output-dir ../../output/sabio_rk --format yaml
"""

import argparse
import re
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from enzymeml_utils import (
    SchemaLoader,
    TypeConverter,
    FilenameUtils,
)

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False
    print("Warning: PyYAML not installed. Install: pip install pyyaml")


# ============================================================
# SCHEMA COMPLETION (EXACT COPY FROM BRENDA)
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
    """Complete document to match full schema (EXACT COPY FROM BRENDA)"""
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
    Remove unwanted fields to match BRENDA output exactly
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
    
    return doc


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def sanitize_id(s: str) -> str:
    """Sanitize string for use as ID"""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def clean_unit(unit_str) -> str:
    """
    Clean unit string - fix common encoding issues
    
    Examples:
        s?1 -> 1/s
        M?1s?1 -> 1/(M*s)
        M -> M
    """
    if not unit_str or pd.isna(unit_str):
        return None
    
    unit = str(unit_str).strip()
    
    # Fix common patterns
    # s?1 -> 1/s (or s-1)
    if unit == 's?1':
        return '1/s'
    # M?1s?1 -> 1/(M*s) or M-1s-1
    if unit == 'M?1s?1':
        return '1/(M*s)'
    # M?1 -> 1/M
    if unit == 'M?1':
        return '1/M'
    
    # General pattern: replace ?1 with -1
    unit = unit.replace('?1', '-1')
    unit = unit.replace('?', '')
    
    return unit if unit else None


def extract_enzyme_name(raw_name):
    """Extract enzyme name"""
    if not raw_name or pd.isna(raw_name):
        return None
    name = str(raw_name).strip()
    name = re.sub(r'\([^)]*mutant[^)]*\)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\([^)]*variant[^)]*\)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\([^)]*wild[- ]?type[^)]*\)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\s+(mutant|variant|wild[- ]?type).*$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\s+', ' ', name).strip()
    return name if name else None


def extract_enzyme_variant(raw_variant):
    """Extract enzyme variant description"""
    if not raw_variant or pd.isna(raw_variant):
        return None
    variant = str(raw_variant).strip()
    if variant.lower() in ['wild type', 'wildtype', 'wt', 'wild-type', 'n/a', '']:
        return None
    return variant


def extract_substrates(raw_substrate):
    """Extract list of substrates"""
    if not raw_substrate or pd.isna(raw_substrate):
        return None
    substrate_str = str(raw_substrate).strip()
    if not substrate_str or substrate_str.lower() in ['', 'n/a', 'nan', 'none']:
        return None
    substrates = []
    for item in re.split(r'[;,]', substrate_str):
        item = item.strip()
        if item:
            substrates.append(item)
    return substrates if substrates else None


def extract_ph(row) -> float:
    """Extract pH value from row"""
    for col in ['Measurement.ph', 'Measurement.pH', 'ph', 'pH']:
        if col in row.index:
            val = row[col]
            if val and pd.notna(val):
                try:
                    return float(val)
                except:
                    pass
    return None


def extract_temperature(row) -> float:
    """Extract temperature value from row"""
    for col in ['Measurement.temperature', 'temperature', 'Temperature']:
        if col in row.index:
            val = row[col]
            if val and pd.notna(val):
                try:
                    # Remove units like °C, C, celsius
                    temp_str = str(val).replace('°C', '').replace('C', '').replace('celsius', '').strip()
                    return float(temp_str)
                except:
                    pass
    return None


def extract_buffer(row) -> str:
    """Extract buffer composition from row"""
    for col in ['EnzymeMLDocument.buffer', 'Buffer.buffer', 'buffer', 'Buffer']:
        if col in row.index:
            val = row[col]
            if val and pd.notna(val):
                buffer_str = str(val).strip()
                if buffer_str and buffer_str.lower() not in ['', 'n/a', 'nan', 'none']:
                    return buffer_str
    return None


def extract_buffer_concentration(row) -> str:
    """Extract first buffer concentration from row"""
    for col in ['Buffer.concentration', 'buffer.concentration', 'buffer_concentration']:
        if col in row.index:
            val = row[col]
            if val and pd.notna(val):
                conc_str = str(val).strip()
                if conc_str and conc_str.lower() not in ['', 'n/a', 'nan', 'none']:
                    # If multiple concentrations separated by semicolon, take the first one
                    if ';' in conc_str:
                        first_conc = conc_str.split(';')[0].strip()
                        return first_conc
                    return conc_str
    return None


def parse_reaction_participants(raw_str) -> list:
    """
    Parse reactants or products string into structured list
    Examples:
        "1 Plasminogen" -> [{"species": "Plasminogen", "stoichiometry": 1.0}]
        "Glucose; ATP" -> [{"species": "Glucose", "stoichiometry": 1.0}, ...]
    """
    if not raw_str or pd.isna(raw_str):
        return None
    
    participants = []
    items = re.split(r'[;,]', str(raw_str))
    
    for item in items:
        item = item.strip()
        if not item or item.lower() in ['n/a', 'nan', 'none', '']:
            continue
        
        # Try to extract stoichiometry (e.g., "2 ATP" -> stoich=2, species="ATP")
        match = re.match(r'^(\d+(?:\.\d+)?)\s+(.+)$', item)
        if match:
            stoich = float(match.group(1))
            species = match.group(2).strip()
        else:
            stoich = 1.0
            species = item
        
        participants.append({
            "species": species,
            "stoichiometry": stoich
        })
    
    return participants if participants else None


# ============================================================
# DOCUMENT COMPOSER (EXACT BRENDA FORMAT - UNIFIED)
# ============================================================

class DocumentComposer:
    """Compose EnzymeML documents in EXACT BRENDA format"""
    
    def __init__(self, converter: TypeConverter):
        self.converter = converter
    
    def compose(
        self,
        row: pd.Series,
        enzyme_name: str,
        enzyme_variant: str,
        substrates: list,
        parameters: list,
        entry_id: int,
        ph: float = None,
        temperature: float = None,
        buffer_composition: str = None,
        buffer_concentration: str = None,
        reactants_parsed: list = None,
        products_parsed: list = None
    ) -> dict:
        """Compose document in EXACT BRENDA format"""
        
        # Extract and clean pubmed_id FIRST (for use in parameters)
        pubmed_cols = ['Reference.pubmed_id', 'PubMedID', 'pubmed_id']
        pubmed_id_raw = None
        for col in pubmed_cols:
            if col in row.index:
                pubmed_id_raw = row[col]
                if pubmed_id_raw and pd.notna(pubmed_id_raw):
                    break
        
        # Clean pubmed_id: remove .0 suffix
        pubmed_id_clean = None
        if pubmed_id_raw and pd.notna(pubmed_id_raw):
            pubmed_str = str(pubmed_id_raw)
            # Remove .0 suffix if present
            if pubmed_str.endswith('.0'):
                pubmed_str = pubmed_str[:-2]
            pubmed_id_clean = pubmed_str
        
        # Build document with EXACT BRENDA structure
        doc = {
            "name": f"SABIO-RK Entry {entry_id}",
            "version": "2.0",
            "description": f"EnzymeML for SABIO-RK entry {entry_id}",
            "created": datetime.now().isoformat(),
            "modified": None,
            "creators": [{
                "given_name": "SABIO-RK",
                "family_name": "Database",
                "mail": "sabiork@h-its.org"
            }]
        }
        
        # Extract EC number for protein
        ec_number = None
        for ec_col in ['Protein.ecnumber', 'ecnumber', 'EC', 'ec', 'Enzyme.ecnumber']:
            if ec_col in row.index:
                ec_val = row[ec_col]
                if ec_val and pd.notna(ec_val) and str(ec_val).strip() not in ['', 'n/a', 'nan', 'none']:
                    ec_number = str(ec_val).strip()
                    break
        
        # Extract UniProt ID
        uniprot_id = None
        for up_col in ['Protein.uniprotid', 'uniprotid', 'UniProt', 'uniprot']:
            if up_col in row.index:
                up_val = row[up_col]
                if up_val and pd.notna(up_val) and str(up_val).strip() not in ['', 'n/a', 'nan', 'none']:
                    uniprot_id = str(up_val).strip()
                    break
        
        # Extract organism
        organism = None
        for org_col in ['Protein.organism', 'organism', 'Organism']:
            if org_col in row.index:
                org_val = row[org_col]
                if org_val and pd.notna(org_val) and str(org_val).strip() not in ['', 'n/a', 'nan', 'none']:
                    organism = str(org_val).strip()
                    break
        
        # Proteins - EXACT BRENDA field order
        ec_safe = ec_number.replace('.', '_') if ec_number else "unknown"
        protein_idx = str(entry_id)
        
        # Build protein ID: include uniprotid if available (like BRENDA)
        if uniprot_id:
            prot_id = f"protein_{ec_safe}_{sanitize_id(uniprot_id)}_p{protein_idx}"
        else:
            prot_id = f"protein_{ec_safe}_p{protein_idx}"
        
        # Determine variant_type and variant_description like BRENDA
        variant_type = None
        variant_description = None
        if enzyme_variant:
            # If there's a variant description, it's a variant/mutant
            variant_type = "mutant"
            variant_description = enzyme_variant
        else:
            # No variant means wildtype
            variant_type = "wildtype"
            variant_description = "wildtype"
        
        protein = {
            "id": prot_id,
            "name": enzyme_name or f"protein#{entry_id}",
            "ecnumber": ec_number,
            "uniprotid": uniprot_id,
            "organism": organism,
            "variant_type": variant_type,
            "variant_description": variant_description,
            "sequence": None,
            "organism_tax_id": None,
            "synonymous_names": None
        }
        doc["proteins"] = [protein]
        
        # Small molecules - EXACT BRENDA field order
        small_molecules = []
        species_map = {}
        if substrates:
            for substrate in substrates:
                sid = f"substrate_{sanitize_id(substrate)}"
                small_molecules.append({
                    "id": sid,
                    "name": substrate,
                    "synonymous_names": None,
                    "inchi": None,
                    "inchikey": None,
                    "canonical_smiles": None,
                    "CHEbIid": None,
                    "KEGGid": None,
                    "notes": None
                })
                species_map[substrate] = sid
        
        # Add reactants and products to small_molecules
        if reactants_parsed:
            for r in reactants_parsed:
                species_name = r['species']
                if species_name not in species_map:
                    mol_id = f"reactant_{sanitize_id(species_name)}"
                    small_molecules.append({
                        "id": mol_id,
                        "name": species_name,
                        "synonymous_names": None,
                        "inchi": None,
                        "inchikey": None,
                        "canonical_smiles": None,
                        "CHEbIid": None,
                        "KEGGid": None,
                        "notes": None
                    })
                    species_map[species_name] = mol_id
        
        if products_parsed:
            for p in products_parsed:
                species_name = p['species']
                if species_name not in species_map:
                    mol_id = f"product_{sanitize_id(species_name)}"
                    small_molecules.append({
                        "id": mol_id,
                        "name": species_name,
                        "synonymous_names": None,
                        "inchi": None,
                        "inchikey": None,
                        "canonical_smiles": None,
                        "CHEbIid": None,
                        "KEGGid": None,
                        "notes": None
                    })
                    species_map[species_name] = mol_id
        
        doc["small_molecules"] = small_molecules if small_molecules else None
        
        # Reactions - EXACT BRENDA field order
        # organism_part: should include BOTH organism and uniprotid if available
        # Format: ec + organism + uniprotid + substrate
        organism_parts = []
        if organism:
            # Take first two words (genus species) from organism
            org_words = organism.split()[:2]
            organism_parts.append("_".join([sanitize_id(w) for w in org_words]))
        if uniprot_id:
            organism_parts.append(sanitize_id(uniprot_id))
        
        # Build organism_part: organism_uniprotid or fallback
        if organism_parts:
            organism_part = "_".join(organism_parts)
        else:
            organism_part = f"p{protein_idx}"
        
        substrate_safe = sanitize_id(substrates[0]) if substrates else "substrate"
        reaction_id = f"reaction_{ec_safe}_{organism_part}_{substrate_safe}"
        
        # Build reactants list
        reactants_list = None
        if reactants_parsed:
            reactants_list = []
            for r in reactants_parsed:
                reactant_obj = {"species_id": species_map[r['species']]}
                # Only include stoichiometry if it's not 1.0
                if r['stoichiometry'] != 1.0:
                    reactant_obj["stoichiometry"] = r['stoichiometry']
                reactants_list.append(reactant_obj)
        
        # Build products list
        products_list = None
        if products_parsed:
            products_list = []
            for p in products_parsed:
                product_obj = {"species_id": species_map[p['species']]}
                # Only include stoichiometry if it's not 1.0
                if p['stoichiometry'] != 1.0:
                    product_obj["stoichiometry"] = p['stoichiometry']
                products_list.append(product_obj)
        
        # Build equation text like BRENDA
        equation_text = None
        if reactants_parsed or products_parsed:
            # Try to build equation from reactants/products
            reactant_names = [r['species'] for r in (reactants_parsed or [])]
            product_names = [p['species'] for p in (products_parsed or [])]
            if reactant_names or product_names:
                left = " + ".join(reactant_names) if reactant_names else "substrate"
                right = " + ".join(product_names) if product_names else "product"
                equation_text = f"{left} → {right}"
        
        if not equation_text:
            equation_text = f"{substrates[0] if substrates else 'substrate'} → product"
        
        reaction = {
            "id": reaction_id,
            "name": f"Reaction catalyzed by {enzyme_name or 'enzyme'} on {substrates[0] if substrates else 'substrate'}",
            "equation": equation_text,
            "smiles_equation": None,
            "reversible": True,
            "reactants": reactants_list,
            "products": products_list,
            "modifiers": [{"species_id": prot_id, "role": "BIOCATALYST"}],
            "selectivity": None,
            "Rheaid": None,
            "KEGGid": None
        }
        doc["reactions"] = [reaction]
        
        # Parameters - EXACT BRENDA field order (WITHOUT fit, stderr, constant)
        param_list = []
        if parameters:
            param_counter = {}
            for param in parameters:
                ptype = param['type']
                param_name = param['name']
                param_type_enum = param['parameter_type']
                
                # Counter for unique IDs
                if ptype not in param_counter:
                    param_counter[ptype] = 0
                param_counter[ptype] += 1
                
                # Build ID like BRENDA: EC_organism_substrate_parameter_counter
                param_id = f"{ec_safe}_{organism_part}_{substrate_safe}_{sanitize_id(param_name)}_{param_counter[ptype]}"
                
                # EXACT BRENDA parameter structure (15 fields, NOT 18)
                # Removed: fit, stderr, constant
                param_obj = {
                    "id": param_id,
                    "name": param_name,
                    "symbol": None,
                    "value": None,
                    "unit": param.get('unit'),
                    "initial_value": self.converter.convert(param.get('value'), 'float'),
                    "upper_bound": None,
                    "lower_bound": None,
                    "associated_species": species_map.get(substrates[0]) if substrates else None,
                    "original_name": None,
                    "parameter.type": param_type_enum,
                    "pH": ph,
                    "temperature": temperature,
                    "temperature_unit": "°C" if temperature is not None else None,
                    "notes": None,
                    "pubmed_id": pubmed_id_clean  # Use cleaned pubmed_id from references
                }
                param_list.append(param_obj)
        
        doc["parameters"] = param_list if param_list else None
        
        # Measurements - EXACT BRENDA field order (use ph_opt and temperature_opt like BRENDA)
        if ph is not None or temperature is not None:
            measurement_id = f"measurement_{ec_safe}_{organism_part}_{substrate_safe}"
            doc["measurements"] = [{
                "id": measurement_id,
                "name": "Experimental conditions (TN)",
                "species_data": None,
                "ph_opt": ph,  # Like BRENDA: ph_opt instead of ph
                "temperature_opt": temperature,  # Like BRENDA: temperature_opt
                "temperature_unit": "°C" if temperature is not None else None
            }]
        else:
            doc["measurements"] = None
        
        # Buffer - wrap as BUFFER object per schema
        if buffer_composition or buffer_concentration:
            doc["buffer"] = {
                "buffer": buffer_composition,
                "buffer_concentration": buffer_concentration
            }
        else:
            doc["buffer"] = None
        
        # References - EXACT BRENDA structure (using pubmed_id_clean extracted at start)
        if pubmed_id_clean:
            doc["references"] = [{
                "title": f"Reference for SABIO-RK entry {entry_id}",
                "authors": None,
                "journal": "SABIO-RK Database",
                "year": None,
                "pubmed_id": pubmed_id_clean,
                "uri": f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id_clean}/",
                "doi": None
            }]
        else:
            doc["references"] = None
        
        return doc


# ============================================================
# YAML GENERATOR
# ============================================================

class YAMLGenerator:
    """Generate YAML/JSON files in EXACT BRENDA format"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
        self.converter = TypeConverter()
        self.composer = DocumentComposer(self.converter)
    
    def generate(self, df: pd.DataFrame, output_dir: str, format: str = 'yaml', limit: int = 0, group_by: str = 'entry'):
        """
        Generate files in EXACT BRENDA format
        
        Args:
            df: Input dataframe
            output_dir: Output directory
            format: Output format ('yaml', 'json', or 'both')
            limit: Limit number of entries (0=no limit)
            group_by: Grouping strategy:
                - 'entry': One file per EntryID (default, separate experiments)
                - 'protein_substrate': Merge entries with same protein+substrate
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"\nGenerating BRENDA-format files...")
        print(f"  Output format: {format}")
        print(f"  Grouping by: {group_by}")
        print(f"  All fields match BRENDA exactly (UNIFIED VERSION)")
        
        # Group data based on strategy
        if group_by == 'protein_substrate':
            # Merge entries with same EC, organism, UniProt, and substrate
            grouping_cols = []
            for col in ['Protein.ecnumber', 'Protein.organism', 'Protein.uniprotid', 'Reaction.reactants']:
                if col in df.columns:
                    grouping_cols.append(col)
            if grouping_cols:
                grouped = list(df.groupby(grouping_cols))
            else:
                print("Warning: Required columns for protein_substrate grouping not found, falling back to EntryID")
                grouped = list(df.groupby("EntryID"))
        else:
            # Default: one file per EntryID
            grouped = list(df.groupby("EntryID"))
        total_entries = len(grouped)
        
        # Apply limit if specified
        if limit > 0 and limit < total_entries:
            grouped = grouped[:limit]
            print(f"  Limiting to {limit} entries (out of {total_entries} total)")
        else:
            print(f"  Processing all {total_entries} entries")
        
        count = 0
        for entry_id, group in grouped:
            row = group.iloc[0]
            
            # Extract data
            enzyme_name = extract_enzyme_name(row.get('Protein.name'))
            enzyme_variant = extract_enzyme_variant(row.get('Protein.variant_description'))
            substrates = extract_substrates(row.get('Reaction.reactants'))
            
            # Extract pH, temperature, and buffer
            ph = extract_ph(row)
            temperature = extract_temperature(row)
            buffer_comp = extract_buffer(row)
            buffer_conc = extract_buffer_concentration(row)
            
            # Extract reactants and products
            reactants_raw = row.get('Reaction.reactants')
            products_raw = row.get('Reaction.products')
            reactants_parsed = parse_reaction_participants(reactants_raw)
            products_parsed = parse_reaction_participants(products_raw)
            
            # Extract parameters with correct type detection and standard units
            params = []
            for _, r in group.iterrows():
                pname = r.get('Parameter.name') or r.get('parameter.name') or ''
                if not pname or pd.isna(pname):
                    continue
                
                pname_str = str(pname).strip()
                pname_lower = pname_str.lower()
                pval = r.get('Parameter.initial_value') or r.get('parameter.startValue')
                
                if not pd.notna(pval):
                    continue
                
                # Determine parameter type and standard unit based on parameter name
                # Priority: kcat/Km > Km > kcat > Vmax
                
                if 'kcat/km' in pname_lower or 'kcatkm' in pname_lower or 'kcat_km' in pname_lower:
                    # Catalytic efficiency: kcat/Km
                    params.append({
                        'type': 'kcat_km',
                        'name': 'kcat/Km',
                        'value': pval,
                        'unit': '1/(M*s)',  # Standard unit for catalytic efficiency
                        'parameter_type': 'SPECIFICITY_CONSTANT'
                    })
                elif 'km' in pname_lower and 'kcat' not in pname_lower:
                    # Michaelis constant: Km
                    params.append({
                        'type': 'km',
                        'name': 'Km',
                        'value': pval,
                        'unit': 'M',  # Standard unit for Km
                        'parameter_type': 'KINETIC_CONSTANT'
                    })
                elif 'kcat' in pname_lower and '/' not in pname_lower:
                    # Turnover number: kcat
                    params.append({
                        'type': 'kcat',
                        'name': 'kcat',
                        'value': pval,
                        'unit': '1/s',  # Standard unit for kcat
                        'parameter_type': 'KINETIC_CONSTANT'
                    })
                elif 'vmax' in pname_lower:
                    # Maximum velocity: Vmax
                    # Note: Vmax unit can vary, use from CSV if available
                    punit = clean_unit(r.get('Parameter.unit') or r.get('parameter.unit'))
                    params.append({
                        'type': 'vmax',
                        'name': 'Vmax',
                        'value': pval,
                        'unit': punit,  # Use unit from CSV for Vmax
                        'parameter_type': 'KINETIC_CONSTANT'
                    })
            
            # Compose document with all extracted information
            doc = self.composer.compose(
                row=row,
                enzyme_name=enzyme_name,
                enzyme_variant=enzyme_variant,
                substrates=substrates,
                parameters=params,
                entry_id=entry_id,
                ph=ph,
                temperature=temperature,
                buffer_composition=buffer_comp,
                buffer_concentration=buffer_conc,
                reactants_parsed=reactants_parsed,
                products_parsed=products_parsed
            )
            
            # Complete to full schema (fills all missing fields with null)
            completed = complete_doc_to_schema(doc, self.schema_classes)
            
            # Clean unwanted fields to match BRENDA exactly
            completed = clean_unwanted_fields(completed)
            
            # Generate filename: SABIO_EC_organism_uniprot_variant_substrate (BRENDA-style)
            filename_parts = ['SABIO']
            
            # Extract EC number (try various column names)
            ec_number = None
            for ec_col in ['Protein.ecnumber', 'ecnumber', 'EC', 'ec', 'Enzyme.ecnumber']:
                if ec_col in row.index:
                    ec_val = row[ec_col]
                    if ec_val and pd.notna(ec_val) and str(ec_val).strip() not in ['', 'n/a', 'nan', 'none']:
                        ec_number = str(ec_val).strip()
                        break
            
            # Add EC number
            if ec_number:
                filename_parts.append(f"EC{sanitize_id(ec_number)}")
            
            # Extract organism (try various column names)
            organism = None
            for org_col in ['Protein.organism', 'organism', 'Organism']:
                if org_col in row.index:
                    org_val = row[org_col]
                    if org_val and pd.notna(org_val) and str(org_val).strip() not in ['', 'n/a', 'nan', 'none']:
                        organism = str(org_val).strip()
                        break
            
            # Extract UniProt ID
            uniprot_id = None
            for up_col in ['Protein.uniprotid', 'uniprotid', 'UniProt', 'uniprot']:
                if up_col in row.index:
                    up_val = row[up_col]
                    if up_val and pd.notna(up_val) and str(up_val).strip() not in ['', 'n/a', 'nan', 'none']:
                        uniprot_id = str(up_val).strip()
                        break
            
            # Add organism (first two words: genus species)
            if organism:
                org_parts = organism.split()[:2]  # Take genus species
                org_name = "_".join([sanitize_id(p) for p in org_parts])
                filename_parts.append(org_name)
            
            # Add UniProt ID (like BRENDA: always add if present)
            if uniprot_id:
                filename_parts.append(f"_{sanitize_id(uniprot_id)}")
            
            # Add variant if exists
            if enzyme_variant:
                filename_parts.append(sanitize_id(enzyme_variant))
            
            # Add substrate (first one, simplified)
            if substrates:
                substrate_name = substrates[0]
                # Take first word or up to 20 characters
                if ' ' in substrate_name:
                    substrate_name = substrate_name.split()[0]
                substrate_name = substrate_name[:20]  # Limit length
                filename_parts.append("s_" + sanitize_id(substrate_name))  # Add 's_' prefix like BRENDA
            
            # Join parts
            filename = "_".join(filename_parts)
            
            # Ensure uniqueness
            original_filename = filename
            counter = 1
            while True:
                test_name = f"{filename}.yaml" if format != 'json' else f"{filename}.json"
                test_path = output_path / test_name
                if not test_path.exists():
                    break
                counter += 1
                filename = f"{original_filename}_{counter}"
            
            # Write file
            if format in ['json', 'both']:
                import json
                json_file = output_path / f"{filename}.json"
                with open(json_file, 'w', encoding='utf-8') as f:
                    json.dump(completed, f, ensure_ascii=False, indent=2)
            
            if format in ['yaml', 'both']:
                if HAVE_YAML:
                    yaml_file = output_path / f"{filename}.yaml"
                    with open(yaml_file, 'w', encoding='utf-8') as f:
                        yaml.safe_dump(completed, f, allow_unicode=True, sort_keys=False, default_flow_style=False)
            
            count += 1
            if count % 100 == 0:
                print(f"  Generated {count} files...")
        
        print(f"\n✅ Complete! Generated {count} files")
        print(f"  Output: {output_path.resolve()}")
        print(f"  Format: EXACT BRENDA match (UNIFIED)")


# ============================================================
# MAIN
# ============================================================

def load_schema_classes(schema_path: Path) -> dict:
    """Load schema classes from YAML file"""
    if not HAVE_YAML:
        raise RuntimeError("PyYAML is required. Install: pip install pyyaml")
    with open(schema_path, "r", encoding="utf-8") as f:
        sch = yaml.safe_load(f)
    return sch.get("classes", {})


def main():
    parser = argparse.ArgumentParser(
        description="SABIO-RK to EnzymeML - EXACT BRENDA Format (UNIFIED VERSION)"
    )
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--schema', required=True, help='Schema YAML file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--format', choices=['yaml', 'json', 'both'],
                       default='yaml', help='Output format')
    parser.add_argument('--group-by', choices=['entry', 'protein_substrate'],
                       default='entry', 
                       help='Grouping strategy: "entry" (one file per EntryID) or "protein_substrate" (merge same protein+substrate)')
    parser.add_argument('--limit', type=int, default=0, 
                       help='Limit number of entries to process (0=no limit)')
    args = parser.parse_args()
    
    if args.format in ['yaml', 'both'] and not HAVE_YAML:
        print("Error: PyYAML not installed. Install: pip install pyyaml")
        args.format = 'json'
    
    # Load schema
    schema_registry = SchemaLoader(Path(args.schema))
    schema_classes = load_schema_classes(Path(args.schema))
    print(f"✅ Loaded schema with {len(schema_registry.classes)} classes")
    
    # Load data
    print(f"\nLoading data: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Loaded {len(df)} rows")
    
    # Check required columns
    if 'EntryID' not in df.columns:
        raise ValueError("ERROR: 'EntryID' column not found")
    
    # Generate files
    generator = YAMLGenerator(schema_classes)
    generator.generate(df, args.output_dir, args.format, args.limit, args.group_by)


if __name__ == '__main__':
    main()
