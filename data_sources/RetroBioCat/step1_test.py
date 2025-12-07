#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RetroBioCat Step1: Convert to EnzymeML Format (BRENDA-compatible output)
======================================================================
Modified from SABIO-RK to output EXACTLY the same format as BRENDA step1_join.py
UNIFIED VERSION - Output format identical to BRENDA/SABIO-RK

Usage:
python step1_test.py --input ../../data/RetroBioCat/trial_activity_data.xlsx --schema ../../schemas/enzymeml-v2-extended.yaml --output-dir ../../output/retrobiocat --format yaml
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
# SCHEMA COMPLETION (EXACT COPY FROM BRENDA/SABIO-RK)
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
    """Clean unit string - fix common encoding issues"""
    if not unit_str or pd.isna(unit_str):
        return None
    
    unit = str(unit_str).strip()
    
    # Fix common patterns
    if unit == 's?1':
        return '1/s'
    if unit == 'M?1s?1':
        return '1/(M*s)'
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
    # Remove mutant/variant/wild-type info
    name = re.sub(r'\([^)]*mutant[^)]*\)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\([^)]*variant[^)]*\)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\([^)]*wild[- ]?type[^)]*\)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\s+(mutant|variant|wild[- ]?type).*$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\s+', ' ', name).strip()
    return name if name else None


def extract_substrates_from_smiles(smiles_str):
    """Extract substrate name from SMILES (use SMILES as identifier)"""
    if not smiles_str or pd.isna(smiles_str):
        return None
    return str(smiles_str).strip()


def extract_products_from_smiles(smiles_str):
    """Extract product name from SMILES (use SMILES as identifier)"""
    if not smiles_str or pd.isna(smiles_str):
        return None
    return str(smiles_str).strip()


def parse_temperature(temp_str):
    """Parse temperature string to float"""
    if not temp_str or pd.isna(temp_str):
        return None
    
    temp = str(temp_str).strip()
    # Extract numeric value
    match = re.search(r'(\d+(?:\.\d+)?)', temp)
    if match:
        try:
            return float(match.group(1))
        except:
            return None
    return None


def parse_buffer(solvent_str, other_conditions_str):
    """Parse buffer composition from solvent and other conditions"""
    parts = []
    
    if solvent_str and pd.notna(solvent_str) and str(solvent_str).strip():
        parts.append(str(solvent_str).strip())
    
    if other_conditions_str and pd.notna(other_conditions_str) and str(other_conditions_str).strip():
        parts.append(str(other_conditions_str).strip())
    
    if parts:
        return "; ".join(parts)
    return None


# ============================================================
# DOCUMENT COMPOSER
# ============================================================

class RetroBioCatDocumentComposer:
    """Compose EnzymeML documents from RetroBioCat data"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
    
    def compose(
        self,
        row: dict,
        enzyme_name: str,
        enzyme_type: str,
        substrates: list,
        products: list,
        parameters: list,
        entry_id: str,
        ph: float,
        temperature: float,
        buffer_composition: str,
        citation: str,
        doi: str,
        **kwargs
    ) -> dict:
        """Compose a complete EnzymeML document"""
        
        # Create document
        doc = {
            'name': f'RetroBioCat_{entry_id}',
            'pubmedid': None,
            'doi': doi,
            'url': None,
            'created': datetime.now().strftime('%Y-%m-%d'),
            'modified': datetime.now().strftime('%Y-%m-%d')
        }
        
        # Add proteins (enzymes)
        proteins = []
        if enzyme_name:
            protein = {
                'id': f'p0',
                'name': enzyme_name,
                'sequence': None,
                'ecnumber': None,
                'organism': None,
                'organism_tax_id': None,
                'uniprotid': None,
                'constant': None,
                'meta_id': None
            }
            proteins.append(protein)
        doc['proteins'] = proteins if proteins else None
        
        # Add reactants (substrates and products from SMILES)
        reactants_list = []
        reactant_counter = 0
        
        # Add substrates
        substrate_ids = []
        if substrates:
            for sub in substrates:
                if sub:
                    reactant = {
                        'id': f's{reactant_counter}',
                        'name': f'Substrate_{reactant_counter}',
                        'smiles': sub,
                        'inchi': None,
                        'inchikey': None,
                        'chebi': None,
                        'constant': None,
                        'meta_id': None
                    }
                    reactants_list.append(reactant)
                    substrate_ids.append(f's{reactant_counter}')
                    reactant_counter += 1
        
        # Add products
        product_ids = []
        if products:
            for prod in products:
                if prod:
                    reactant = {
                        'id': f'p{reactant_counter}',
                        'name': f'Product_{reactant_counter}',
                        'smiles': prod,
                        'inchi': None,
                        'inchikey': None,
                        'chebi': None,
                        'constant': None,
                        'meta_id': None
                    }
                    reactants_list.append(reactant)
                    product_ids.append(f'p{reactant_counter}')
                    reactant_counter += 1
        
        doc['reactants'] = reactants_list if reactants_list else None
        
        # Add reactions
        reactions = []
        if enzyme_name or substrates or products:
            # Create reaction
            reaction = {
                'id': 'r0',
                'name': row.get('reaction', 'Biocatalytic reaction'),
                'reversible': False,
                'temperature_value': temperature,
                'temperature_unit': 'C' if temperature else None,
                'ph': ph,
                'ontology': None,
                'meta_id': None
            }
            
            # Add educts (substrates)
            educts = []
            for sid in substrate_ids:
                educt = {
                    'species_id': sid,
                    'stoichiometry': 1.0,
                    'constant': None,
                    'ontology': None
                }
                educts.append(educt)
            reaction['educts'] = educts if educts else None
            
            # Add products
            prods = []
            for pid in product_ids:
                prod = {
                    'species_id': pid,
                    'stoichiometry': 1.0,
                    'constant': None,
                    'ontology': None
                }
                prods.append(prod)
            reaction['products'] = prods if prods else None
            
            # Add modifiers (enzyme)
            modifiers = []
            if enzyme_name:
                modifier = {
                    'species_id': 'p0',
                    'stoichiometry': None,
                    'constant': None,
                    'ontology': None
                }
                modifiers.append(modifier)
            reaction['modifiers'] = modifiers if modifiers else None
            
            reactions.append(reaction)
        
        doc['reactions'] = reactions if reactions else None
        
        # Add conditions
        conditions = []
        
        # pH condition
        if ph and pd.notna(ph):
            ph_condition = {
                'id': 'c0',
                'species_id': None,
                'temperature_value': None,
                'temperature_unit': None,
                'ph': float(ph),
                'created': None,
                'creator_family_name': None,
                'creator_given_name': None,
                'modified': None,
                'meta_id': None
            }
            conditions.append(ph_condition)
        
        # Temperature condition
        if temperature and pd.notna(temperature):
            temp_condition = {
                'id': 'c1',
                'species_id': None,
                'temperature_value': float(temperature),
                'temperature_unit': 'C',
                'ph': None,
                'created': None,
                'creator_family_name': None,
                'creator_given_name': None,
                'modified': None,
                'meta_id': None
            }
            conditions.append(temp_condition)
        
        # Buffer condition
        if buffer_composition:
            buffer_condition = {
                'id': 'c2',
                'species_id': None,
                'temperature_value': None,
                'temperature_unit': None,
                'ph': None,
                'created': None,
                'creator_family_name': None,
                'creator_given_name': None,
                'modified': None,
                'meta_id': None
            }
            conditions.append(buffer_condition)
        
        doc['conditions'] = conditions if conditions else None
        
        # Add parameters
        if parameters:
            doc['parameters'] = parameters
        else:
            doc['parameters'] = None
        
        # Add measurements (empty for now, but include structure)
        doc['measurements'] = None
        
        # Add files and references
        doc['files'] = None
        doc['references'] = None
        
        # Add creator info
        doc['creator_family_name'] = None
        doc['creator_given_name'] = None
        
        # Add citation if available
        if citation:
            doc['references'] = [{
                'id': 'ref0',
                'doi': doi,
                'pubmed_id': None,
                'url': None
            }]
        
        return doc


# ============================================================
# YAML GENERATOR
# ============================================================

class YAMLGenerator:
    """Generate YAML files from RetroBioCat data"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
        self.composer = RetroBioCatDocumentComposer(schema_classes)
    
    def generate(self, df: pd.DataFrame, output_dir: str, format: str = 'yaml', 
                limit: int = 0, group_by: str = 'entry'):
        """Generate YAML/JSON files"""
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"\n{'='*60}")
        print(f"Generating RetroBioCat EnzymeML files")
        print(f"{'='*60}")
        print(f"  Output: {output_path.resolve()}")
        print(f"  Format: {format}")
        print(f"  Total rows: {len(df)}")
        if limit > 0:
            print(f"  Limit: {limit} entries")
        print(f"{'='*60}\n")
        
        # Apply limit if specified
        if limit > 0:
            df = df.head(limit)
        
        count = 0
        for idx, row in df.iterrows():
            # Generate unique entry ID
            entry_id = f"RBC_{idx:06d}"
            
            # Extract enzyme information
            enzyme_name = extract_enzyme_name(row.get('enzyme_name'))
            enzyme_type = row.get('enzyme_type')
            
            # Extract substrates and products from SMILES
            substrates = []
            products = []
            
            # Substrate 1
            sub1_smiles = extract_substrates_from_smiles(row.get('substrate_1_smiles'))
            if sub1_smiles:
                substrates.append(sub1_smiles)
            
            # Substrate 2
            sub2_smiles = extract_substrates_from_smiles(row.get('substrate_2_smiles'))
            if sub2_smiles:
                substrates.append(sub2_smiles)
            
            # Product 1
            prod1_smiles = extract_products_from_smiles(row.get('product_1_smiles'))
            if prod1_smiles:
                products.append(prod1_smiles)
            
            # Extract conditions
            ph = row.get('ph')
            temperature = parse_temperature(row.get('temperature'))
            buffer_composition = parse_buffer(row.get('solvent'), row.get('other_conditions'))
            
            # Extract citation
            citation = row.get('short_citation')
            doi = row.get('html_doi')
            
            # Extract kinetic parameters
            params = []
            
            # kcat
            kcat_val = row.get('kcat')
            if kcat_val and pd.notna(kcat_val):
                try:
                    kcat_float = float(kcat_val)
                    params.append({
                        'type': 'kcat',
                        'name': 'kcat',
                        'value': kcat_float,
                        'unit': '1/s',
                        'parameter_type': 'KINETIC_CONSTANT'
                    })
                except:
                    pass
            
            # Km
            km_val = row.get('km')
            if km_val and pd.notna(km_val):
                try:
                    km_float = float(km_val)
                    params.append({
                        'type': 'km',
                        'name': 'Km',
                        'value': km_float,
                        'unit': 'M',
                        'parameter_type': 'KINETIC_CONSTANT'
                    })
                except:
                    pass
            
            # Specific activity
            spec_act = row.get('specific_activity')
            if spec_act and pd.notna(spec_act):
                try:
                    spec_act_float = float(spec_act)
                    params.append({
                        'type': 'specific_activity',
                        'name': 'Specific Activity',
                        'value': spec_act_float,
                        'unit': 'U/mg',
                        'parameter_type': 'ACTIVITY'
                    })
                except:
                    pass
            
            # Conversion
            conversion = row.get('conversion')
            if conversion and pd.notna(conversion):
                try:
                    conversion_float = float(conversion)
                    params.append({
                        'type': 'conversion',
                        'name': 'Conversion',
                        'value': conversion_float,
                        'unit': '%',
                        'parameter_type': 'YIELD'
                    })
                except:
                    pass
            
            # Compose document with all extracted information
            doc = self.composer.compose(
                row=row,
                enzyme_name=enzyme_name,
                enzyme_type=enzyme_type,
                substrates=substrates,
                products=products,
                parameters=params,
                entry_id=entry_id,
                ph=ph,
                temperature=temperature,
                buffer_composition=buffer_composition,
                citation=citation,
                doi=doi
            )
            
            # Complete to full schema (fills all missing fields with null)
            completed = complete_doc_to_schema(doc, self.schema_classes)
            
            # Clean unwanted fields to match BRENDA exactly
            completed = clean_unwanted_fields(completed)
            
            # Generate filename: RBC_enzyme_type_enzyme_name_reaction
            filename_parts = ['RBC']
            
            # Add enzyme type
            if enzyme_type:
                filename_parts.append(sanitize_id(enzyme_type))
            
            # Add enzyme name (simplified)
            if enzyme_name:
                name_short = enzyme_name.split()[0] if ' ' in enzyme_name else enzyme_name
                name_short = name_short[:20]  # Limit length
                filename_parts.append(sanitize_id(name_short))
            
            # Add reaction type
            reaction_type = row.get('reaction')
            if reaction_type:
                reaction_short = reaction_type.split()[0] if ' ' in reaction_type else reaction_type
                reaction_short = reaction_short[:20]
                filename_parts.append(sanitize_id(reaction_short))
            
            # Add entry ID for uniqueness
            filename_parts.append(entry_id)
            
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
        print(f"  Format: EXACT BRENDA/SABIO-RK match (UNIFIED)")


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
        description="RetroBioCat to EnzymeML - EXACT BRENDA/SABIO-RK Format (UNIFIED VERSION)"
    )
    parser.add_argument('--input', required=True, help='Input Excel file')
    parser.add_argument('--schema', required=True, help='Schema YAML file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--format', choices=['yaml', 'json', 'both'],
                       default='yaml', help='Output format')
    parser.add_argument('--group-by', choices=['entry', 'enzyme_substrate'],
                       default='entry', 
                       help='Grouping strategy: "entry" (one file per entry) or "enzyme_substrate" (merge same enzyme+substrate)')
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
    df = pd.read_excel(args.input)
    print(f"  Loaded {len(df)} rows")
    
    # Generate files
    generator = YAMLGenerator(schema_classes)
    generator.generate(df, args.output_dir, args.format, args.limit, args.group_by)


if __name__ == '__main__':
    main()
