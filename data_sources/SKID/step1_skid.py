#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SKiD Step1: Convert to EnzymeML Format (BRENDA-compatible output)
======================================================================
Modified from SABIO-RK to output EXACTLY the same format as BRENDA step1_join.py
UNIFIED VERSION - Output format identical to BRENDA/SABIO-RK

SKiD (Structure-Oriented Kinetics Dataset) Main Dataset Processing
- Intelligently merges kcat_dataset and Km_dataset sheets
- Uses (UniProt_ID, Substrate, Mutation) as the merge key
- Generates one file per enzyme-substrate-mutation combination

Usage:
python step1_skid_main.py --input Main_dataset_v1.xlsx --schema ../../schemas/enzymeml-v2-extended.yaml --output-dir ../../output/skid --format yaml
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


def extract_enzyme_name_from_data(ec_number, organism, uniprot_id):
    """Generate enzyme name from available data"""
    if not ec_number or pd.isna(ec_number):
        if uniprot_id and pd.notna(uniprot_id):
            name = f"Enzyme ({uniprot_id})"
        else:
            name = "Enzyme"
    else:
        name = f"EC {ec_number}"
    
    # Add organism if available
    if organism and pd.notna(organism):
        organism_str = str(organism).strip()
        # Take genus and species (first two words)
        org_parts = organism_str.split()[:2]
        if org_parts:
            name += f" ({' '.join(org_parts)})"
    
    return name


def parse_mutation(mutation_str, mutant_flag):
    """
    Parse mutation information
    Returns: (is_wildtype, mutation_str)
    """
    if pd.isna(mutation_str) or str(mutation_str).strip() in ['-----', '', 'nan']:
        return True, None
    
    if pd.isna(mutant_flag):
        return True, None
        
    mutant_flag_str = str(mutant_flag).strip().lower()
    if mutant_flag_str == 'no':
        return True, None
    
    mutation = str(mutation_str).strip()
    if mutation in ['-----', '']:
        return True, None
        
    return False, mutation


# ============================================================
# DATA MERGER
# ============================================================

class SKiDDataMerger:
    """Intelligently merge kcat and Km datasets"""
    
    def __init__(self, kcat_df: pd.DataFrame, km_df: pd.DataFrame):
        self.kcat_df = kcat_df
        self.km_df = km_df
        
    def create_merge_key(self, df: pd.DataFrame) -> pd.Series:
        """Create merge key: UniProt_ID|Substrate|Mutation"""
        return (
            df['UniProt_ID'].astype(str) + '|' +
            df['Substrate'].astype(str) + '|' +
            df['Mutation'].fillna('-----').astype(str)
        )
    
    def merge_data(self):
        """
        Merge kcat and Km data using (UniProt_ID, Substrate, Mutation) as key
        Returns: merged dataframe with both kcat and Km where available
        """
        # Add merge keys
        self.kcat_df['merge_key'] = self.create_merge_key(self.kcat_df)
        self.km_df['merge_key'] = self.create_merge_key(self.km_df)
        
        # Merge on the key
        merged = pd.merge(
            self.kcat_df,
            self.km_df,
            on='merge_key',
            how='outer',
            suffixes=('_kcat', '_km')
        )
        
        # Consolidate common columns (prefer kcat version if both exist)
        common_cols = ['EC_number', 'Substrate', 'UniProt_ID', 'Organism_name', 
                      'Mutant', 'Mutation', 'pH', 'Temperature', 'Substrate_SMILES']
        
        for col in common_cols:
            col_kcat = f"{col}_kcat"
            col_km = f"{col}_km"
            if col_kcat in merged.columns and col_km in merged.columns:
                # Coalesce: use kcat value if available, otherwise km value
                merged[col] = merged[col_kcat].fillna(merged[col_km])
                # Drop the suffixed columns
                merged = merged.drop(columns=[col_kcat, col_km])
            elif col_kcat in merged.columns:
                merged = merged.rename(columns={col_kcat: col})
            elif col_km in merged.columns:
                merged = merged.rename(columns={col_km: col})
        
        return merged


# ============================================================
# DOCUMENT COMPOSER
# ============================================================

class SKiDDocumentComposer:
    """Compose EnzymeML documents from merged SKiD data"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
    
    def compose(
        self,
        row: dict,
        enzyme_name: str,
        ec_number: str,
        organism: str,
        uniprot_id: str,
        substrate: str,
        substrate_smiles: str,
        mutation: str,
        parameters: list,
        entry_id: str,
        ph: float,
        temperature: float,
        is_wildtype: bool = False,
        **kwargs
    ) -> dict:
        """Compose a complete EnzymeML document"""
        
        # Create document
        doc = {
            'name': f'SKiD_{entry_id}',
            'pubmedid': None,
            'doi': None,
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
                'ecnumber': ec_number if ec_number and pd.notna(ec_number) else None,
                'organism': organism if organism and pd.notna(organism) else None,
                'organism_tax_id': None,
                'uniprotid': uniprot_id if uniprot_id and pd.notna(uniprot_id) else None,
                'constant': None,
                'meta_id': None
            }
            proteins.append(protein)
        doc['proteins'] = proteins if proteins else None
        
        # Add reactants (substrate)
        reactants_list = []
        substrate_id = None
        if substrate:
            reactant = {
                'id': 's0',
                'name': substrate,
                'smiles': substrate_smiles if substrate_smiles and pd.notna(substrate_smiles) else None,
                'inchi': None,
                'inchikey': None,
                'chebi': None,
                'constant': None,
                'meta_id': None
            }
            reactants_list.append(reactant)
            substrate_id = 's0'
        
        doc['reactants'] = reactants_list if reactants_list else None
        
        # Add reactions
        reactions = []
        if enzyme_name or substrate:
            # Create reaction
            reaction = {
                'id': 'r0',
                'name': f'Enzymatic conversion of {substrate}' if substrate else 'Enzymatic reaction',
                'reversible': False,
                'temperature_value': float(temperature) if temperature and pd.notna(temperature) else None,
                'temperature_unit': 'C' if temperature and pd.notna(temperature) else None,
                'ph': float(ph) if ph and pd.notna(ph) else None,
                'ontology': None,
                'meta_id': None
            }
            
            # Add educts (substrates)
            educts = []
            if substrate_id:
                educt = {
                    'species_id': substrate_id,
                    'stoichiometry': 1.0,
                    'constant': None,
                    'ontology': None
                }
                educts.append(educt)
            reaction['educts'] = educts if educts else None
            
            # Products - not specified in SKiD data
            reaction['products'] = None
            
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
        
        return doc


# ============================================================
# YAML GENERATOR
# ============================================================

class YAMLGenerator:
    """Generate YAML files from merged SKiD data"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
        self.composer = SKiDDocumentComposer(schema_classes)
    
    def generate(self, merged_df: pd.DataFrame, output_dir: str, format: str = 'yaml', 
                limit: int = 0, output_mode: str = 'both'):
        """
        Generate YAML/JSON files from merged data
        
        Args:
            output_mode: 'both' (WT and mutant), 'wt_only' (wild-type only), 'mutant_only' (mutant only)
        """
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"\n{'='*60}")
        print(f"Generating SKiD EnzymeML files (Main Dataset)")
        print(f"{'='*60}")
        print(f"  Output: {output_path.resolve()}")
        print(f"  Format: {format}")
        print(f"  Output mode: {output_mode}")
        print(f"  Total rows: {len(merged_df)}")
        if limit > 0:
            print(f"  Limit: {limit} entries")
        print(f"{'='*60}\n")
        
        # Apply limit if specified
        if limit > 0:
            merged_df = merged_df.head(limit)
        
        count = 0
        for idx, row in merged_df.iterrows():
            # Generate unique entry ID
            entry_id = f"SKID_M_{count:06d}"
            
            # Extract enzyme information
            ec_number = row.get('EC_number')
            organism = row.get('Organism_name')
            uniprot_id = row.get('UniProt_ID')
            enzyme_name = extract_enzyme_name_from_data(ec_number, organism, uniprot_id)
            
            # Extract substrate information
            substrate = row.get('Substrate')
            substrate_smiles = row.get('Substrate_SMILES')
            
            # Extract mutation information
            is_wildtype, mutation = parse_mutation(row.get('Mutation'), row.get('Mutant'))
            
            # Check output mode
            if output_mode == 'wt_only' and not is_wildtype:
                continue
            if output_mode == 'mutant_only' and is_wildtype:
                continue
            
            # Extract conditions
            ph = row.get('pH')
            temperature = row.get('Temperature')
            
            # Extract kinetic parameters
            params = []
            
            # kcat (from kcat_dataset)
            kcat_val = row.get('kcat_value (1/s)')
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
            
            # Km (from Km_dataset)
            km_val = row.get('Km_value (mM)')
            if km_val and pd.notna(km_val):
                try:
                    km_float = float(km_val)
                    params.append({
                        'type': 'km',
                        'name': 'Km',
                        'value': km_float,
                        'unit': 'mM',
                        'parameter_type': 'KINETIC_CONSTANT'
                    })
                except:
                    pass
            
            # Calculate kcat/Km if both are available
            if kcat_val and pd.notna(kcat_val) and km_val and pd.notna(km_val):
                try:
                    kcat_float = float(kcat_val)
                    km_float = float(km_val)
                    if km_float != 0:
                        # kcat/Km (catalytic efficiency)
                        # Convert Km from mM to M: mM * 0.001 = M
                        km_in_M = km_float * 0.001
                        kcat_km = kcat_float / km_in_M
                        params.append({
                            'type': 'kcat_km',
                            'name': 'kcat/Km',
                            'value': kcat_km,
                            'unit': '1/(M*s)',
                            'parameter_type': 'SPECIFICITY_CONSTANT'
                        })
                except:
                    pass
            
            # Skip if no parameters
            if not params:
                continue
            
            # Compose document with all extracted information
            doc = self.composer.compose(
                row=row,
                enzyme_name=enzyme_name,
                ec_number=ec_number,
                organism=organism,
                uniprot_id=uniprot_id,
                substrate=substrate,
                substrate_smiles=substrate_smiles,
                mutation=mutation,
                parameters=params,
                entry_id=entry_id,
                ph=ph,
                temperature=temperature,
                is_wildtype=is_wildtype
            )
            
            # Complete to full schema (fills all missing fields with null)
            completed = complete_doc_to_schema(doc, self.schema_classes)
            
            # Clean unwanted fields to match BRENDA exactly
            completed = clean_unwanted_fields(completed)
            
            # Generate filename: SKID_EC_uniprot_organism_substrate_mutation
            filename_parts = ['SKID']
            
            # Add EC number
            if ec_number and pd.notna(ec_number):
                filename_parts.append(f"EC{sanitize_id(str(ec_number))}")
            
            # Add UniProt ID (shorter and more specific than organism)
            if uniprot_id and pd.notna(uniprot_id):
                filename_parts.append(sanitize_id(str(uniprot_id)))
            
            # Add organism (genus species) if no UniProt
            elif organism and pd.notna(organism):
                org_parts = str(organism).split()[:2]  # Take genus species
                org_name = "_".join([sanitize_id(p) for p in org_parts])
                filename_parts.append(org_name)
            
            # Add substrate (simplified)
            if substrate and pd.notna(substrate):
                substrate_short = str(substrate).split()[0] if ' ' in str(substrate) else str(substrate)
                substrate_short = substrate_short[:20]  # Limit length
                filename_parts.append("s_" + sanitize_id(substrate_short))
            
            # Add mutation or WT marker
            if mutation:
                filename_parts.append(sanitize_id(mutation))
            else:
                filename_parts.append("WT")
            
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
            if count % 500 == 0:
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
        description="SKiD Main Dataset to EnzymeML - EXACT BRENDA/SABIO-RK Format (UNIFIED VERSION)"
    )
    parser.add_argument('--input', required=True, help='Input Excel file (Main_dataset_v1.xlsx)')
    parser.add_argument('--schema', required=True, help='Schema YAML file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--format', choices=['yaml', 'json', 'both'],
                       default='yaml', help='Output format')
    parser.add_argument('--output-mode', choices=['both', 'wt_only', 'mutant_only'],
                       default='both',
                       help='Output mode: both (WT and mutant), wt_only, or mutant_only')
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
    
    # Load data from both sheets
    print(f"\nLoading data: {args.input}")
    print("  Reading kcat_dataset sheet...")
    df_kcat = pd.read_excel(args.input, sheet_name='kcat_dataset')
    print(f"  Loaded {len(df_kcat)} kcat records")
    
    print("  Reading Km_dataset sheet...")
    df_km = pd.read_excel(args.input, sheet_name='Km_dataset')
    print(f"  Loaded {len(df_km)} Km records")
    
    # Merge data intelligently
    print("\n  Merging kcat and Km data using (UniProt_ID, Substrate, Mutation) as key...")
    merger = SKiDDataMerger(df_kcat, df_km)
    merged_df = merger.merge_data()
    print(f"  Merged dataset: {len(merged_df)} combinations")
    
    # Count how many have both kcat and Km
    has_both = merged_df[
        merged_df['kcat_value (1/s)'].notna() & 
        merged_df['Km_value (mM)'].notna()
    ]
    print(f"  Combinations with both kcat and Km: {len(has_both)}")
    print(f"  Combinations with only kcat: {merged_df['kcat_value (1/s)'].notna().sum() - len(has_both)}")
    print(f"  Combinations with only Km: {merged_df['Km_value (mM)'].notna().sum() - len(has_both)}")
    
    # Generate files
    generator = YAMLGenerator(schema_classes)
    generator.generate(merged_df, args.output_dir, args.format, args.limit, args.output_mode)


if __name__ == '__main__':
    main()
