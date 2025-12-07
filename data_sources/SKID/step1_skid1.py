#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SKiD Step1: Convert to EnzymeML Format - GROUP BY (UniProt, Substrate, Mutation)
=================================================================================
CORRECT GROUPING LOGIC:
- Same UniProt + Substrate + Mutation → SAME YAML file
- Different Mutation → DIFFERENT YAML file
- Each file can contain multiple kcat/Km values for the same combination

IMPROVEMENTS:
1. Group by (UniProt_ID, Substrate, Mutation) - correct three-way grouping
2. Load Unique_substrates sheet to get InChI, InChIKey, and External_identifiers
3. Merge substrate info into kcat/Km datasets based on substrate name
4. External_identifiers are stored in small_molecules.notes field
5. Reference parsing: Extract title, authors, year, journal, pubmed_id, doi
6. Add pubmed_id to each kcat/Km parameter
7. Use organism_name for proteins.name
8. Add product as placeholder (6 small_molecules total)
9. Shortened filenames to avoid Windows 260-character path limit

Usage:
python step1_skid_with_substrate_info.py --input Main_dataset_v1.xlsx --schema ../../schemas/enzymeml-v2-extended.yaml --output-dir ../../output/skid --format yaml
"""

import argparse
import re
import sys
import ast
from pathlib import Path
from datetime import datetime
from collections import defaultdict

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
# SCHEMA COMPLETION
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
    """Complete document to match full schema"""
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
    """Keep fields but set to null as requested"""
    if 'equations' in doc:
        doc['equations'] = None  # <- 不删除，置空
    # 保险：确保每个 reaction 的 equation 也是 null
    if isinstance(doc.get('reactions'), list):
        for rxn in doc['reactions']:
            if isinstance(rxn, dict):
                rxn['equation'] = None
    return doc



# ============================================================
# HELPER FUNCTIONS
# ============================================================

def sanitize_id(s: str) -> str:
    """Sanitize string for use as ID"""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def parse_reference_citation(citation_str):
    """
    Parse reference citation string to extract structured information
    
    Format: "Authors, Title, Journal, Year, Pages, Volume, PMID: xxxxx"
    Example: "Huang, L.,Sayoga, G., Title, ACS Catal., 2018, 8680-8684, 8, PMID: 12345"
    
    Returns: dict with title, authors, year, journal, pages, volume, pubmed_id
    """
    if not citation_str or pd.isna(citation_str):
        return None
    
    citation = str(citation_str).strip()
    
    # Extract PMID first
    pmid = None
    pmid_match = re.search(r'PMID:\s*(\d+)', citation)
    if pmid_match:
        pmid = pmid_match.group(1)
        # Remove PMID part from citation
        citation = re.sub(r',?\s*PMID:\s*\d+', '', citation)
    
    # Split by comma
    parts = [p.strip() for p in citation.split(',')]
    
    if len(parts) < 4:
        # Not enough parts, return basic info
        return {
            'title': None,
            'authors': None,
            'year': None,
            'journal': None,
            'pages': None,
            'volume': None,
            'pubmed_id': pmid,
            'doi': None,
            'uri': None
        }
    
    # Try to parse structure
    # Common format: Author1, Author2, ..., Title, Journal, Year, Pages, Volume
    
    # Find year (4-digit number)
    year_idx = None
    year = None
    for i, part in enumerate(parts):
        if re.match(r'^\d{4}$', part):
            year_idx = i
            year = part
            break
    
    if year_idx is None:
        # No clear year found
        return {
            'title': None,
            'authors': ', '.join(parts[:-2]) if len(parts) > 2 else None,
            'year': None,
            'journal': parts[-1] if parts else None,
            'pages': None,
            'volume': None,
            'pubmed_id': pmid,
            'doi': None,
            'uri': None
        }
    
    # Extract components
    # Before year: authors, title, journal
    # After year: pages, volume
    before_year = parts[:year_idx]
    after_year = parts[year_idx+1:]
    
    # Heuristic: journal is likely 1-2 parts before year
    journal_idx = max(0, year_idx - 2)
    journal = parts[journal_idx] if journal_idx < year_idx else None
    
    # Title is likely 1 part before journal
    title_idx = max(0, journal_idx - 1)
    title = parts[title_idx] if title_idx < journal_idx else None
    
    # Authors are everything before title
    authors = ', '.join(parts[:title_idx]) if title_idx > 0 else None
    
    # Pages and volume
    pages = after_year[0] if len(after_year) > 0 else None
    volume = after_year[1] if len(after_year) > 1 else None
    
    return {
        'title': title,
        'authors': authors,
        'year': year,
        'journal': journal,
        'pages': pages,
        'volume': volume,
        'pubmed_id': pmid,
        'doi': None,
        'uri': None
    }


def parse_references(ref_str):
    """
    Parse References string to extract structured citation information
    
    Format: "{'161': 'citation text'}"
    Returns: list of reference dictionaries with parsed fields
    """
    if not ref_str or pd.isna(ref_str):
        return None
    
    try:
        ref_dict = ast.literal_eval(str(ref_str))
        if isinstance(ref_dict, dict):
            references = []
            for ref_id, ref_text in ref_dict.items():
                # Parse the citation
                parsed = parse_reference_citation(ref_text)
                if parsed:
                    ref_entry = {
                        'id': str(ref_id),
                        'title': parsed.get('title'),
                        'authors': parsed.get('authors'),
                        'year': parsed.get('year'),
                        'journal': parsed.get('journal'),
                        'pages': parsed.get('pages'),
                        'volume': parsed.get('volume'),
                        'pubmed_id': parsed.get('pubmed_id'),
                        'doi': parsed.get('doi'),
                        'uri': parsed.get('uri'),
                        'citation': ref_text  # Keep full citation for reference
                    }
                else:
                    ref_entry = {
                        'id': str(ref_id),
                        'title': None,
                        'authors': None,
                        'year': None,
                        'journal': None,
                        'pages': None,
                        'volume': None,
                        'pubmed_id': None,
                        'doi': None,
                        'uri': None,
                        'citation': ref_text
                    }
                references.append(ref_entry)
            return references if references else None
    except:
        return [{
            'id': 'ref1',
            'title': None,
            'authors': None,
            'year': None,
            'journal': None,
            'pages': None,
            'volume': None,
            'pubmed_id': None,
            'doi': None,
            'uri': None,
            'citation': str(ref_str)
        }]
    
    return None


def extract_pmid_from_references(ref_str):
    """Extract PMID from references string"""
    if not ref_str or pd.isna(ref_str):
        return None
    
    # Try to extract PMID
    pmid_match = re.search(r'PMID:\s*(\d+)', str(ref_str))
    if pmid_match:
        return pmid_match.group(1)
    
    return None


def extract_enzyme_name_from_data(ec_number, organism, uniprot_id):
    """Generate enzyme name from organism"""
    # Use organism name if available
    if organism and pd.notna(organism):
        return str(organism).strip()
    
    # Fallback to EC or UniProt
    if ec_number and pd.notna(ec_number):
        return f"protein#{ec_number}"
    
    if uniprot_id and pd.notna(uniprot_id):
        return f"protein#{uniprot_id}"
    
    return "protein"


def parse_mutation(mutation_str, mutant_flag):
    """Parse mutation information"""
    if pd.isna(mutation_str) or str(mutation_str).strip() in ['-----', '', 'nan']:
        return True, None, 'wildtype', 'wildtype'
    
    if pd.isna(mutant_flag):
        return True, None, 'wildtype', 'wildtype'
        
    mutant_flag_str = str(mutant_flag).strip().lower()
    if mutant_flag_str == 'no':
        return True, None, 'wildtype', 'wildtype'
    
    mutation = str(mutation_str).strip()
    if mutation in ['-----', '']:
        return True, None, 'wildtype', 'wildtype'
        
    return False, mutation, 'mutant', mutation


def parse_temperature(temp_str) -> float:
    """Parse temperature value from string"""
    if not temp_str or pd.isna(temp_str):
        return None
    
    temp_str = str(temp_str).strip()
    
    if temp_str in ['-----', '', 'nan', 'none', 'n/a']:
        return None
    
    match = re.search(r'[-+]?\d*\.?\d+', temp_str)
    if match:
        try:
            return float(match.group())
        except:
            return None
    
    return None


def parse_ph(ph_str) -> float:
    """Parse pH value from string"""
    if not ph_str or pd.isna(ph_str):
        return None
    
    ph_str = str(ph_str).strip()
    
    if ph_str in ['-----', '', 'nan', 'none', 'n/a']:
        return None
    
    match = re.search(r'\d+\.?\d*', ph_str)
    if match:
        try:
            return float(match.group())
        except:
            return None
    
    return None


# ============================================================
# DATA MERGER CLASS
# ============================================================

class SKiDDataMerger:
    """Merge kcat and Km datasets"""
    
    def __init__(self, df_kcat: pd.DataFrame, df_km: pd.DataFrame):
        self.df_kcat = df_kcat
        self.df_km = df_km
    
    def merge_data(self) -> pd.DataFrame:
        """Merge kcat and Km data, keeping all records"""
        # Prepare kcat data
        df_kcat = self.df_kcat.copy()
        df_kcat['has_kcat'] = True
        df_kcat['has_km'] = False
        df_kcat['data_source'] = 'kcat'
        df_kcat = df_kcat.rename(columns={
            'Temperature': 'T_oC',
            'Substrate_SMILES': 'Smiles',
            'Organism_name': 'Organism'
        })
        
        # Prepare Km data
        df_km = self.df_km.copy()
        df_km['has_kcat'] = False
        df_km['has_km'] = True
        df_km['data_source'] = 'Km'
        df_km = df_km.rename(columns={
            'Temperature': 'T_oC',
            'Substrate_SMILES': 'Smiles',
            'Organism_name': 'Organism'
        })
        
        # Align columns
        all_cols = set(df_kcat.columns) | set(df_km.columns)
        for col in all_cols:
            if col not in df_kcat.columns:
                df_kcat[col] = None
            if col not in df_km.columns:
                df_km[col] = None
        
        # Ensure same column order
        col_order = sorted(all_cols)
        df_kcat = df_kcat[col_order]
        df_km = df_km[col_order]
        
        # Concatenate
        merged = pd.concat([df_kcat, df_km], ignore_index=True)
        
        return merged


# ============================================================
# DOCUMENT COMPOSER (SABIO-RK FORMAT - FINAL)
# ============================================================

class SABIORKFormatComposer:
    """Compose EnzymeML documents in SABIO-RK format"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
        self.doc_counter = 0
    
    def compose(self, group_data, uniprot_id, substrate, substrate_smiles):
        """Compose a single EnzymeML document for a protein-substrate group"""
        doc_id = f"SKiD_{uniprot_id}_{sanitize_id(substrate)}"
        
        # Get representative row for basic info
        first_row = group_data.iloc[0]
        
        # Extract basic info
        ec_number = first_row.get('EC_number')
        organism = first_row.get('Organism')
        sequence = first_row.get('Sequence') if 'Sequence' in first_row.index else None
        
        # Extract substrate additional info
        substrate_inchi = first_row.get('InChI') if 'InChI' in first_row.index else None
        substrate_inchikey = first_row.get('InChIKey') if 'InChIKey' in first_row.index else None
        substrate_external_ids = first_row.get('External_Identifiers') if 'External_Identifiers' in first_row.index else None
        
        # For variant info
        has_any_mutant = group_data['Mutant'].apply(
            lambda x: str(x).strip().lower() == 'yes' if pd.notna(x) else False
        ).any()
        
        if has_any_mutant:
            mutations = group_data[
                group_data['Mutant'].apply(
                    lambda x: str(x).strip().lower() == 'yes' if pd.notna(x) else False
                )
            ]['Mutation'].unique()
            variant_type = 'mutant'
            variant_description = '; '.join([str(m) for m in mutations if pd.notna(m)])
        else:
            variant_type = 'wildtype'
            variant_description = None
        
        # Generate enzyme name from organism
        enzyme_name = extract_enzyme_name_from_data(ec_number, organism, uniprot_id)
        
        # Collect and parse all references
        all_references = []
        for idx, row in group_data.iterrows():
            refs = parse_references(row.get('References'))
            if refs:
                for ref in refs:
                    # Avoid duplicates
                    if ref not in all_references:
                        all_references.append(ref)
        
        # Create document structure (SABIO-RK format with product placeholder)
        doc = {
            'name': f'SKiD Entry {doc_id}',
            'version': '2.0',
            'description': f'EnzymeML for SKiD entry {doc_id}',
            'created': datetime.now().isoformat(),
            'modified': None,
            'creators': [
                {
                    'given_name': 'SKiD',
                    'family_name': 'Database',
                    'mail': 'skid@database.org'
                }
            ],
            'vessels': None,
            'proteins': [
                {
                    'id': f'protein_{uniprot_id}',
                    'name': enzyme_name,  # Use organism name
                    'sequence': str(sequence) if sequence and pd.notna(sequence) else None,
                    'ecnumber': str(ec_number) if ec_number and pd.notna(ec_number) else None,
                    'organism': str(organism) if organism and pd.notna(organism) else None,
                    'organism_tax_id': None,
                    'synonymous_names': None,
                    'variant_type': variant_type,
                    'variant_description': variant_description,
                    'uniprotid': str(uniprot_id) if uniprot_id and pd.notna(uniprot_id) else None
                }
            ],
            'complexes': None,
            'small_molecules': [
                {
                    'id': f'substrate_{sanitize_id(substrate)}',
                    'name': substrate,
                    'canonical_smiles': substrate_smiles if substrate_smiles and pd.notna(substrate_smiles) else None,
                    'inchi': str(substrate_inchi) if substrate_inchi and pd.notna(substrate_inchi) else None,
                    'inchikey': str(substrate_inchikey) if substrate_inchikey and pd.notna(substrate_inchikey) else None,
                    'synonymous_names': None,
                    'KEGGid': None,
                    'CHEbIid': None,
                    'notes': str(substrate_external_ids) if substrate_external_ids and pd.notna(substrate_external_ids) else None
                },
                {
                    'id': 'product',  # Product as placeholder
                    'name': None,
                    'canonical_smiles': None,
                    'inchi': None,
                    'inchikey': None,
                    'synonymous_names': None,
                    'KEGGid': None,
                    'CHEbIid': None,
                    'notes': None
                },
                {
                    'id': 'activator',
                    'name': None,
                    'canonical_smiles': None,
                    'inchi': None,
                    'inchikey': None,
                    'synonymous_names': None,
                    'KEGGid': None,
                    'CHEbIid': None,
                    'notes': None
                },
                {
                    'id': 'inhibitor',
                    'name': None,
                    'canonical_smiles': None,
                    'inchi': None,
                    'inchikey': None,
                    'synonymous_names': None,
                    'KEGGid': None,
                    'CHEbIid': None,
                    'notes': None
                },
                {
                    'id': 'cofactor',
                    'name': None,
                    'canonical_smiles': None,
                    'inchi': None,
                    'inchikey': None,
                    'synonymous_names': None,
                    'KEGGid': None,
                    'CHEbIid': None,
                    'notes': None
                },
                {
                    'id': 'metal',
                    'name': None,
                    'canonical_smiles': None,
                    'inchi': None,
                    'inchikey': None,
                    'synonymous_names': None,
                    'KEGGid': None,
                    'CHEbIid': None,
                    'notes': None
                }
            ],
            'reactions': [
                {
                    'id': f'reaction_{uniprot_id}_{sanitize_id(substrate)}',
                    'name': f'Reaction catalyzed by enzyme on {substrate}',
                    'equation': None,
                    'smiles_equation': None,
                    'reversible': True,
                    'reactants': [
                        {
                            'species_id': f'substrate_{sanitize_id(substrate)}'
                        }
                    ],
                    'products': [
                        {
                            'species_id': 'product'
                        }
                    ],
                    'modifiers': [
                        {
                            'species_id': f'protein_{uniprot_id}',
                            'role': 'BIOCATALYST'
                        }
                    ],
                    'selectivity': None,
                    'Rheaid': None,
                    'KEGGid': None,
                    'notes': None
                }
            ],
            'measurements': [],
            'parameters': [],
            'buffer': {
                'buffer': None,
                'buffer_concentration': None
            },
            'references': all_references if all_references else None
        }
        
        # Deduplicate measurements by (pH, temperature)
        measurements_added = set()
        
        # Process all rows in the group
        for idx, row in group_data.iterrows():
            # Extract parameters
            kcat_val = row.get('kcat_value (1/s)')
            km_val = row.get('Km_value (mM)')
            ph = parse_ph(row.get('pH'))
            temperature = parse_temperature(row.get('T_oC'))
            mutation = row.get('Mutation')
            data_source = row.get('data_source', 'unknown')
            
            # Extract PMID from references
            pmid = extract_pmid_from_references(row.get('References'))
            
            # Create measurement (deduplicated)
            measurement_key = (ph, temperature)
            if measurement_key not in measurements_added and (ph or temperature):
                measurement_id = f'measurement_{len(doc["measurements"])}'
                measurement_entry = {
                    'id': measurement_id,
                    'name': f'Experimental conditions (pH={ph}, T={temperature}°C)',
                    'species_data': None,
                    'ph_opt': ph,
                    'temperature_opt': temperature,
                    'temperature_unit': '°C' if temperature else None
                }
                doc['measurements'].append(measurement_entry)
                measurements_added.add(measurement_key)
            
            # Add kcat parameter
            if kcat_val and pd.notna(kcat_val):
                try:
                    kcat_float = float(kcat_val)
                    param_id = f'{uniprot_id}_{sanitize_id(substrate)}_kcat_{idx}'
                    
                    if mutation and pd.notna(mutation) and str(mutation).strip() not in ['-----', '']:
                        param_id = f'{uniprot_id}_{sanitize_id(substrate)}_{sanitize_id(str(mutation))}_kcat_{idx}'
                    
                    notes_parts = [f'Source: {data_source}']
                    if mutation and pd.notna(mutation) and str(mutation).strip() not in ['-----', '']:
                        notes_parts.insert(0, f'Mutation: {mutation}')
                    
                    doc['parameters'].append({
                        'id': param_id,
                        'name': 'kcat',
                        'symbol': None,
                        'value': None,
                        'unit': '1/s',
                        'initial_value': kcat_float,
                        'upper_bound': None,
                        'lower_bound': None,
                        'associated_species': f'substrate_{sanitize_id(substrate)}',
                        'original_name': None,
                        'parameter.type': 'KINETIC_CONSTANT',
                        'pH': ph,
                        'temperature': temperature,
                        'temperature_unit': '°C' if temperature else None,
                        'notes': '; '.join(notes_parts),
                        'pubmed_id': pmid  # Add PMID to parameter
                    })
                except:
                    pass
            
            # Add Km parameter
            if km_val and pd.notna(km_val):
                try:
                    km_float = float(km_val)
                    param_id = f'{uniprot_id}_{sanitize_id(substrate)}_Km_{idx}'
                    
                    if mutation and pd.notna(mutation) and str(mutation).strip() not in ['-----', '']:
                        param_id = f'{uniprot_id}_{sanitize_id(substrate)}_{sanitize_id(str(mutation))}_Km_{idx}'
                    
                    notes_parts = [f'Source: {data_source}']
                    if mutation and pd.notna(mutation) and str(mutation).strip() not in ['-----', '']:
                        notes_parts.insert(0, f'Mutation: {mutation}')
                    
                    doc['parameters'].append({
                        'id': param_id,
                        'name': 'Km',
                        'symbol': None,
                        'value': None,
                        'unit': 'mM',  # Keep mM
                        'initial_value': km_float,  # Keep original value
                        'upper_bound': None,
                        'lower_bound': None,
                        'associated_species': f'substrate_{sanitize_id(substrate)}',
                        'original_name': None,
                        'parameter.type': 'KINETIC_CONSTANT',
                        'pH': ph,
                        'temperature': temperature,
                        'temperature_unit': '°C' if temperature else None,
                        'notes': '; '.join(notes_parts),
                        'pubmed_id': pmid  # Add PMID to parameter
                    })
                except:
                    pass
        
        # If no measurements were added, set to None
        if not doc['measurements']:
            doc['measurements'] = None
        
        # If no parameters were added, set to None
        if not doc['parameters']:
            doc['parameters'] = None
        
        self.doc_counter += 1
        return doc


# ============================================================
# YAML GENERATOR
# ============================================================

class YAMLGenerator:
    """Generate YAML files with grouping by protein+substrate"""
    
    def __init__(self, schema_classes: dict):
        self.schema_classes = schema_classes
        self.composer = SABIORKFormatComposer(schema_classes)
    
    def generate(self, df: pd.DataFrame, output_dir: str, format: str = 'yaml', 
                 limit: int = 0, output_mode: str = 'both'):
        """Generate YAML files - ONE FILE PER (UniProt, Substrate, Mutation) combination"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"\nGenerating files in SABIO-RK format...")
        print(f"  Total records: {len(df)}")
        print(f"  Output mode: {output_mode}")
        print(f"  Grouping by: (UniProt_ID, Substrate, Mutation)")
        
        # Statistics
        has_kcat = df['kcat_value (1/s)'].notna().sum()
        has_km = df['Km_value (mM)'].notna().sum()
        has_both = df[(df['kcat_value (1/s)'].notna()) & (df['Km_value (mM)'].notna())].shape[0]
        
        print(f"\n  Records with kcat: {has_kcat}")
        print(f"  Records with Km: {has_km}")
        print(f"  Records with both: {has_both}")
        
        # Group by (UniProt_ID, Substrate, Mutation)
        # Handle wildtype: if Mutant != 'yes', treat Mutation as 'WT'
        df['Mutation_key'] = df.apply(
            lambda row: row['Mutation'] if str(row.get('Mutant', '')).strip().lower() == 'yes' 
            else 'WT', axis=1
        )
        
        grouped = df.groupby(['UniProt_ID', 'Substrate', 'Mutation_key'], dropna=False)
        print(f"  Unique (UniProt, Substrate, Mutation) combinations: {len(grouped)}")
        
        count = 0
        file_counter = {}  # Track filename usage to avoid conflicts
        
        # Process each group
        for (uniprot_id, substrate, mutation_key), group in grouped:
            if limit > 0 and count >= limit:
                break
            
            # Skip if no valid uniprot_id or substrate
            if pd.isna(uniprot_id) or pd.isna(substrate):
                continue
            
            # Get substrate SMILES
            substrate_smiles = group['Smiles'].dropna().iloc[0] if not group['Smiles'].dropna().empty else None
            
            # Get EC number and organism for filename
            ec_number = group['EC_number'].dropna().iloc[0] if not group['EC_number'].dropna().empty else None
            organism = group['Organism'].dropna().iloc[0] if not group['Organism'].dropna().empty else None
            
            # Get mutation info
            mutation = mutation_key if mutation_key != 'WT' else None
            
            # Compose document with all records in this group
            doc = self.composer.compose(group, uniprot_id, substrate, substrate_smiles)
            
            # Complete to full schema
            completed = complete_doc_to_schema(doc, self.schema_classes)
            
            # Clean unwanted fields
            completed = clean_unwanted_fields(completed)
            
            # Generate filename: SKID_EC_organism_uniprot_(mutation)_substrate.yaml
            filename_parts = ['SKID']
            
            # Add EC number (shortened)
            if ec_number and pd.notna(ec_number):
                ec_short = sanitize_id(str(ec_number))[:15]  # Limit EC to 15 chars
                filename_parts.append(f"EC{ec_short}")
            
            # Add organism (genus only, first 10 chars)
            if organism and pd.notna(organism):
                org_first = str(organism).split()[0][:10]
                filename_parts.append(sanitize_id(org_first))
            
            # Add UniProt ID (max 10 chars)
            if uniprot_id and pd.notna(uniprot_id):
                uniprot_short = sanitize_id(str(uniprot_id))[:10]
                filename_parts.append(uniprot_short)
            
            # Add mutation if present (single mutation, max 15 chars)
            # For wildtype, do NOT add 'WT' to filename
            if mutation and pd.notna(mutation):
                mutation_short = sanitize_id(str(mutation))[:15]
                filename_parts.append(mutation_short)
            
            # Add substrate (max 15 chars)
            substrate_short = sanitize_id(str(substrate))[:15]
            filename_parts.append(substrate_short)
            
            # Create base filename
            base_filename = "_".join(filename_parts)
            
            # Ensure filename is unique and not too long
            # Windows has 260 char path limit, reserve some for directory
            max_filename_len = 150
            if len(base_filename) > max_filename_len:
                # Use a hash for very long names
                import hashlib
                hash_suffix = hashlib.md5(base_filename.encode()).hexdigest()[:8]
                base_filename = base_filename[:max_filename_len-9] + "_" + hash_suffix
            
            # Add counter if filename exists
            if base_filename not in file_counter:
                file_counter[base_filename] = 0
            file_counter[base_filename] += 1
            
            if file_counter[base_filename] > 1:
                filename = f"{base_filename}_{file_counter[base_filename]}"
            else:
                filename = base_filename
            
            # Write file
            try:
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
                if count % 1000 == 0:
                    print(f"  Generated {count} files...")
            
            except Exception as e:
                print(f"  Warning: Failed to write file {filename}: {e}")
                continue
        
        print(f"\n✅ Complete! Generated {count} files")
        print(f"  Output: {output_path.resolve()}")
        print(f"  Format: SABIO-RK compatible (One file per UniProt+Substrate+Mutation)")


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
        description="SKiD Main Dataset to EnzymeML - SABIO-RK Format (Final)"
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
    
    # Load data from all sheets
    print(f"\nLoading data: {args.input}")
    
    # Load Unique_substrates for InChI and External_identifiers
    print("  Reading Unique_substrates sheet...")
    df_substrates = pd.read_excel(args.input, sheet_name='Unique_substrates')
    print(f"  Loaded {len(df_substrates)} substrate records")
    
    # Create substrate info mapping: {substrate_name: {InChI, InChIKey, External_identifiers}}
    substrate_info = {}
    for _, row in df_substrates.iterrows():
        ligand_name = row.get('Ligand_name')
        if pd.notna(ligand_name):
            substrate_info[str(ligand_name).strip()] = {
                'InChI': row.get('InChI') if pd.notna(row.get('InChI')) else None,
                'InChIKey': row.get('InChIKey') if pd.notna(row.get('InChIKey')) else None,
                'External_identifiers': row.get('External_identifiers') if pd.notna(row.get('External_identifiers')) else None
            }
    print(f"  Created substrate info mapping for {len(substrate_info)} substrates")
    
    print("  Reading kcat_dataset sheet...")
    df_kcat = pd.read_excel(args.input, sheet_name='kcat_dataset')
    print(f"  Loaded {len(df_kcat)} kcat records")
    
    print("  Reading Km_dataset sheet...")
    df_km = pd.read_excel(args.input, sheet_name='Km_dataset')
    print(f"  Loaded {len(df_km)} Km records")
    
    # Merge substrate info into kcat and Km datasets
    print("\n  Merging substrate information (InChI, InChIKey, External_identifiers)...")
    for df in [df_kcat, df_km]:
        df['InChI'] = df['Substrate'].apply(
            lambda x: substrate_info.get(str(x).strip(), {}).get('InChI') if pd.notna(x) else None
        )
        df['InChIKey'] = df['Substrate'].apply(
            lambda x: substrate_info.get(str(x).strip(), {}).get('InChIKey') if pd.notna(x) else None
        )
        df['External_Identifiers'] = df['Substrate'].apply(
            lambda x: substrate_info.get(str(x).strip(), {}).get('External_identifiers') if pd.notna(x) else None
        )
    
    # Merge kcat and Km data
    print("\n  Merging kcat and Km data (preserving all records)...")
    merger = SKiDDataMerger(df_kcat, df_km)
    merged_df = merger.merge_data()
    print(f"  Combined dataset: {len(merged_df)} total records")
    
    # Generate files
    generator = YAMLGenerator(schema_classes)
    generator.generate(merged_df, args.output_dir, args.format, args.limit, args.output_mode)


if __name__ == '__main__':
    main()
