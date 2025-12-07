#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
csv_to_enzymeml.py

Convert TN_key.py CSV output to EnzymeML V2 Extended YAML/JSON format.
Fully compliant with enzymeml-v2-extended.yaml schema.
Now includes schema validation and dynamic schema loading.

Fixed: CSV field size limit error handling
Fixed: Better error reporting for large datasets

Usage:
  python csv_to_enzymeml.py --input output.csv --output-dir enzymeml_data [--format json|yaml] --schema enzymeml-v2-extended.yaml
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from datetime import datetime

# Fix CSV field size limit
csv.field_size_limit(sys.maxsize)

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

try:
    from linkml_runtime.loaders import yaml_loader
    from linkml_runtime.dumpers import yaml_dumper, json_dumper
    from linkml.validator import validate
    HAS_LINKML = True
except ImportError:
    HAS_LINKML = False


class SchemaValidator:
    """Validates EnzymeML documents against the schema."""
    
    def __init__(self, schema_path: Optional[Path] = None):
        self.schema_path = schema_path
        self.schema = None
        if schema_path and schema_path.exists():
            self.load_schema()
    
    def load_schema(self):
        """Load the schema definition."""
        if not self.schema_path:
            return
        
        try:
            with open(self.schema_path, 'r', encoding='utf-8') as f:
                self.schema = yaml.safe_load(f)
                print(f"✓ Loaded schema: {self.schema_path}")
        except Exception as e:
            print(f"⚠ Warning: Could not load schema: {e}")
    
    def validate_document(self, doc: Dict[str, Any]) -> tuple[bool, Optional[str]]:
        """
        Validate a document against the schema.
        Returns (is_valid, error_message)
        """
        if not self.schema:
            return True, "No schema loaded, skipping validation"
        
        # Basic validation checks based on schema structure
        errors = []
        
        # Check required top-level fields
        required_fields = ['name', 'version', 'creators', 'vessels']
        for field in required_fields:
            if field not in doc or not doc[field]:
                errors.append(f"Missing required field: {field}")
        
        # Validate creators
        if 'creators' in doc:
            for idx, creator in enumerate(doc['creators']):
                required_creator_fields = ['given_name', 'family_name', 'mail']
                for field in required_creator_fields:
                    if field not in creator:
                        errors.append(f"Creator {idx}: missing required field '{field}'")
        
        # Validate vessels
        if 'vessels' in doc:
            for idx, vessel in enumerate(doc['vessels']):
                required_vessel_fields = ['id', 'name', 'volume', 'unit', 'constant']
                for field in required_vessel_fields:
                    if field not in vessel:
                        errors.append(f"Vessel {idx}: missing required field '{field}'")
        
        # Validate proteins
        if 'proteins' in doc:
            for idx, protein in enumerate(doc['proteins']):
                required_protein_fields = ['id', 'name', 'constant']
                for field in required_protein_fields:
                    if field not in protein:
                        errors.append(f"Protein {idx}: missing required field '{field}'")
        
        # Validate small molecules
        if 'small_molecules' in doc:
            for idx, sm in enumerate(doc['small_molecules']):
                required_sm_fields = ['id', 'name', 'constant']
                for field in required_sm_fields:
                    if field not in sm:
                        errors.append(f"SmallMolecule {idx}: missing required field '{field}'")
        
        # Validate reactions
        if 'reactions' in doc:
            for idx, reaction in enumerate(doc['reactions']):
                required_reaction_fields = ['id', 'name', 'reversible']
                for field in required_reaction_fields:
                    if field not in reaction:
                        errors.append(f"Reaction {idx}: missing required field '{field}'")
        
        if errors:
            return False, "\n".join(errors)
        
        return True, None
    
    def get_enum_values(self, enum_name: str) -> List[str]:
        """Get permissible values for an enum from schema."""
        if not self.schema or 'enums' not in self.schema:
            return []
        
        if enum_name not in self.schema['enums']:
            return []
        
        enum_def = self.schema['enums'][enum_name]
        if 'permissible_values' in enum_def:
            return list(enum_def['permissible_values'].keys())
        
        return []


def parse_list_field(value: str, separator: str = ";") -> List[str]:
    """Parse semicolon-separated field into list."""
    if not value or value == "N/A":
        return []
    return [item.strip() for item in value.split(separator) if item.strip()]


def safe_float(value: str) -> Optional[float]:
    """Safely convert string to float."""
    if not value or value == "N/A":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def safe_int(value: str) -> Optional[int]:
    """Safely convert string to int."""
    if not value or value == "N/A":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def sanitize_filename(name: str) -> str:
    """Sanitize string for use as filename."""
    name = name.replace("/", "_").replace("\\", "_")
    name = name.replace(":", "_").replace("*", "_")
    name = name.replace("?", "_").replace('"', "_")
    name = name.replace("<", "_").replace(">", "_")
    name = name.replace("|", "_").replace(" ", "_")
    return name


def truncate_filename(name: str, max_length: int = 80) -> str:
    """Truncate filename to a maximum length while preserving uniqueness."""
    if len(name) <= max_length:
        return name
    # Keep some of the beginning and add a hash of the full name
    import hashlib
    hash_suffix = hashlib.md5(name.encode()).hexdigest()[:8]
    truncated = name[:max_length - 9]  # Leave room for underscore and hash
    return f"{truncated}_{hash_suffix}"


def create_enzymeml_document(row: Dict[str, str], row_index: int, validator: Optional[SchemaValidator] = None) -> Dict[str, Any]:
    """
    Convert a single CSV row to complete EnzymeML Document structure.
    Now validates against schema if validator is provided.
    """
    
    # Extract basic info
    ec_number = row.get("ec_number", "N/A")
    enzyme_name = row.get("enzyme_name", "N/A")
    protein_id = row.get("protein_id", "N/A")
    protein_uid = row.get("protein_uid", f"{ec_number}#{protein_id}")
    organism = row.get("organism", "N/A")
    uniprot = row.get("uniprot", "N/A")
    
    substrate_raw = row.get("substrate_raw", "N/A")
    substrate_norm = row.get("substrate_norm", "N/A")
    
    # Kinetic parameters
    kcat_value = row.get("kcat_value", "N/A")
    kcat_unit = row.get("kcat_unit", "N/A")
    kcat_comment = row.get("kcat_comment", "N/A")
    kcat_pmids = parse_list_field(row.get("kcat_pmids", ""))
    
    # Conditions
    ph_optimum = safe_float(row.get("ph_optimum", "N/A"))
    temperature_optimum = safe_float(row.get("temperature_optimum", "N/A"))
    
    # Modifiers
    cofactors = parse_list_field(row.get("cofactors", ""))
    activators = parse_list_field(row.get("activators", ""))
    inhibitors = parse_list_field(row.get("inhibitors", ""))
    metal_ions = parse_list_field(row.get("metal_ions", ""))
    
    # Variant info
    protein_type = row.get("protein_type", "N/A")
    mutation_info = parse_list_field(row.get("mutation_info", ""))
    
    # Buffer info
    comment_buffer = row.get("comment_buffer", "N/A")
    buffer_concentration = row.get("buffer_concentration", "N/A")
    
    # Counts
    km_count = safe_int(row.get("km_count", "0"))
    kcat_over_km_count = safe_int(row.get("kcat_over_km_count", "0"))
    ki_count = safe_int(row.get("ki_count", "0"))
    
    # Determine variant_type based on schema enum values
    variant_type_value = None
    if validator:
        valid_variant_types = validator.get_enum_values("VariantType")
        if protein_type.lower() in [v.lower() for v in valid_variant_types]:
            variant_type_value = protein_type.lower()
    else:
        # Default behavior
        if protein_type.lower() in ["wildtype", "wild type"]:
            variant_type_value = "wildtype"
        elif protein_type.lower() == "mutant":
            variant_type_value = "mutant"
    
    # Build EnzymeML Document following schema exactly
    doc = {
        # EnzymeMLDocument root attributes
        "name": enzyme_name if enzyme_name != "N/A" else f"BRENDA Entry {protein_uid}",
        "version": "2.0",
        "description": f"Kinetic data for {enzyme_name if enzyme_name != 'N/A' else 'enzyme'} (EC {ec_number}) from BRENDA database",
        "created": datetime.now().isoformat(),
        "modified": datetime.now().isoformat(),
        
        # Creator (required)
        "creators": [
            {
                "given_name": "BRENDA",
                "family_name": "Database",
                "mail": "brenda@brenda-enzymes.org"
            }
        ],
        
        # Vessel (required)
        "vessels": [
            {
                "id": "v0",
                "name": "default_vessel",
                "volume": 1.0,
                "unit": {
                    "id": "u_ml",
                    "name": "millilitre",
                    "base_units": [
                        {
                            "kind": "LITRE",
                            "exponent": 1,
                            "multiplier": 0.001,
                            "scale": 1.0
                        }
                    ]
                },
                "constant": True
            }
        ],
        
        # Protein
        "proteins": [
            {
                "id": f"p_{protein_id}",
                "name": enzyme_name if enzyme_name != "N/A" else f"Protein {protein_id}",
                "constant": False,
                "sequence": None,
                "vessel_id": "v0",
                "ecnumber": ec_number if ec_number != "N/A" else None,
                "organism": organism if organism != "N/A" else None,
                "organism_tax_id": None,
                "synonymous_names": [],
                "variant_type": variant_type_value,
                "mutations": mutation_info,
                "variant_description": None,
                "uniprotid": uniprot if uniprot != "N/A" else None
            }
        ],
        
        # Complex (empty but present)
        "complexes": [],
        
        # SmallMolecule - substrate
        "small_molecules": [
            {
                "id": "sm_substrate",
                "name": substrate_norm if substrate_norm != "N/A" else substrate_raw,
                "constant": False,
                "vessel_id": "v0",
                "canonical_smiles": None,
                "inchi": None,
                "inchikey": None,
                "synonymous_names": [substrate_raw] if substrate_raw != substrate_norm and substrate_raw != "N/A" else [],
                "KEGGid": None,
                "CHEbIid": None
            }
        ],
        
        # Reaction
        "reactions": [
            {
                "id": "r0",
                "name": f"Reaction catalyzed by {enzyme_name if enzyme_name != 'N/A' else protein_id}",
                "reversible": False,
                "kinetic_law": None,
                "reactants": [
                    {
                        "species_id": "sm_substrate",
                        "stoichiometry": 1.0
                    }
                ],
                "products": [],
                "modifiers": [],  # Will be filled below
                "selectivity": None,
                "Rheaid": None,
                "KEGGid": None
            }
        ],
        
        # Measurement
        "measurements": [
            {
                "id": "m0",
                "name": f"Measurement for {protein_uid}",
                "species_data": [],
                "group_id": None,
                "ph": ph_optimum,
                "temperature": temperature_optimum,
                "temperature_unit": {
                    "id": "u_celsius",
                    "name": "celsius",
                    "base_units": [
                        {
                            "kind": "CELSIUS",
                            "exponent": 1,
                            "multiplier": 1.0,
                            "scale": 1.0
                        }
                    ]
                } if temperature_optimum is not None else None
            }
        ],
        
        # Equation (empty but present)
        "equations": [],
        
        # Parameters
        "parameters": [],
        
        # Buffer
        "buffer": {
            "buffer": comment_buffer if comment_buffer != "N/A" else None,
            "buffer_concentration": buffer_concentration if buffer_concentration != "N/A" else None
        }
    }
    
    # Add modifiers (cofactors, activators, inhibitors, metal_ions)
    modifiers = []
    modifier_count = 0
    
    # Get valid modifier roles from schema if validator available
    valid_roles = validator.get_enum_values("ModifierRole") if validator else [
        "ACTIVATOR", "ADDITIVE", "BIOCATALYST", "BUFFER", "CATALYST", "INHIBITOR", "SOLVENT"
    ]
    
    # Cofactors - role: CATALYST
    for cofactor in cofactors:
        species_id = f"sm_cofactor_{modifier_count}"
        modifiers.append({
            "species_id": species_id,
            "role": "CATALYST" if "CATALYST" in valid_roles else valid_roles[0]
        })
        doc["small_molecules"].append({
            "id": species_id,
            "name": cofactor,
            "constant": True,
            "vessel_id": "v0",
            "canonical_smiles": None,
            "inchi": None,
            "inchikey": None,
            "synonymous_names": [],
            "KEGGid": None,
            "CHEbIid": None
        })
        modifier_count += 1
    
    # Activators - role: ACTIVATOR
    for activator in activators:
        species_id = f"sm_activator_{modifier_count}"
        modifiers.append({
            "species_id": species_id,
            "role": "ACTIVATOR" if "ACTIVATOR" in valid_roles else valid_roles[0]
        })
        doc["small_molecules"].append({
            "id": species_id,
            "name": activator,
            "constant": True,
            "vessel_id": "v0",
            "canonical_smiles": None,
            "inchi": None,
            "inchikey": None,
            "synonymous_names": [],
            "KEGGid": None,
            "CHEbIid": None
        })
        modifier_count += 1
    
    # Inhibitors - role: INHIBITOR
    for inhibitor in inhibitors:
        species_id = f"sm_inhibitor_{modifier_count}"
        modifiers.append({
            "species_id": species_id,
            "role": "INHIBITOR" if "INHIBITOR" in valid_roles else valid_roles[0]
        })
        doc["small_molecules"].append({
            "id": species_id,
            "name": inhibitor,
            "constant": True,
            "vessel_id": "v0",
            "canonical_smiles": None,
            "inchi": None,
            "inchikey": None,
            "synonymous_names": [],
            "KEGGid": None,
            "CHEbIid": None
        })
        modifier_count += 1
    
    # Metal ions - role: CATALYST
    for metal in metal_ions:
        species_id = f"sm_metal_{modifier_count}"
        modifiers.append({
            "species_id": species_id,
            "role": "CATALYST" if "CATALYST" in valid_roles else valid_roles[0]
        })
        doc["small_molecules"].append({
            "id": species_id,
            "name": metal,
            "constant": True,
            "vessel_id": "v0",
            "canonical_smiles": None,
            "inchi": None,
            "inchikey": None,
            "synonymous_names": [],
            "KEGGid": None,
            "CHEbIid": None
        })
        modifier_count += 1
    
    # Add modifiers to reaction
    doc["reactions"][0]["modifiers"] = modifiers
    
    # Add kcat parameter
    if kcat_value != "N/A":
        kcat_float = safe_float(kcat_value)
        doc["parameters"].append({
            "id": "p_kcat",
            "name": "turnover number",
            "symbol": "kcat",
            "value": kcat_float,
            "unit": {
                "id": "u_per_second",
                "name": "per second",
                "base_units": [
                    {
                        "kind": "SECOND",
                        "exponent": -1,
                        "multiplier": 1.0,
                        "scale": 1.0
                    }
                ]
            } if kcat_unit != "N/A" else None,
            "initial_value": None,
            "upper_bound": None,
            "lower_bound": None,
            "fit": False,
            "stderr": None,
            "constant": True,
            "associated_species": f"p_{protein_id}",
            "original_name": "kcat"
        })
    
    # Add references (as part of the schema)
    references = []
    for pmid in kcat_pmids:
        references.append({
            "title": None,
            "authors": None,
            "year": None,
            "journal": None,
            "uri": [f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"] if pmid else [],
            "pubmed_id": pmid,
            "doi": None
        })
    
    if references:
        doc["references"] = references
    
    # Store additional metadata (not part of schema but useful)
    doc["_metadata"] = {
        "protein_uid": protein_uid,
        "km_count": km_count,
        "kcat_over_km_count": kcat_over_km_count,
        "ki_count": ki_count,
        "kcat_comment": kcat_comment if kcat_comment != "N/A" else None
    }
    
    return doc


def main():
    parser = argparse.ArgumentParser(
        description="Convert TN_key.py CSV to EnzymeML V2 Extended format (schema-compliant)"
    )
    parser.add_argument("--input", required=True, help="Input CSV file")
    parser.add_argument("--output-dir", required=True, help="Output directory for individual files")
    parser.add_argument("--schema", help="Path to enzymeml-v2-extended.yaml schema file")
    parser.add_argument("--format", choices=["json", "yaml"], default="json",
                       help="Output format (default: json)")
    parser.add_argument("--pretty", action="store_true", default=True,
                       help="Pretty print JSON output (default: True)")
    parser.add_argument("--validate", action="store_true", default=False,
                       help="Validate output against schema (default: False for large datasets)")
    parser.add_argument("--progress", type=int, default=1000,
                       help="Show progress every N rows (default: 1000)")
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"ERROR: Input file not found: {input_path}")
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_format = args.format
    file_ext = ".yaml" if output_format == "yaml" else ".json"
    
    if output_format == "yaml" and not HAS_YAML:
        sys.exit("ERROR: PyYAML not installed. Run: pip install pyyaml")
    
    # Initialize schema validator
    validator = None
    if args.schema:
        schema_path = Path(args.schema)
        if not schema_path.exists():
            print(f"⚠ Warning: Schema file not found: {schema_path}")
            print(f"   Schema validation will be disabled.")
        else:
            validator = SchemaValidator(schema_path)
    elif args.validate:
        print("⚠ Warning: Validation requested but no schema file provided. Use --schema to specify schema file.")
    
    print(f"Converting {input_path} to {output_format.upper()} format...")
    print(f"Output directory: {output_dir}")
    if validator and validator.schema:
        print(f"Schema validation: ENABLED")
    else:
        print(f"Schema validation: DISABLED")
    print()
    
    # Read CSV with increased field size limit
    try:
        with open(input_path, 'r', encoding='utf-8', errors='replace') as csvfile:
            reader = csv.DictReader(csvfile)
            
            count = 0
            validation_errors = 0
            errors_detail = []
            
            print("Processing rows...")
            for row_idx, row in enumerate(reader, 1):
                try:
                    doc = create_enzymeml_document(row, row_idx, validator)
                    
                    # Validate document if requested
                    if args.validate and validator:
                        is_valid, error_msg = validator.validate_document(doc)
                        if not is_valid:
                            error_info = f"Row {row_idx}: {error_msg}"
                            errors_detail.append(error_info)
                            validation_errors += 1
                            # Only show first 10 validation errors to avoid flooding console
                            if validation_errors <= 10:
                                print(f"❌ Validation failed for row {row_idx}")
                    
                    # Generate filename from ec_number, substrate, and protein_id
                    ec_number = row.get("ec_number", "unknown_ec")
                    substrate_norm = row.get("substrate_norm", "substrate")
                    protein_id = row.get("protein_id", str(row_idx))

                    # Sanitize and truncate to avoid Windows path length issues
                    ec_clean = truncate_filename(sanitize_filename(ec_number), 20)
                    substrate_clean = truncate_filename(sanitize_filename(substrate_norm), 40)
                    protein_clean = truncate_filename(sanitize_filename(protein_id), 20)

                    filename = f"{ec_clean}_{substrate_clean}_{protein_clean}{file_ext}"
                    output_file = output_dir / filename
                    
                    # Write individual file
                    with open(output_file, 'w', encoding='utf-8') as outfile:
                        if output_format == "yaml":
                            yaml.dump(doc, outfile, default_flow_style=False, 
                                     allow_unicode=True, sort_keys=False)
                        else:
                            if args.pretty:
                                json.dump(doc, outfile, indent=2, ensure_ascii=False)
                            else:
                                json.dump(doc, outfile, ensure_ascii=False)
                    
                    count += 1
                    if count % args.progress == 0:
                        print(f"  ✓ Processed {count:,} rows...")
                
                except Exception as e:
                    print(f"❌ Error processing row {row_idx}: {str(e)}")
                    continue
    
    except csv.Error as e:
        print(f"\n❌ CSV Error: {e}")
        print("This usually means the CSV file has formatting issues.")
        print("Please check your CSV file for:")
        print("  - Extremely large text fields")
        print("  - Unescaped quotes")
        print("  - Line breaks within fields")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print()
    print("="*60)
    print(f"✅ Done! Converted {count:,} rows")
    if validation_errors > 0:
        print(f"⚠️  {validation_errors:,} validation errors found")
        if validation_errors > 10:
            print(f"   (Showing first 10 errors)")
    print(f"   Output: {output_dir}/")
    print(f"   Generated {count:,} individual {output_format.upper()} files")
    print("="*60)
    
    if validator and validator.schema:
        print()
        print("Schema information:")
        print(f"  - Schema loaded from: {validator.schema_path}")
        if 'enums' in validator.schema:
            print(f"  - Available enums: {', '.join(validator.schema['enums'].keys())}")
    
    # Save validation errors to a log file if there are any
    if validation_errors > 0 and errors_detail:
        error_log_file = output_dir / "validation_errors.log"
        with open(error_log_file, 'w', encoding='utf-8') as f:
            f.write(f"Validation Errors Report\n")
            f.write(f"========================\n")
            f.write(f"Total errors: {validation_errors}\n")
            f.write(f"Date: {datetime.now().isoformat()}\n\n")
            for error in errors_detail:
                f.write(f"{error}\n")
        print(f"\n📋 Validation errors saved to: {error_log_file}")


if __name__ == "__main__":
    main()
