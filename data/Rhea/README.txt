================================================================================
Rhea TSV
================================================================================

This directory contains various files in tab-separated values (TSV) format that
can be easily imported into relational databases and spreadsheet applications.

1. Cross-references to other databases

Files named rhea2<db>.tsv contain cross-references to other databases, where
<db> is the name of the cross-referenced database. They all contain 4 fields:

- RHEA_ID:   the unique identifier of the Rhea reaction.
- DIRECTION: the direction of the Rhea reaction (UN, LR, RL or BI).
- MASTER_ID: the unique identifier of the corresponding Rhea reaction with
             undefined direction (UN).
- ID:        the unique identifier of the cross-referenced record.

The file rhea2xrefs.tsv contains all cross-references, except those to
UniProtKB, with this additional field:

- DB:        the name of the cross-referenced database.

2. Small molecule information (ChEBI)

- chebiId_name.tsv
  List of ChEBI identifiers with the name that is used in Rhea.

- chebi_pH7_3_mapping.tsv
  Mapping of ChEBI identifiers to the major microspecies at pH 7.3.
  Fields: 1. ChEBI ID, 2. ChEBI ID of major microspecies at pH 7.3,
  3. origin of mapping ('computation' or 'curation').

3. SMILES

Note: These files are currently experimental and may change in the future
without prior notice.

- rhea-chebi-smiles.tsv
  Canonical SMILES for the subset of ChEBI used in Rhea.
  Computed with RDKit (https://www.rdkit.org/) using the ChEBI Molfile as input.
  Fields: 1. ChEBI ID, 2. SMILES

- rhea-reaction-smiles.tsv
  Reaction SMILES for Rhea directed reactions (left-to-right and right-to-left).
  Computed with RDKit (https://www.rdkit.org/) using the Rhea RXN file as input.
  Fields: 1. Rhea ID, 2. Reaction SMILES

4. Miscellaneous

- rhea-directions.tsv
  Mapping of undirected reactions with their corresponding left-to-right,
  right-to-left and bidirectional reactions.

- rhea-relationships.tsv
  Hierarchical reaction classification.

- rhea-obsoletes.tsv
  List of obsolete reactions.

- rhea-tsv.tar.gz
  Archive of all TSV files, except rhea2uniprot_trembl.tsv.gz
