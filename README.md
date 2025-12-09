# EnzyTransfer

This script is a collection for **standardizing heterogeneous enzyme kinetics datasets** (SABIO-RK, BRENDA, Rhea, SKiD, RetroBioCat, …) into a **unified EnzymeML-style YAML/JSON format**, driven by a shared LinkML schema.

---

## Repository layout

Recommended layout in the project root:

```text
.
├── data/                      # Raw input data (NOT tracked in git)
│   ├── Sabio_rk/
│   ├── BRENDA/
│   ├── Rhea/
│   ├── SKID/
│   └── RetroBioCat/
├── data_sources/              # Dataset-specific pipelines
│   ├── BRENDA/
│   │   ├── step0_llm.py       # LLM-based header standardizer (generic)
│   │   └── step1_join.py      # BRENDA → unified EnzymeML
│   ├── sabio_rk/
│   │   ├── step0_header.py    # (placeholder) header helper
│   │   └── step1_extract.py   # SABIO-RK → unified EnzymeML
│   ├── Rhea/
│   │   └── step1_join_modified.py   # Rhea → unified EnzymeML
│   ├── SKID/
│   │   └── step1_skid.py      # SKiD → unified EnzymeML
│   ├── RetroBioCat/
│   │   ├── step0_pub.py       # Filter YAMLs by target PMIDs
│   │   └── step1_test.py      # RetroBioCat → unified EnzymeML
│   ├── doi_pub/
│   │   └── step0.py           # DOI → PMID helper
│   ├── mutation/
│   │   ├── step0.py           # Fill sequences & categorize WT / mutant
│   │   ├── step1.py           # Fetch/repair sequences from UniProt
│   │   ├── step1_2.py         # Pattern-based retry for mutations
│   │   └── step2_multi.py     # Resolve records with multiple UniProt IDs
│   └── merge/
│       └── merge_sequences.py # Sequence-based deduplication & merge
├── schemas/
│   └── enzymeml-v2-extended.yaml   # Central LinkML schema
├── output/                   # (Recommended) EnzymeML YAML/JSON from each source
│   ├── sabio_rk/
│   ├── brenda/
│   ├── rhea/
│   ├── skid/
│   └── retrobiocat/
├── .gitignore
└── .gitattributes
