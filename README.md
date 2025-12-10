# EnzyTransfer

**EnzyTransfer** is a collection of Python scripts designed to standardize heterogeneous enzyme kinetics datasets (including SABIO-RK, BRENDA, Rhea, SKiD, and RetroBioCat) into a **unified EnzymeML-style YAML/JSON format**.

## 📖 Overview

The core philosophy of EnzyTransfer is to decouple data sources from data structure:

* **Raw Data**: Stored in `data/`, organized by source folders.
* **Processing**: Source-specific pipelines reside in `data_sources/`.
* **Schema**: All converters reference a shared schema (e.g., `schemas/enzymeml-v2-extended.yaml`) to produce a **compatible EnzymeML-like structure**.
* **Result**: Downstream tools can treat all data sources uniformly, regardless of their origin.

---

## ✨ Features

* **Multi-source Standardization**: Converters for SABIO-RK, BRENDA, Rhea, SKiD, and RetroBioCat.
* **Schema-driven Design**: Uses a central extended EnzymeML schema for proteins, small molecules, reactions, kinetic parameters, and measurements.
* **Mutation & Sequence Enrichment**: Tools to fill protein sequences, parse mutation strings, fetch UniProt sequences, and classify wildtype vs. mutant records.
* **Reference Utilities**: Includes DOI-to-PMID converters and RetroBioCat YAML filtering tools.
* **Merge & Deduplicate**: Parallel sequence-based comparison and merging across different data sources.

---

## 📂 Repository Structure

```text
.
├── data/                          # Raw input data (local, usually ignored by git)
│   ├── Sabio_rk/
│   ├── BRENDA/
│   ├── Rhea/
│   ├── SKID/
│   └── RetroBioCat/
├── data_sources/                  # Source-specific and post-processing scripts
│   ├── BRENDA/
│   │   ├── step0_llm.py           # Universal LLM-based header standardizer
│   │   └── step1_join.py          # BRENDA → EnzymeML (unified format)
│   ├── sabio_rk/
│   │   ├── step0_header.py        # (Optional) SABIO header helper
│   │   └── step1_extract.py       # SABIO-RK → EnzymeML
│   ├── Rhea/
│   │   └── step1_join_modified.py # Rhea → EnzymeML
│   ├── SKID/
│   │   └── step1_skid.py          # SKiD → EnzymeML
│   ├── RetroBioCat/
│   │   ├── step0_pub.py           # Filter YAMLs by target PMIDs
│   │   └── step1_test.py          # RetroBioCat → EnzymeML
│   ├── doi_pub/
│   │   └── step0.py               # DOI → PMID converter
│   ├── mutation/                  # Mutation handling pipeline
│   │   ├── step0.py               # Fill sequences & classify WT/Mutant
│   │   ├── step1.py               # UniProt fetch + mutation retry
│   │   ├── step1_2.py             # Pattern-based mutation retry
│   │   └── step2_multi.py         # Resolve records with multiple UniProt IDs
│   └── merge/
│       └── merge_sequences.py     # Parallel sequence comparison & merge
├── schemas/
│   └── enzymeml-v2-extended.yaml  # Central extended EnzymeML schema
├── output/                        # Standardized EnzymeML YAML/JSON (local)
├── enzymeml_utils.py              # Shared helper module (Must be in PYTHONPATH)
└── requirements.txt               # Python dependencies


---

## Installation

Clone the repository and set up a Python environment:

```bash
git clone https://github.com/Xukai-YE/enzytransfer.git
cd enzytransfer

# optional but recommended
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate

# install core dependencies
pip install pandas pyyaml openpyxl requests

# optional: only if you use the LLM header standardizer
pip install openai
