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
```
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


```

---


## Quick start example (SABIO-RK)

A typical end-to-end pipeline for **SABIO-RK**:

1. **Prepare raw data**

   Place your original tables under `data/`, for example:

   ```text
   data/Sabio_rk/raw_kcat.csv
   schemas/enzymeml-v2-extended.yaml
   ```

2. **(Optional) Normalize headers with an LLM**

   ```bash
   python data_sources/BRENDA/step0_llm.py \
     --input data/Sabio_rk/raw_kcat.csv \
     --schema schemas/enzymeml-v2-extended.yaml \
     --output data/Sabio_rk/final_standardized.csv \
     --api-key sk-... \
     -v
   ```

3. **Convert SABIO-RK to unified EnzymeML**

   ```bash
   python data_sources/sabio_rk/step1_extract.py \
     --input data/Sabio_rk/final_standardized.csv \
     --schema schemas/enzymeml-v2-extended.yaml \
     --output-dir output/sabio_rk \
     --format yaml
   ```

4. **(Optional) Run mutation and merge tools**

   ```bash
   # Normalize sequences and classify wildtype vs mutant
   python data_sources/mutation/step0.py \
     --input output/sabio_rk \
     --report-dir mutation_reports \
     --ok-dir mutation_exports/ok \
     --issues-dir mutation_exports/issues

   # Retry failed mutations using UniProt sequences
   python data_sources/mutation/step1.py \
     --input mutation_exports/issues \
     --cache-dir uniprot_cache

   # Merge multiple sources, e.g. SABIO-RK + BRENDA
   python data_sources/merge/merge_sequences.py \
     output/sabio_rk \
     output/brenda \
     -o output/merged \
     -j 8
   ```

---

## Usage by component

Below are **minimal working examples** for each main script.
For full options and explanations, run:

```bash
python path/to/script.py --help
```

and/or read the docstring at the top of each script.

### 1. LLM-based header standardizer

**Script:** `data_sources/BRENDA/step0_llm.py`
**Purpose:** map arbitrary CSV/Excel column names to schema-defined field names using a large language model.

```bash
python data_sources/BRENDA/step0_llm.py \
  --input data/Sabio_rk/raw_kcat.csv \
  --schema schemas/enzymeml-v2-extended.yaml \
  --output data/Sabio_rk/final_standardized.csv \
  --api-key sk-... \
  -v
```

Key arguments:

* `--input` (required): input file (`.csv` / `.xlsx` / `.parquet` / `.json` / `.tsv`).
* `--schema` (optional): LinkML-style schema YAML (recommended).
* `--output`: output file (default: `<input>_standardized.csv`).
* `--api-key`: OpenAI API key (or set `OPENAI_API_KEY`).
* `--model`: OpenAI model (default: `gpt-4o`).
* `--temperature`: LLM temperature (default: `0.0`).
* `--max-tokens`: max LLM response tokens (default: `2500`).
* `--encoding`: force input encoding (auto-detect if not set).
* `--sep`: force CSV separator (auto-detect if not set).
* `--verbose` / `-v`: verbose output.

---

### 2. SABIO-RK → EnzymeML

**Script:** `data_sources/sabio_rk/step1_extract.py`
**Purpose:** convert SABIO-RK data into the same EnzymeML-style format as the BRENDA exporter.

```bash
python data_sources/sabio_rk/step1_extract.py \
  --input data/Sabio_rk/final_standardized.csv \
  --schema schemas/enzymeml-v2-extended.yaml \
  --output-dir output/sabio_rk \
  --format yaml \
  --group-by entry
```

Key arguments:

* `--input` (required): input CSV/TSV/Excel file.
* `--schema` (required): schema YAML file.
* `--output-dir` (required): directory for EnzymeML output.
* `--format`: `yaml`, `json`, or `both` (default: `yaml`).
* `--group-by`:

  * `entry`: one file per SABIO entry.
  * `protein_substrate`: merge rows with the same protein+substrate.
* `--limit`: limit the number of entries to process (default: `0`, meaning *no limit*).

---

### 3. BRENDA → EnzymeML

**Script:** `data_sources/BRENDA/step1_join.py`
**Purpose:** convert BRENDA CSV exports into unified EnzymeML.

```bash
python data_sources/BRENDA/step1_join.py \
  --brenda-dir data/BRENDA/ \
  --schema schemas/enzymeml-v2-extended.yaml \
  --output-dir output/brenda \
  --format yaml \
  --group-by pair
```

Key arguments:

* `--brenda-dir` (required): path to a BRENDA ZIP file or folder with BRENDA CSVs.
* `--schema` (required): schema YAML file.
* `--output-dir` (required): directory for EnzymeML output.
* `--format`: `yaml`, `json`, or `both` (default: `json`).
* `--group-by`: grouping strategy:

  * `pair`
  * `pair-ref`
  * `tn-row`
* `--limit`: maximum number of groups to process (default: `0`, meaning *no limit*).

The script includes several BRENDA-specific normalization rules (for example, kcat units, activator/inhibitor placeholders, reference filtering). See the script docstring for details.

---

### 4. Rhea → EnzymeML

**Script:** `data_sources/Rhea/step1_join_modified.py`
**Purpose:** convert Rhea TSV files into the unified EnzymeML structure, SABIO-RK/BRENDA-compatible.

```bash
python data_sources/Rhea/step1_join_modified.py \
  --input data/Rhea/ \
  --schema schemas/enzymeml-v2-extended.yaml \
  --output-dir output/rhea \
  --format yaml \
  --limit 20000
```

Key arguments:

* `--input` (required): directory containing Rhea TSV files, or a single TSV file (the parent directory will be used).
* `--schema` (required): EnzymeML v2 extended schema YAML file.
* `--output-dir` (required): directory for EnzymeML files.
* `--format`: `yaml`, `json`, or `both` (default: `yaml`).
* `--rhea-ids` / `-r`: specific Rhea IDs to process (if omitted, use all IDs present in the TSV file(s)).
* `--limit` / `-l`: maximum number of reactions to process (default: `20000`; reactions list is truncated to this length).

---

### 5. SKiD → EnzymeML

**Script:** `data_sources/SKID/step1_skid.py`
**Purpose:** convert the SKiD main dataset into EnzymeML, grouped by `(UniProt_ID, Substrate, Mutation)`.

```bash
python data_sources/SKID/step1_skid.py \
  --input data/SKID/Main_dataset_v1.xlsx \
  --schema schemas/enzymeml-v2-extended.yaml \
  --output-dir output/skid \
  --format yaml \
  --output-mode both
```

Key arguments:

* `--input` (required): SKiD Excel file (e.g. `Main_dataset_v1.xlsx`).
* `--schema` (required): schema YAML file.
* `--output-dir` (required): directory for EnzymeML output.
* `--format`: `yaml`, `json`, or `both` (default: `yaml`).
* `--output-mode`:

  * `both`: WT and mutant (default).
  * `wt_only`
  * `mutant_only`
* `--limit`: limit the number of entries to process (default: `0`, meaning *no limit*).

The script:

* Groups rows by `(UniProt_ID, Substrate, Mutation)`.
* Loads the `Unique_substrates` sheet to enrich substrate information (InChI, InChIKey, external IDs).
* Attaches external identifiers and reference metadata.

---

### 6. RetroBioCat → EnzymeML

**Script:** `data_sources/RetroBioCat/step1_test.py`
**Purpose:** convert RetroBioCat activity tables into the unified EnzymeML format.

```bash
python data_sources/RetroBioCat/step1_test.py \
  --input data/RetroBioCat/trial_activity_data.xlsx \
  --schema schemas/enzymeml-v2-extended.yaml \
  --output-dir output/retrobiocat \
  --format yaml \
  --group-by entry
```

Key arguments:

* `--input` (required): RetroBioCat Excel file.
* `--schema` (required): schema YAML file.
* `--output-dir` (required): directory for EnzymeML output.
* `--format`: `yaml`, `json`, or `both` (default: `yaml`).
* `--group-by`:

  * `entry`
  * `enzyme_substrate`
* `--limit`: maximum number of entries to process (default: `0`, meaning *no limit*).

---

### 7. RetroBioCat PMID filter

**Script:** `data_sources/RetroBioCat/step0_pub.py`
**Purpose:** select YAML files that contain any PMID from a built-in target list.

```bash
python data_sources/RetroBioCat/step0_pub.py \
  --input-dir output/retrobiocat \
  --out-dir output/retrobiocat_filtered \
  --move
```

Key arguments:

* `--input-dir` (required): directory containing `.yml` / `.yaml` files.
* `--out-dir` (required): output directory for filtered files and reports.
* `--move`: if set, move files instead of copying.

The script:

* Recursively scans YAMLs under `--input-dir`.
* Extracts PMIDs from multiple fields / URL patterns.
* Copies/moves matched YAMLs to `by_target_pmids/` under `--out-dir`.
* Writes:

  * `pmid_matches.csv`
  * `target_pmids_not_found.csv`

---

### 8. DOI → PMID converter

**Script:** `data_sources/doi_pub/step0.py`
**Purpose:** read an Excel file, map DOIs to PMIDs via NCBI, and write a new Excel with an added `PMID` column.

Usage pattern:

1. Open the script and locate the `if __name__ == "__main__":` block.

2. Set the parameters there, for example:

   ```python
   input_file = "../../data/RetroBioCat/trial_activity_data.xlsx"
   output_file = "../../data/RetroBioCat/trial_activity_data_with_pmid.xlsx"
   your_email = "your_email@example.com"
   # the column that contains DOIs in your sheet:
   doi_column = "html_doi"
   ```

3. Run:

   ```bash
   python data_sources/doi_pub/step0.py
   ```

The script prints progress and statistics and writes the updated Excel file.

---

### 9. Mutation & sequence pipeline

Mutation-related scripts live in `data_sources/mutation/` and operate on EnzymeML YAML files.

#### 9.1 Fill sequences & classify (Step 0)

**Script:** `data_sources/mutation/step0.py`
**Purpose:** write protein sequences into YAMLs, set `variant_type`, and classify records as wildtype or mutant. Also generates CSV reports.

```bash
python data_sources/mutation/step0.py \
  --input output/sabio_rk \
  --report-dir mutation_reports \
  --ok-dir mutation_exports/ok \
  --issues-dir mutation_exports/issues
```

Key arguments:

* `--input` / `-i` (required): YAML file or directory (recursively scanned).
* `--report-dir`: directory for CSV reports (default: `./mutation_reports`).
* `--ok-dir`: directory for “OK” YAML exports (default: `../../output/mutation_exports/ok`).
* `--issues-dir`: directory for “Issue” YAML exports (default: `../../output/mutation_exports/issues`).

The script normalizes `report_dir`, `ok_dir`, and `issues_dir` to absolute paths before processing.

---

#### 9.2 UniProt fetch & mutation retry (Step 1)

**Script:** `data_sources/mutation/step1.py`
**Purpose:** for records with issues, fetch UniProt sequences, retry applying mutations, and update YAMLs.

```bash
python data_sources/mutation/step1.py \
  --input mutation_exports/issues \
  --output mutation_step1_output \
  --report mutation_step1_report.csv \
  --checkpoint mutation_step1_checkpoint.json \
  --cache-dir uniprot_cache
```

Key arguments:

* `--input` / `-i` (required): directory containing issue YAML files.
* `--output` / `-o`: output directory for corrected YAMLs (default: `./mutation_exports/retry_success`).
* `--report` / `-r`: path for retry report CSV (default: `./mutation_reports/retry_report.csv`).
* `--checkpoint` / `-c`: checkpoint JSON file (default: `./mutation_reports/.checkpoint.json`).
* `--cache-dir`: directory for UniProt sequence cache (default: `./uniprot_cache`).
* `--rate-limit`: minimum seconds between UniProt API requests (default: `1.0`).
* `--resume`: resume from last checkpoint.
* `--clear-checkpoint`: clear checkpoint and start fresh.

All paths are normalized to absolute paths before use.

---

#### 9.3 Pattern-based mutation retry (Step 1.2)

**Script:** `data_sources/mutation/step1_2.py`
**Purpose:** use regex / pattern-based alignment to retry remaining failed mutations with checkpoint support.

```bash
python data_sources/mutation/step1_2.py \
  --input mutation_exports/issues \
  --output mutation_step1_2_output \
  --report mutation_step1_2_report.csv
```

Key arguments:

* `--input` / `-i` (required): directory containing issue YAML files.
* `--output` / `-o`: output directory for corrected YAMLs (default: `./mutation_exports/retry_success`).
* `--report` / `-r`: path for retry report CSV (default: `./mutation_reports/retry_report.csv`).
* `--checkpoint` / `-c`: checkpoint JSON file (default: `./mutation_reports/.checkpoint.json`).
* `--resume`: resume from last checkpoint.
* `--clear-checkpoint`: clear checkpoint and start fresh.

---

#### 9.4 Multi-UniProt-ID resolver (Step 2)

**Script:** `data_sources/mutation/step2_multi.py`
**Purpose:** handle records where `uniprotid` includes multiple IDs (separated by `;`, `,`, `|`, etc.), test mutations against each candidate, and mark resolved vs unresolved records.

```bash
python data_sources/mutation/step2_multi.py \
  --input mutation_exports/issues \
  --resolved mutation_step2_resolved \
  --unresolved mutation_step2_unresolved \
  --report mutation_step2_report.csv \
  --cache-dir uniprot_cache
```

Key arguments:

* `--input` / `-i` (required): directory containing YAML files with multiple UniProt IDs.
* `--resolved`: output directory for resolved files (default: `./mutation_exports/multi_id_resolved`).
* `--unresolved`: output directory for unresolved files (default: `./mutation_exports/multi_id_unresolved`).
* `--report` / `-r`: path for report CSV (default: `./mutation_reports/multi_id_report.csv`).
* `--checkpoint` / `-c`: checkpoint JSON file (default: `./mutation_reports/.checkpoint_multi_id.json`).
* `--cache-dir`: directory for UniProt sequence cache (default: `./uniprot_cache`).
* `--rate-limit`: minimum seconds between UniProt API requests (default: `1.0`).
* `--resume`: resume from last checkpoint.
* `--clear-checkpoint`: clear checkpoint and start fresh.

All paths are normalized to absolute paths before use.

---

### 10. Sequence comparison & merge

**Script:** `data_sources/merge/merge_sequences.py`
**Purpose:** scan one or more directories of EnzymeML YAML files, compare sequences/metadata, and merge overlapping entries. Uses multiprocessing for speed.

```bash
python data_sources/merge/merge_sequences.py \
  output/sabio_rk \
  output/brenda \
  output/rhea \
  -o output/merged \
  -j 8
```

Key arguments:

* `directories` (positional, required): one or more directories to scan (at least two, unless you pass `--allow-single-dir`).
* `-o`, `--output-dir`: merged output directory (default: `merged_output`).
* `-j`, `--jobs`: number of worker processes (default: number of CPU cores).
* `--allow-single-dir`: allow using a single directory (for testing / internal duplicate checks).
* `-v`, `--verbose`: verbose output.

---

## Extending the pipeline

* **Adding a new data source**

  * Create `data_sources/NEW_SOURCE/`.
  * Implement `step1_new_source.py` that:

    * Reads your raw tables under `data/NEW_SOURCE/`.
    * Uses the shared schema (`schemas/enzymeml-v2-extended.yaml`).
    * Constructs EnzymeML-style YAML/JSON consistent with existing `step1_*` scripts.

* **Changing the schema**

  * Update `schemas/enzymeml-v2-extended.yaml` (or your schema file).
  * Adjust `enzymeml_utils` helpers and any affected `step1_*` scripts.



