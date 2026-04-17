# Analysis Code — Shivakumar et al.
## *Mechanistic insights into microsomes in a Pichia pastoris cell-free system*

---

### Overview

This repository contains the Python analysis and figure generation code supporting
the proteomic profiling and related sections of the manuscript.
Three scripts are provided, corresponding to distinct stages of the analysis pipeline:

| Script | Purpose |
|--------|---------|
| `normalisation_and_export.py` | Two-step normalisation pipeline; exports all underlying panel data to Excel |
| `baseline_proteome.py` | Generates QC and main figures for control condition proteomics (Figure 2) |
| `BfA_proteome.py` | Generates heatmap and scatter plots for BfA vs control comparison (Figure 5) |

---

### Input files required

All three scripts expect the following files in the same working directory:

- `Proteomics_raw_data.xlsx` — raw LC-MS intensity table (687 protein entries × 6 fractions)
- `pichia_to_cerevisiae_orthologs.csv` — curated K. phaffii → S. cerevisiae ortholog mappings
- `pichia_uniprot_entries.json` — UniProt GO annotations for detected proteins (BfA_proteome.py only)

---

### Installation

Python 3.9 or later is required. All dependencies can be installed with:

```bash
pip install -r requirements.txt
```

---

### Running the scripts

Each script is self-contained and can be run directly from the command line
from the directory containing the input files:

```bash
python normalisation_and_export.py
python baseline_proteome.py
python BfA_proteome.py
```

Output files are written to the working directory. PDF and PNG versions of
all figures are saved automatically.

---

### Functional annotation

Proteins were categorised using a three-tiered approach:

1. Direct lookup against a curated GO Slim map of K. phaffii gene symbols
2. S. cerevisiae ortholog mapping for uncharacterised K. phaffii entries
3. UniProt GO annotation for any remaining unannotated proteins

The full annotation map is embedded in `normalisation_and_export.py` and
`baseline_proteome.py` and can be extended as needed.

---

### Protein identification

Peptide identification and protein quantification were carried out using
ProteinLynx Global SERVER v3.01 (Waters), with data-independent acquisition
(LC-MS^E^) and searches against the UniProt K. phaffii (CBS 7435) reference
proteome. Gene-level intensities were obtained by summing across all accessions
sharing the same gene symbol (paralogs and strain variants).

---

### Dependencies

See `requirements.txt`. Key packages:

- pandas, numpy — data handling and normalisation
- matplotlib, seaborn — figure generation
- scipy — hierarchical clustering, correlation
- openpyxl — Excel export
- matplotlib-venn — Venn diagram rendering (optional; falls back to built-in renderer)

---

### Contact

Tejasvi Shivakumar — Department of Chemical Engineering, Imperial College London
Karen M. Polizzi (corresponding) — k.polizzi@imperial.ac.uk
