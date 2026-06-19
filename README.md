# Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but Fails to Predict Effect Direction

**A Magnitude–Direction Dissociation in Individual-Level Expression Prediction**

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2FXXXXXXXX-b31b1b)](https://doi.org/10.1101/XXXXXXXX)
[![Zenodo](https://img.shields.io/badge/Zenodo-DOI%3A10.5281%2Fzenodo.20754856-blue)](https://doi.org/10.5281/zenodo.20754856)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](docs/LICENSE)

**Author:** DongKoo Lee — Independent Researcher, Gwangju-si, Republic of Korea — ceo@nrootm.com — ORCID: [to be assigned]

---

## Abstract

We cross-reference Enformer's individual-level per-gene Pearson R (Sasse et al. 2023, *Nature Genetics*, ROSMAP cohort n=839) with GTEx v8 Brain Cortex eGenes (n=9,082) across 6,808 genes, and report a fundamental dissociation: prediction *magnitude* (|R|) monotonically predicts cis-eQTL gene identity (eGene rate 31.6% → 82.5%; |R| model AIC=8,274 vs. R+R² AIC=8,378, ΔAIC=104), but prediction *direction* for eGenes is near-random (54.0% correct; max D10=62.5%). Direction errors are irreversible by fine-tuning on ROSMAP (Mann–Whitney p≈10⁻¹⁰⁰), yet PrediXcan recovers correct direction for R<0 eGenes (PrediXcanR=0.295), implicating Enformer's non-brain ENCODE training corpus as the mechanistic source of directional failure.

---

## Key results

| Result | Value |
|--------|-------|
| eGene rate, lowest \|R\| decile (D1) | 31.6% |
| eGene rate, highest \|R\| decile (D10) | 82.5% |
| \|R\| model AIC advantage over R+R² | ΔAIC = 104 |
| Direction accuracy for eGenes (all deciles) | 54.0% |
| PrediXcanR for R<0 eGenes | 0.295 |

---

## Quick start

```bash
git clone https://github.com/[username]/eqtl-enformer-magnitude-direction
cd eqtl-enformer-magnitude-direction
pip install -r code/requirements.txt
python code/01_download_data.py
python code/02_match_gene_ids.py
python code/03_decile_analysis.py
python code/04_predixcan_analysis.py
python code/05_make_figures.py
```

Expected runtime: ~5 minutes on a laptop (download-limited). No GPU required.

Alternatively, use the conda environment:

```bash
conda env create -f code/environment.yml
conda activate eqtl-magnitude-direction
python code/01_download_data.py
# ... continue as above
```

---

## Repository structure

```
eqtl_enformer_ushape/
├── README.md                        <- This file
├── CITATION.cff                     <- Machine-readable citation
├── .zenodo.json                     <- Zenodo metadata for DOI assignment
│
├── manuscript/
│   ├── manuscript.html              <- Full manuscript (print to PDF)
│   └── figures/
│       ├── fig1_magnitude_direction.{png,pdf}
│       ├── fig2_finetune_direction.{png,pdf}
│       └── fig3_model_comparison.{png,pdf}
│
├── code/
│   ├── 01_download_data.py          <- Fetch Enformer correlations + GTEx eGenes
│   ├── 02_match_gene_ids.py         <- Merge by ENSEMBL ID, assign deciles
│   ├── 03_decile_analysis.py        <- Logistic regression, 2x2 factorial, CIs
│   ├── 04_predixcan_analysis.py     <- PrediXcan vs Enformer direction comparison
│   ├── 05_make_figures.py           <- Reproduce all manuscript figures
│   ├── requirements.txt
│   └── environment.yml
│
├── data/
│   └── processed/
│       ├── gene_level_results.csv   <- 6,808 genes, all metrics
│       ├── classA_D10_eGenes.csv    <- D10 eGenes (legacy compatibility)
│       └── classB_D1_eGenes.csv     <- D1 eGenes (legacy compatibility)
│
└── docs/
    ├── README.md                    <- Docs index and publication workflow
    ├── cover_letter.md              <- Journal submission cover letter
    ├── DATA_AVAILABILITY.md         <- Full data provenance statement
    ├── METHODS_REPRODUCIBILITY.md   <- Step-by-step methods guide
    ├── SUBMISSION_CHECKLIST.md      <- Pre-submission checklist
    ├── researcher_email.md          <- Outreach email templates
    └── LICENSE                      <- MIT license text
```

---

## Data availability

No controlled-access or individual-level data are used. All inputs are publicly available:

| Dataset | Source | Access |
|---------|--------|--------|
| Enformer per-gene Pearson R | [mostafavilabuw/EnformerAssessment](https://github.com/mostafavilabuw/EnformerAssessment) (Sasse et al. 2023) | Public |
| GTEx v8 Brain Cortex eGenes | [GTEx Portal API v2](https://gtexportal.org/api/v2) | Public |

See `docs/DATA_AVAILABILITY.md` for the full provenance statement.

---

## Citation

```bibtex
@article{lee2026enformer,
  author    = {Lee, DongKoo},
  title     = {Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but
               Fails to Predict Effect Direction: A Magnitude--Direction Dissociation
               in Individual-Level Expression Prediction},
  journal   = {bioRxiv},
  year      = {2026},
  doi       = {10.1101/XXXXXXXX},
  url       = {https://doi.org/10.1101/XXXXXXXX},
  note      = {Preprint}
}
```

---

## License

Code: MIT License (see `docs/LICENSE`).
Text, figures, and tables: CC BY 4.0.

---

> **Not peer reviewed.** This is a preprint. Please cite the bioRxiv version until peer-reviewed publication.
