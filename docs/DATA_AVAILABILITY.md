# Data Availability Statement

## For manuscript submission

> All analysis code and processed data tables are freely available at GitHub
> (https://github.com/[username]/eqtl-enformer-ushape) and archived on Zenodo
> (DOI: [to be assigned upon release]). Raw input data were obtained from two
> publicly available sources: (1) the EnformerAssessment GitHub repository
> (github.com/mostafavilabuw/EnformerAssessment; Sasse et al., 2023, Nature Genetics),
> which provides Enformer per-gene Pearson R correlations and supplementary tables
> derived from the ROSMAP cohort without individual-level data; and (2) the GTEx Portal
> REST API v2 (gtexportal.org/api/v2), which provides GTEx v8 Brain Cortex eGene
> summary statistics. No controlled-access or individual-level data are used in this study.

---

## Data sources

### 1. Enformer per-gene correlations

- **Source:** Sasse A, et al. "Benchmarking of deep neural networks for predicting
  personal gene expression from DNA sequence highlights the necessity for multiomic input."
  *Nat Genet.* 2023;55(9):1517–1527.
- **Repository:** https://github.com/mostafavilabuw/EnformerAssessment
- **Files used:**
  - `Data/Prediction_correlationsCageAdultBrain_Allstats.txt` — per-gene Pearson R
  - `Data/SupplementaryTable1.tsv` — gene coordinates, fine-tuned R, PrediXcan R
- **License:** MIT (per repository)
- **Access:** Public, no registration required

### 2. GTEx v8 Brain Cortex eGenes

- **Source:** GTEx Consortium. "The GTEx Consortium atlas of genetic regulatory effects
  across human tissues." *Science.* 2020;369(6509):1318–1330.
- **API:** https://gtexportal.org/api/v2/association/egene
- **Parameters:** tissueSiteDetailId=Brain_Cortex, datasetId=gtex_v8
- **n:** 9,082 eGenes from 205 Brain Cortex samples
- **License:** GTEx data use policy (summary statistics are publicly available without
  individual-level data access requirements)
- **Access:** Public API, no registration required

---

## Processed data files

The following processed files are available in `data/processed/`:

| File | n rows | Description |
|------|--------|-------------|
| `gene_level_results.csv` | 6,808 | All genes with Enformer R, is_egene flag, covariates, decile |
| `classA_D10_eGenes.csv` | 502 | Class A eGenes (D10, best Enformer performance) |
| `classB_D1_eGenes.csv` | 453 | Class B eGenes (D1, worst Enformer performance) |
| `decile_summary.csv` | 10 | Per-decile eGene rate, bootstrap CI, mean R |
| `../tables/table1_classA_classB.csv` | 2 | Class A vs B summary statistics |

These files are derived from publicly available summary statistics and contain no
individual-level genomic or expression data. No controlled-access approval is needed
to use them.

---

## What is NOT available here

- Individual-level ROSMAP WGS genotype data (controlled access via Rush Alzheimer's
  Disease Center: rush.edu/research/our-researchers/rosmap)
- Individual-level ROSMAP RNA-seq expression data (same controlled-access route)
- GTEx individual-level genotype or expression data (controlled access via dbGaP:
  phs000424.v9.p2)

These data were not used in this study. All analyses used only the publicly available
summary outputs from Sasse et al. (2023) and the GTEx Portal API.
