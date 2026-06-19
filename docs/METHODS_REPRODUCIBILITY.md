# Methods and Reproducibility Guide

This document describes the exact procedures used in the analysis and how to verify
each step independently.

---

## Software environment

| Package | Version used | Purpose |
|---------|-------------|---------|
| Python | 3.11+ | Runtime |
| pandas | 2.x | Data loading, merging |
| numpy | 1.26+ | Numerical operations |
| scipy | 1.12+ | Statistical tests (Fisher, Mann-Whitney, Spearman) |
| statsmodels | 0.14+ | Logistic regression (Logit MLE) |
| scikit-learn | 1.4+ | StandardScaler for covariate z-scoring |
| matplotlib | 3.8+ | Figures |

Install: `pip install -r code/requirements.txt`

---

## Step-by-step method

### Step 1. Data acquisition (code/01_download_data.py)

**Enformer correlations:**
- URL: `https://raw.githubusercontent.com/mostafavilabuw/EnformerAssessment/main/Data/Prediction_correlationsCageAdultBrain_Allstats.txt`
- Format: space-delimited, lines beginning with `#` are comments/headers
- Columns: `gene PearsonR Pvalue MeanObs StdObs MeanEnf StdEnf`

**Supplementary Table 1:**
- URL: `https://raw.githubusercontent.com/mostafavilabuw/EnformerAssessment/main/Data/SupplementaryTable1.tsv`
- Format: tab-delimited, first line is `#`-commented header
- 25 columns including genomic coordinates, fine-tuned R, PrediXcan R

**GTEx eGenes:**
- Endpoint: `https://gtexportal.org/api/v2/association/egene`
- Parameters: `tissueSiteDetailId=Brain_Cortex`, `datasetId=gtex_v8`, `itemsPerPage=2000`
- Paginated: fetch pages 0–4 until `len(data) < itemsPerPage`
- Field extracted: `gencodeId` (ENSEMBL ID with version, e.g. `ENSG00000123.5`)
- Version suffix stripped: `.5` removed to yield `ENSG00000123`

### Step 2. Gene ID matching (code/02_match_gene_ids.py)

- `corr['gene']` contains versioned ENSEMBL IDs → strip version: `.str.split('.').str[0]`
- `supp['gene_id']` same procedure
- Inner merge on stripped ID → 6,808 genes retained after removing NaN rows
- `is_egene` flag: 1 if gene_id in GTEx eGene set, 0 otherwise
- Gene length: `|end_hg38 - start_hg38|`
- `log10_gene_length`: `log10(gene_length + 1)` (pseudocount avoids log(0))

### Step 3. Decile analysis (code/03_decile_analysis.py)

- `pd.qcut(df['PearsonR'], 10, labels=False)` → decile 0–9 (equal gene count per decile)
- Per-decile eGene rate: `sum(is_egene) / n`
- Bootstrap CI: 2,000 resamples per decile, `seed=42`, percentile method [2.5, 97.5]

**Logistic regression:**
```python
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
df[['MeanObs_z','StdObs_z','len_z']] = scaler.fit_transform(
    df[['MeanObs','StdObs','log10_gene_length']])
df['R_sq'] = df['PearsonR']**2

X = sm.add_constant(df[['PearsonR','R_sq','MeanObs_z','StdObs_z','len_z']].values)
result = sm.Logit(df['is_egene'].values, X).fit()
```

Key output: `result.pvalues[2]` → p-value for R² term (index 2 = third predictor).

### Step 4. PrediXcan analysis (code/04_predixcan_analysis.py)

- Class A = genes with `decile == 9` AND `is_egene == 1` (n = 502)
- Class B = genes with `decile == 0` AND `is_egene == 1` (n = 453)
- PrediXcan comparison: `scipy.stats.mannwhitneyu(A_px, B_px, alternative='two-sided')`
- Covariate-adjusted: `statsmodels.api.OLS(PrediXcanR ~ MeanObs_z + len_z + is_classB)`
- Fine-tuning delta: `PearsonRfineTuned - PearsonR` (per gene)

### Step 5. Figures (code/05_make_figures.py)

- All figures use `matplotlib` with `backend='Agg'` (no display required)
- DPI: 300 for PNG, vector for PDF
- Bootstrap CIs from Step 3 are passed as yerr arrays
- Random seeds fixed for scatter sample points

---

## Key statistical decisions

| Decision | Rationale |
|----------|-----------|
| Mann-Whitney U (not t-test) | PrediXcan R distributions are non-normal; MWU is robust |
| Logistic regression with quadratic term | Tests U-shape formally vs. monotone alternative |
| Bootstrap CI (not binomial CI) | Accounts for possible non-independence from decile assignment |
| Standardized covariates | Allows coefficient comparison across predictors on different scales |
| `pd.qcut` for deciles | Ensures equal gene counts per decile; prevents edge-bin artifacts |

---

## Verification

To verify the key result (logistic R² term p-value) independently:

```python
import pandas as pd, numpy as np, statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

df = pd.read_csv('data/processed/gene_level_results.csv')
df['R_sq'] = df['PearsonR']**2
scaler = StandardScaler()
df[['m','s','l']] = scaler.fit_transform(df[['MeanObs','StdObs','log10_gene_length']])
X = sm.add_constant(df[['PearsonR','R_sq','m','s','l']].values)
result = sm.Logit(df['is_egene'].values, X).fit(disp=0)
print(result.pvalues)  # index 2 should be ~5e-75
```

Expected output (R² term): p ≈ 5.0×10⁻⁷⁵

---

## Notes for reviewers

- All randomness is seeded (`np.random.seed(42)`)
- No data transformation beyond log10(gene_length) and z-scoring
- No model selection or hyperparameter tuning
- All p-values are two-sided unless explicitly stated
- No multiple comparison correction applied to the primary logistic regression
  (single pre-specified test for the quadratic term)
