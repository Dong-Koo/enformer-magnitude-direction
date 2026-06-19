# Raw Data Sources

This directory contains files downloaded by `code/01_download_data.py`.

## Files

| Filename | Source | Downloaded from |
|----------|--------|----------------|
| `enformer_correlations.txt` | Sasse et al. 2023 | github.com/mostafavilabuw/EnformerAssessment |
| `supp_table1.tsv` | Sasse et al. 2023 | github.com/mostafavilabuw/EnformerAssessment |
| `gtex_brain_cortex_egenes.txt` | GTEx v8 | gtexportal.org/api/v2 |

## Why raw data is not committed to git

Raw data files (particularly `enformer_correlations.txt` at ~600 KB and
`supp_table1.tsv` at ~2 MB) are excluded from git to avoid bloating the repository.
Run `python code/01_download_data.py` to download them.

Both sources are public and stable. The EnformerAssessment GitHub repository
is associated with a published Nature Genetics paper and has been stable since 2023.
The GTEx Portal API v2 has been available since 2021.

## Data format notes

**enformer_correlations.txt**: Space-delimited. Lines starting with `#` are comments
(including the header line `# gene PearsonR Pvalue MeanObs StdObs MeanEnf StdEnf`).
Read with `header=None`, `comment='#'`.

**supp_table1.tsv**: Tab-delimited. First line is a `#`-prefixed header with 25 columns.
Read with `header=None`, `comment='#'`, and manually specify column names.

**gtex_brain_cortex_egenes.txt**: One ENSEMBL ID per line, version suffix stripped.
E.g.: `ENSG00000000419` (not `ENSG00000000419.12`).
