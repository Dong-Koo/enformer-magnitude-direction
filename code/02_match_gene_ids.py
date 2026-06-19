"""
02_match_gene_ids.py
--------------------
Match gene IDs across data sources and build a unified gene-level table.

Inputs (data/raw/):
  enformer_correlations.txt  — space-delimited, #-commented
  supp_table1.tsv            — tab-delimited, 25 columns
  gtex_brain_cortex_egenes.txt — one ENSEMBL ID per line

Output:
  data/processed/gene_level_results.csv

Columns in output:
  gene_id, gene_name, chr_hg38,
  PearsonR, Pvalue, MeanObs, StdObs, MeanEnf, StdEnf,
  PearsonRfineTuned, PrediXcanR, SuSieCausalBrainhg38,
  is_egene, gene_length, log10_gene_length
"""

import os

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
PROC_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
os.makedirs(PROC_DIR, exist_ok=True)

CORR_FILE = os.path.join(RAW_DIR, "enformer_correlations.txt")
SUPP_FILE = os.path.join(RAW_DIR, "supp_table1.tsv")
EGENE_FILE = os.path.join(RAW_DIR, "gtex_brain_cortex_egenes.txt")
OUTPUT_FILE = os.path.join(PROC_DIR, "gene_level_results.csv")


def strip_version(gene_id: str) -> str:
    """Remove Ensembl version suffix: ENSG00000XXXXX.Y -> ENSG00000XXXXX."""
    if isinstance(gene_id, str):
        return gene_id.split(".")[0]
    return gene_id


# ---------------------------------------------------------------------------
# 1. Load Enformer correlations
#    Columns: gene PearsonR Pvalue MeanObs StdObs MeanEnf StdEnf
# ---------------------------------------------------------------------------
def load_enformer_correlations(path: str) -> pd.DataFrame:
    print(f"Loading enformer correlations: {path}")
    corr = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["gene_id_raw", "PearsonR", "Pvalue", "MeanObs", "StdObs", "MeanEnf", "StdEnf"],
    )
    corr["gene_id"] = corr["gene_id_raw"].apply(strip_version)
    corr = corr.drop(columns=["gene_id_raw"])
    print(f"  Loaded {len(corr)} rows, {corr['gene_id'].nunique()} unique gene IDs")
    return corr


# ---------------------------------------------------------------------------
# 2. Load Supplementary Table 1 (25 columns)
#    Expected header order per specification:
#    gene_name gene_id chr_hg38 start_hg38 end_hg38 strand_hg38 tss_hg38
#    MeanObs StdObs PearsonRfineTuned PvalueRfineTuned MeanFineTuned StdFineTuned
#    PearsonRCAGE PvalueRCAGE MeanCAGE StdCAGE
#    RmeanRandom Rstdrandom TstatCAGE Log10Pvalue Log10BHPvalue
#    PrediXcanR SuSieCausalBrainhg38 SusieCausalCortexhg38
# ---------------------------------------------------------------------------
SUPP_COLS = [
    "gene_name", "gene_id_raw", "chr_hg38", "start_hg38", "end_hg38",
    "strand_hg38", "tss_hg38",
    "MeanObs_supp", "StdObs_supp",
    "PearsonRfineTuned", "PvalueRfineTuned", "MeanFineTuned", "StdFineTuned",
    "PearsonRCAGE", "PvalueRCAGE", "MeanCAGE", "StdCAGE",
    "RmeanRandom", "Rstdrandom", "TstatCAGE", "Log10Pvalue", "Log10BHPvalue",
    "PrediXcanR", "SuSieCausalBrainhg38", "SusieCausalCortexhg38",
]


def load_supp_table1(path: str) -> pd.DataFrame:
    print(f"Loading supplementary table 1: {path}")
    # Try with header first; if column count mismatches, fall back to named cols
    raw = pd.read_csv(path, sep="\t", comment="#", header=None, low_memory=False)
    n_cols = raw.shape[1]
    print(f"  Detected {n_cols} columns")

    if n_cols == 25:
        raw.columns = SUPP_COLS
        # Drop any row that looks like a text header (first row may be header)
        if raw.iloc[0]["gene_id_raw"] in ("gene_id", "GeneID", "gene"):
            raw = raw.iloc[1:].reset_index(drop=True)
    else:
        # Try reading with header row
        raw = pd.read_csv(path, sep="\t", comment="#", header=0, low_memory=False)
        # Map known column names
        rename_map = {}
        for col in raw.columns:
            if col.lower() in ("gene_id", "geneid"):
                rename_map[col] = "gene_id_raw"
            elif col.lower() in ("gene_name", "genename", "gene"):
                rename_map[col] = "gene_name"
        raw = raw.rename(columns=rename_map)
        if "gene_id_raw" not in raw.columns and "gene_id" in raw.columns:
            raw = raw.rename(columns={"gene_id": "gene_id_raw"})
        print(f"  Columns after rename: {list(raw.columns)}")

    supp = raw.copy()
    supp["gene_id"] = supp["gene_id_raw"].apply(strip_version)

    # Coerce numeric columns
    for col in ["start_hg38", "end_hg38", "PearsonRfineTuned", "PrediXcanR",
                "SuSieCausalBrainhg38"]:
        if col in supp.columns:
            supp[col] = pd.to_numeric(supp[col], errors="coerce")

    print(f"  Loaded {len(supp)} rows, {supp['gene_id'].nunique()} unique gene IDs")
    return supp


# ---------------------------------------------------------------------------
# 3. Load GTEx eGene list
# ---------------------------------------------------------------------------
def load_egenes(path: str) -> set:
    print(f"Loading eGene list: {path}")
    with open(path) as fh:
        ids = {strip_version(line.strip()) for line in fh if line.strip()}
    print(f"  Loaded {len(ids)} eGene IDs")
    return ids


# ---------------------------------------------------------------------------
# 4. Merge and save
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 60)
    print("02_match_gene_ids.py  —  Merging gene ID tables")
    print("=" * 60)

    corr = load_enformer_correlations(CORR_FILE)
    supp = load_supp_table1(SUPP_FILE)
    egenes = load_egenes(EGENE_FILE)

    # Keep only needed columns from supp
    supp_cols_keep = [
        "gene_id", "gene_name", "chr_hg38", "start_hg38", "end_hg38",
        "PearsonRfineTuned", "PrediXcanR", "SuSieCausalBrainhg38",
    ]
    supp_cols_keep = [c for c in supp_cols_keep if c in supp.columns]
    supp_slim = supp[supp_cols_keep].drop_duplicates(subset="gene_id")

    # Inner merge on gene_id
    merged = pd.merge(corr, supp_slim, on="gene_id", how="inner")
    print(f"\nAfter merge (corr ∩ supp): {len(merged)} genes")

    # eGene flag
    merged["is_egene"] = merged["gene_id"].isin(egenes).astype(int)

    # Gene length from hg38 coordinates
    if "start_hg38" in merged.columns and "end_hg38" in merged.columns:
        merged["gene_length"] = (merged["end_hg38"] - merged["start_hg38"]).abs()
    else:
        merged["gene_length"] = np.nan

    merged["log10_gene_length"] = np.log10(merged["gene_length"].clip(lower=1))

    # Final column order
    out_cols = [
        "gene_id", "gene_name", "chr_hg38",
        "PearsonR", "Pvalue", "MeanObs", "StdObs", "MeanEnf", "StdEnf",
        "PearsonRfineTuned", "PrediXcanR", "SuSieCausalBrainhg38",
        "is_egene", "gene_length", "log10_gene_length",
    ]
    out_cols = [c for c in out_cols if c in merged.columns]
    result = merged[out_cols].reset_index(drop=True)

    result.to_csv(OUTPUT_FILE, index=False)
    print(f"\nSaved: {OUTPUT_FILE}")

    # Summary stats
    n_total = len(result)
    n_egene = result["is_egene"].sum()
    n_enformer = len(corr["gene_id"].unique())
    n_supp_uniq = supp["gene_id"].nunique()
    overlap = result["gene_id"].isin(egenes).sum()
    overlap_rate = overlap / n_egene * 100 if n_egene > 0 else 0.0

    print()
    print("Summary")
    print("-------")
    print(f"  Enformer genes       : {n_enformer}")
    print(f"  Supp Table 1 genes   : {n_supp_uniq}")
    print(f"  Merged gene set      : {n_total}")
    print(f"  eGenes in merged set : {n_egene}")
    print(f"  eGene overlap rate   : {overlap_rate:.1f}%  ({overlap}/{n_egene})")
    print(f"  Missing PrediXcanR   : {result['PrediXcanR'].isna().sum()}")
    print("Done.")


if __name__ == "__main__":
    main()
