"""
04_predixcan_analysis.py
------------------------
Factorial analysis of |R| stratum x sign(R) for eGenes vs. non-eGenes, and
PrediXcan R comparison across four gene categories.

Design:
  - Primary 2x2 factorial: {high|R|, low|R|} x {R>0, R<0} for ALL genes
    Tests whether R sign predicts eGene status after |R| is held constant.
  - PrediXcan R across 4 categories (manuscript Figure 2B):
      (1) Non-eGene, low |R|
      (2) Non-eGene, high |R|
      (3) eGene, R < 0  (Enformer wrong direction)
      (4) eGene, R > 0  (Enformer correct direction)
  - Fine-tuning stability: FinetuneR by base R sign (eGenes only, Figure 2A)

Inputs:
  data/processed/gene_level_results.csv  (must contain abs_R, absR_decile
                                          columns written by 03_decile_analysis.py)

Outputs:
  data/processed/factorial_direction_analysis.csv
  data/processed/high_absR_Rpos_egenes.csv
  data/processed/high_absR_Rneg_egenes.csv
  data/processed/low_absR_Rpos_egenes.csv
  data/processed/low_absR_Rneg_egenes.csv
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
PROC_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

INPUT_FILE = os.path.join(PROC_DIR, "gene_level_results.csv")
FACTORIAL_OUT = os.path.join(PROC_DIR, "factorial_direction_analysis.csv")
HI_POS_OUT = os.path.join(PROC_DIR, "high_absR_Rpos_egenes.csv")
HI_NEG_OUT = os.path.join(PROC_DIR, "high_absR_Rneg_egenes.csv")
LO_POS_OUT = os.path.join(PROC_DIR, "low_absR_Rpos_egenes.csv")
LO_NEG_OUT = os.path.join(PROC_DIR, "low_absR_Rneg_egenes.csv")

N_DECILES = 10


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def assign_deciles(series: pd.Series, n: int = 10) -> pd.Series:
    """Assign 0-based decile index using pd.qcut."""
    labels = list(range(n))
    decile, _ = pd.qcut(series, q=n, labels=labels, retbins=True, duplicates="drop")
    return decile.astype(int)


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=1)


def print_section(title: str) -> None:
    print()
    print("=" * 60)
    print(title)
    print("=" * 60)


def summarize(label: str, arr: np.ndarray) -> None:
    if len(arr) == 0:
        print(f"  {label}: no data")
        return
    q25, q75 = np.percentile(arr, [25, 75])
    print(f"  {label}:")
    print(f"    n={len(arr)}, mean={np.mean(arr):.4f}, SD={np.std(arr, ddof=1):.4f}")
    print(f"    median={np.median(arr):.4f}, IQR=[{q25:.4f}, {q75:.4f}]")


def mw_test(a: np.ndarray, b: np.ndarray, label: str) -> float:
    if len(a) < 2 or len(b) < 2:
        print(f"  {label}: insufficient data")
        return np.nan
    stat, pval = stats.mannwhitneyu(a, b, alternative="two-sided")
    n1, n2 = len(a), len(b)
    r_rb = 1 - (2 * stat) / (n1 * n2)
    print(f"  {label}: U={stat:.1f}, p={pval:.3e}, r_rb={r_rb:.3f}")
    return pval


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print_section("04_predixcan_analysis.py  --  Factorial Direction Analysis")

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} genes, {int(df['is_egene'].sum())} eGenes")

    df = df.dropna(subset=["PearsonR"]).reset_index(drop=True)

    # Ensure abs_R and absR_decile are present (computed by 03_decile_analysis.py)
    if "abs_R" not in df.columns:
        print("  abs_R missing -- computing from PearsonR ...")
        df["abs_R"] = df["PearsonR"].abs()
    if "absR_decile" not in df.columns:
        print("  absR_decile missing -- computing from |PearsonR| ...")
        df["absR_decile"] = assign_deciles(df["abs_R"], n=N_DECILES)
    if "R_positive" not in df.columns:
        df["R_positive"] = (df["PearsonR"] > 0).astype(int)

    # -----------------------------------------------------------------------
    # Section 1: 2x2 Factorial — |R| stratum x sign(R) -> eGene rate
    # -----------------------------------------------------------------------
    print_section("1. 2x2 Factorial: |R| stratum x sign(R) -> eGene rate")
    print("(Tests whether prediction direction predicts eGene status after |R| fixed)")

    # Use top (D10) and bottom (D1) |R| deciles as the two |R| strata
    high_absR = df[df["absR_decile"] == 9]   # top |R| decile
    low_absR  = df[df["absR_decile"] == 0]   # bottom |R| decile

    factorial_rows = []
    for stratum_label, sub_stratum in [("High|R| (D10)", high_absR), ("Low|R| (D1)", low_absR)]:
        for sign_label, sub_sign in [("R>0", sub_stratum[sub_stratum["PearsonR"] > 0]),
                                      ("R<0", sub_stratum[sub_stratum["PearsonR"] < 0])]:
            n_total = len(sub_sign)
            n_egene = int(sub_sign["is_egene"].sum())
            rate = n_egene / n_total if n_total > 0 else np.nan
            print(f"  {stratum_label}, {sign_label}: n={n_total}, eGenes={n_egene}, rate={rate*100:.1f}%")
            factorial_rows.append({
                "stratum": stratum_label,
                "R_sign": sign_label,
                "n": n_total,
                "n_egene": n_egene,
                "egene_rate": rate,
            })

    # Chi-squared test within high |R| stratum (key test from manuscript)
    hi_pos = high_absR[high_absR["PearsonR"] > 0]
    hi_neg = high_absR[high_absR["PearsonR"] < 0]
    ct_hi = np.array([
        [int((hi_pos["is_egene"] == 1).sum()), int((hi_pos["is_egene"] == 0).sum())],
        [int((hi_neg["is_egene"] == 1).sum()), int((hi_neg["is_egene"] == 0).sum())],
    ])
    chi2_hi, p_hi, _, _ = stats.chi2_contingency(ct_hi)
    print(f"\n  Chi-squared (High|R| stratum, R>0 vs R<0): chi2={chi2_hi:.3f}, p={p_hi:.3f}")

    # Chi-squared test within low |R| stratum
    lo_pos = low_absR[low_absR["PearsonR"] > 0]
    lo_neg = low_absR[low_absR["PearsonR"] < 0]
    ct_lo = np.array([
        [int((lo_pos["is_egene"] == 1).sum()), int((lo_pos["is_egene"] == 0).sum())],
        [int((lo_neg["is_egene"] == 1).sum()), int((lo_neg["is_egene"] == 0).sum())],
    ])
    chi2_lo, p_lo, _, _ = stats.chi2_contingency(ct_lo)
    print(f"  Chi-squared (Low|R| stratum, R>0 vs R<0): chi2={chi2_lo:.3f}, p={p_lo:.3f}")

    # Save factorial summary
    factorial_df = pd.DataFrame(factorial_rows)
    factorial_df["chi2_within_stratum"] = np.nan
    factorial_df["chi2_pval_within_stratum"] = np.nan
    for stratum, chi2, pval in [("High|R| (D10)", chi2_hi, p_hi), ("Low|R| (D1)", chi2_lo, p_lo)]:
        mask = factorial_df["stratum"] == stratum
        factorial_df.loc[mask, "chi2_within_stratum"] = chi2
        factorial_df.loc[mask, "chi2_pval_within_stratum"] = pval
    factorial_df.to_csv(FACTORIAL_OUT, index=False)
    print(f"\n  Saved factorial summary: {FACTORIAL_OUT}")

    # Define the 4 eGene groups for downstream analysis
    hi_pos_eg = hi_pos[hi_pos["is_egene"] == 1].copy()
    hi_neg_eg = hi_neg[hi_neg["is_egene"] == 1].copy()
    lo_pos_eg = lo_pos[lo_pos["is_egene"] == 1].copy()
    lo_neg_eg = lo_neg[lo_neg["is_egene"] == 1].copy()

    print(f"\n  eGene group sizes:")
    print(f"    High|R|, R>0 eGenes: n={len(hi_pos_eg)}")
    print(f"    High|R|, R<0 eGenes: n={len(hi_neg_eg)}")
    print(f"    Low|R|,  R>0 eGenes: n={len(lo_pos_eg)}")
    print(f"    Low|R|,  R<0 eGenes: n={len(lo_neg_eg)}")

    # -----------------------------------------------------------------------
    # Section 2: PrediXcan R across 4 gene categories (Figure 2B)
    # -----------------------------------------------------------------------
    print_section("2. PrediXcan R across 4 gene categories")
    print("(Figure 2B: shows genotype-based model recovers direction where Enformer fails)")

    # Threshold for 'high |R|' vs 'low |R|' in non-eGene comparison:
    # use median |R| across all genes as threshold
    absR_median = df["abs_R"].median()
    print(f"  Median |R| threshold: {absR_median:.4f}")

    non_egenes = df[df["is_egene"] == 0]
    egenes_all = df[df["is_egene"] == 1]

    cat1 = non_egenes[non_egenes["abs_R"] <= absR_median]   # non-eGene, low |R|
    cat2 = non_egenes[non_egenes["abs_R"] > absR_median]    # non-eGene, high |R|
    cat3 = egenes_all[egenes_all["PearsonR"] < 0]            # eGene, R<0 (Enformer wrong direction)
    cat4 = egenes_all[egenes_all["PearsonR"] > 0]            # eGene, R>0 (Enformer correct direction)

    categories = [
        ("Non-eGene, low|R|",  cat1),
        ("Non-eGene, high|R|", cat2),
        ("eGene, R<0 (Enformer wrong direction)",  cat3),
        ("eGene, R>0 (Enformer correct direction)", cat4),
    ]

    pred_stats = []
    for label, sub in categories:
        if "PrediXcanR" in sub.columns:
            pred = sub["PrediXcanR"].dropna().values
        else:
            pred = np.array([])
        n_pred = len(pred)
        mean_pred = np.mean(pred) if n_pred > 0 else np.nan
        median_pred = np.median(pred) if n_pred > 0 else np.nan
        sem_pred = stats.sem(pred) if n_pred > 1 else np.nan
        print(f"\n  {label}:")
        print(f"    n(with PrediXcanR)={n_pred}, mean={mean_pred:.4f}, "
              f"median={median_pred:.4f}, SEM={sem_pred:.4f}")
        pred_stats.append({
            "category": label,
            "n_genes": len(sub),
            "n_with_PrediXcanR": n_pred,
            "PrediXcanR_mean": mean_pred,
            "PrediXcanR_median": median_pred,
            "PrediXcanR_sem": sem_pred,
        })

    # Key contrast: eGene R<0 vs non-eGene low|R|
    pred_cat3 = cat3["PrediXcanR"].dropna().values
    pred_cat1 = cat1["PrediXcanR"].dropna().values
    if len(pred_cat3) > 1 and len(pred_cat1) > 1:
        print()
        mw_test(pred_cat3, pred_cat1, "eGene R<0 vs Non-eGene low|R| PrediXcanR")

    # eGene R<0 vs eGene R>0
    pred_cat4 = cat4["PrediXcanR"].dropna().values
    if len(pred_cat3) > 1 and len(pred_cat4) > 1:
        mw_test(pred_cat3, pred_cat4, "eGene R<0 vs eGene R>0 PrediXcanR")

    # -----------------------------------------------------------------------
    # Section 3: Fine-tuning stability (Figure 2A)
    # -----------------------------------------------------------------------
    print_section("3. Fine-tuning stability: FinetuneR by base R sign (eGenes only)")

    if "PearsonRfineTuned" in df.columns:
        eg_rpos = egenes_all[egenes_all["PearsonR"] > 0]
        eg_rneg = egenes_all[egenes_all["PearsonR"] < 0]

        ft_rpos = eg_rpos["PearsonRfineTuned"].dropna().values
        ft_rneg = eg_rneg["PearsonRfineTuned"].dropna().values

        summarize("eGene R>0: FinetuneR", ft_rpos)
        summarize("eGene R<0: FinetuneR", ft_rneg)

        if len(ft_rpos) > 1 and len(ft_rneg) > 1:
            mw_pval = mw_test(ft_rpos, ft_rneg, "eGene R>0 vs R<0 FinetuneR")

        # Key result: what fraction of R<0 eGenes remain negative after fine-tuning?
        n_rneg_eg = len(eg_rneg)
        n_still_neg = int((eg_rneg["PearsonRfineTuned"].dropna() < 0).sum())
        pct_still_neg = n_still_neg / n_rneg_eg * 100 if n_rneg_eg > 0 else np.nan
        print(f"\n  R<0 eGenes where FinetuneR also < 0: {n_still_neg}/{n_rneg_eg} = {pct_still_neg:.1f}%")

        base_rpos = eg_rpos["PearsonR"].mean()
        base_rneg = eg_rneg["PearsonR"].mean()
        fine_rpos = np.mean(ft_rpos) if len(ft_rpos) > 0 else np.nan
        fine_rneg = np.mean(ft_rneg) if len(ft_rneg) > 0 else np.nan
        print(f"\n  Summary (mean values):")
        print(f"    eGene R>0: base PearsonR={base_rpos:.4f}, FinetuneR={fine_rpos:.4f}")
        print(f"    eGene R<0: base PearsonR={base_rneg:.4f}, FinetuneR={fine_rneg:.4f}")
    else:
        print("  PearsonRfineTuned column not available -- skipping")

    # -----------------------------------------------------------------------
    # Section 4: Per-|R|-decile summary (sanity check)
    # -----------------------------------------------------------------------
    print_section("4. Per-|R|-decile eGene rates (sanity check)")
    decile_rates = (
        df.groupby("absR_decile")
        .agg(n=("is_egene", "count"), n_egene=("is_egene", "sum"))
        .assign(rate=lambda x: x["n_egene"] / x["n"] * 100)
    )
    print(decile_rates.to_string())

    # -----------------------------------------------------------------------
    # Save eGene group files
    # -----------------------------------------------------------------------
    hi_pos_eg.to_csv(HI_POS_OUT, index=False)
    hi_neg_eg.to_csv(HI_NEG_OUT, index=False)
    lo_pos_eg.to_csv(LO_POS_OUT, index=False)
    lo_neg_eg.to_csv(LO_NEG_OUT, index=False)

    print(f"\nSaved High|R|,R>0 eGenes ({len(hi_pos_eg)}): {HI_POS_OUT}")
    print(f"Saved High|R|,R<0 eGenes ({len(hi_neg_eg)}): {HI_NEG_OUT}")
    print(f"Saved Low|R|,R>0 eGenes ({len(lo_pos_eg)}): {LO_POS_OUT}")
    print(f"Saved Low|R|,R<0 eGenes ({len(lo_neg_eg)}): {LO_NEG_OUT}")
    print("\nDone.")


if __name__ == "__main__":
    main()
