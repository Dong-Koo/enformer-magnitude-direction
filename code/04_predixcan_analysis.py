"""
04_predixcan_analysis.py
------------------------
Compare PrediXcan R and expression statistics between Class A (D10 eGenes)
and Class B (D1 eGenes), and fit a covariate-adjusted OLS model.

Inputs:
  data/processed/gene_level_results.csv   (must already contain 'decile' column
                                           written by 03_decile_analysis.py)

Outputs:
  data/processed/classA_D10_eGenes.csv
  data/processed/classB_D1_eGenes.csv
  (prints all statistical results to stdout)
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
PROC_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

INPUT_FILE = os.path.join(PROC_DIR, "gene_level_results.csv")
CLASSA_OUT = os.path.join(PROC_DIR, "classA_D10_eGenes.csv")
CLASSB_OUT = os.path.join(PROC_DIR, "classB_D1_eGenes.csv")

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
    """Print mean ± SD, median [IQR] for an array."""
    if len(arr) == 0:
        print(f"  {label}: no data")
        return
    q25, q75 = np.percentile(arr, [25, 75])
    print(f"  {label}:")
    print(f"    n={len(arr)}, mean={np.mean(arr):.4f}, SD={np.std(arr, ddof=1):.4f}")
    print(f"    median={np.median(arr):.4f}, IQR=[{q25:.4f}, {q75:.4f}]")


def mw_test(a: np.ndarray, b: np.ndarray, label: str) -> None:
    """Mann-Whitney U test between arrays a and b."""
    if len(a) < 2 or len(b) < 2:
        print(f"  {label}: insufficient data for test")
        return
    stat, pval = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Effect size: rank-biserial correlation
    n1, n2 = len(a), len(b)
    r_rb = 1 - (2 * stat) / (n1 * n2)
    print(f"  {label}: U={stat:.1f}, p={pval:.3e}, r_rb={r_rb:.3f}")


def cohend(a: np.ndarray, b: np.ndarray) -> float:
    """Cohen's d effect size."""
    pooled_sd = np.sqrt((np.std(a, ddof=1) ** 2 + np.std(b, ddof=1) ** 2) / 2)
    if pooled_sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled_sd


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print_section("04_predixcan_analysis.py  —  Class A vs Class B eGene Analysis")

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} genes, {df['is_egene'].sum()} eGenes")

    # Ensure decile column is present
    if "decile" not in df.columns:
        print("  'decile' column missing — computing now ...")
        df = df.dropna(subset=["PearsonR"]).reset_index(drop=True)
        df["decile"] = assign_deciles(df["PearsonR"], n=N_DECILES)
    else:
        df = df.dropna(subset=["PearsonR"]).reset_index(drop=True)

    # -----------------------------------------------------------------------
    # Define classes
    # -----------------------------------------------------------------------
    class_a = df[(df["decile"] == 9) & (df["is_egene"] == 1)].copy()  # D10
    class_b = df[(df["decile"] == 0) & (df["is_egene"] == 1)].copy()  # D1

    print(f"\n  Class A (D10 eGenes): n={len(class_a)}")
    print(f"  Class B (D1  eGenes): n={len(class_b)}")

    # -----------------------------------------------------------------------
    # 1. PrediXcan R comparison
    # -----------------------------------------------------------------------
    print_section("1. PrediXcan R Comparison")

    a_pred = class_a["PrediXcanR"].dropna().values
    b_pred = class_b["PrediXcanR"].dropna().values

    summarize("Class A PrediXcanR", a_pred)
    summarize("Class B PrediXcanR", b_pred)
    if len(a_pred) > 0 and len(b_pred) > 0:
        mw_test(a_pred, b_pred, "Class A vs B PrediXcanR")
        d = cohend(a_pred, b_pred)
        print(f"  Cohen's d (A - B): {d:.3f}")

    # -----------------------------------------------------------------------
    # 2. Expression statistics comparison
    # -----------------------------------------------------------------------
    print_section("2. Expression Statistics (MeanObs, StdObs)")

    for col in ["MeanObs", "StdObs"]:
        if col not in df.columns:
            continue
        a_vals = class_a[col].dropna().values
        b_vals = class_b[col].dropna().values
        summarize(f"Class A {col}", a_vals)
        summarize(f"Class B {col}", b_vals)
        mw_test(a_vals, b_vals, f"Class A vs B {col}")
        print()

    # -----------------------------------------------------------------------
    # 3. Fine-tuning delta: PearsonRfineTuned - PearsonR
    # -----------------------------------------------------------------------
    print_section("3. Fine-Tuning Delta (PearsonRfineTuned - PearsonR)")

    if "PearsonRfineTuned" in df.columns:
        class_a["ft_delta"] = class_a["PearsonRfineTuned"] - class_a["PearsonR"]
        class_b["ft_delta"] = class_b["PearsonRfineTuned"] - class_b["PearsonR"]

        a_delta = class_a["ft_delta"].dropna().values
        b_delta = class_b["ft_delta"].dropna().values

        summarize("Class A fine-tuning delta", a_delta)
        summarize("Class B fine-tuning delta", b_delta)
        if len(a_delta) > 0 and len(b_delta) > 0:
            mw_test(a_delta, b_delta, "Class A vs B fine-tuning delta")
            d = cohend(a_delta, b_delta)
            print(f"  Cohen's d (A - B): {d:.3f}")
    else:
        print("  PearsonRfineTuned not available — skipping")

    # -----------------------------------------------------------------------
    # 4. Covariate-adjusted OLS: PrediXcanR ~ MeanObs_z + log10_gl_z + is_classB
    # -----------------------------------------------------------------------
    print_section("4. Covariate-Adjusted OLS: PrediXcanR ~ covariates + is_classB")

    ols_df = pd.concat([class_a, class_b], ignore_index=True).copy()
    ols_df["is_classB"] = (ols_df["decile"] == 0).astype(int)

    needed = ["PrediXcanR", "MeanObs", "log10_gene_length"]
    ols_df = ols_df.dropna(subset=needed)
    print(f"  OLS sample size (A+B eGenes with complete data): {len(ols_df)}")

    if len(ols_df) >= 10:
        ols_df["MeanObs_z"] = zscore(ols_df["MeanObs"])
        ols_df["log10_gl_z"] = zscore(ols_df["log10_gene_length"])

        X = sm.add_constant(ols_df[["MeanObs_z", "log10_gl_z", "is_classB"]])
        y = ols_df["PrediXcanR"]

        ols_model = sm.OLS(y, X).fit()
        print(ols_model.summary())

        coef = ols_model.params.get("is_classB", np.nan)
        pval = ols_model.pvalues.get("is_classB", np.nan)
        ci = ols_model.conf_int().loc["is_classB"] if "is_classB" in ols_model.params else [np.nan, np.nan]
        print(f"\n  is_classB coef: {coef:.4f}, p={pval:.3e}, 95%CI [{ci[0]:.4f}, {ci[1]:.4f}]")
    else:
        print("  Too few observations for OLS — skipping")

    # -----------------------------------------------------------------------
    # 5. Decile eGene rates summary (for sanity check)
    # -----------------------------------------------------------------------
    print_section("5. Per-Decile eGene Rates (sanity check)")
    decile_rates = (
        df.groupby("decile")
        .agg(n=("is_egene", "count"), n_egene=("is_egene", "sum"))
        .assign(rate=lambda x: x["n_egene"] / x["n"] * 100)
    )
    print(decile_rates.to_string())

    # -----------------------------------------------------------------------
    # Save class files
    # -----------------------------------------------------------------------
    class_a.to_csv(CLASSA_OUT, index=False)
    class_b.to_csv(CLASSB_OUT, index=False)
    print(f"\nSaved Class A eGenes ({len(class_a)} genes): {CLASSA_OUT}")
    print(f"Saved Class B eGenes ({len(class_b)} genes): {CLASSB_OUT}")
    print("\nDone.")


if __name__ == "__main__":
    main()
