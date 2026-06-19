"""
03_decile_analysis.py
---------------------
Decile analysis of Enformer per-gene Pearson R vs. GTEx Brain Cortex eGene status.

Key analyses:
  - Assign deciles (D1–D10) by PearsonR via pd.qcut
  - Bootstrap 95% CIs on per-decile eGene rates (2000 resamples)
  - Logistic regression: is_egene ~ intercept + PearsonR + PearsonR² +
      MeanObs_z + StdObs_z + log10_gene_length_z
  - Save summary tables

Outputs:
  data/processed/decile_summary.csv
  manuscript/tables/table1_classA_classB.csv
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
TABLE_DIR = os.path.join(PROJECT_ROOT, "manuscript", "tables")
os.makedirs(TABLE_DIR, exist_ok=True)

INPUT_FILE = os.path.join(PROC_DIR, "gene_level_results.csv")
DECILE_OUT = os.path.join(PROC_DIR, "decile_summary.csv")
CLASS_TABLE_OUT = os.path.join(TABLE_DIR, "table1_classA_classB.csv")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N_BOOT = 2000
BOOT_SEED = 42
N_DECILES = 10


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def assign_deciles(series: pd.Series) -> pd.Series:
    """
    Assign decile labels 0–9 (D1=0, D10=9) using pd.qcut on PearsonR.
    Returns integer decile index starting at 0.
    """
    labels = list(range(N_DECILES))
    decile, _ = pd.qcut(series, q=N_DECILES, labels=labels, retbins=True, duplicates="drop")
    return decile.astype(int)


def bootstrap_ci(data: np.ndarray, n_boot: int = 2000, seed: int = 42,
                 alpha: float = 0.05) -> tuple[float, float]:
    """Bootstrap 95% CI for the mean of *data* (proportion 0/1 for eGene rate)."""
    rng = np.random.default_rng(seed)
    boot_means = np.array([
        rng.choice(data, size=len(data), replace=True).mean()
        for _ in range(n_boot)
    ])
    lo = np.percentile(boot_means, 100 * alpha / 2)
    hi = np.percentile(boot_means, 100 * (1 - alpha / 2))
    return float(lo), float(hi)


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 60)
    print("03_decile_analysis.py  —  Decile analysis")
    print("=" * 60)

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} genes, {df['is_egene'].sum()} eGenes")

    # Drop rows with missing PearsonR
    df = df.dropna(subset=["PearsonR"]).reset_index(drop=True)
    print(f"After dropping missing PearsonR: {len(df)} genes")

    # Assign deciles (0-based index: 0 = D1 = lowest R)
    df["decile"] = assign_deciles(df["PearsonR"])

    # -----------------------------------------------------------------------
    # Per-decile summary with bootstrap CIs
    # -----------------------------------------------------------------------
    print("\nComputing per-decile statistics and bootstrap CIs ...")
    records = []
    for d in range(N_DECILES):
        sub = df[df["decile"] == d]
        n = len(sub)
        n_egene = int(sub["is_egene"].sum())
        egene_rate = n_egene / n if n > 0 else np.nan
        mean_r = sub["PearsonR"].mean()
        median_r = sub["PearsonR"].median()
        ci_lo, ci_hi = bootstrap_ci(sub["is_egene"].values, n_boot=N_BOOT, seed=BOOT_SEED + d)
        records.append({
            "decile": d + 1,          # 1-based for reporting
            "decile_idx": d,          # 0-based internal
            "n": n,
            "n_egene": n_egene,
            "egene_rate": egene_rate,
            "ci_lo": ci_lo,
            "ci_hi": ci_hi,
            "mean_R": mean_r,
            "median_R": median_r,
        })

    decile_df = pd.DataFrame(records)

    overall_rate = df["is_egene"].mean()

    print("\nDecile Summary Table")
    print("-" * 80)
    header = (
        f"{'D':>3} {'n':>6} {'n_eG':>6} {'eGene%':>8} "
        f"{'CI_lo':>7} {'CI_hi':>7} {'mean_R':>8} {'median_R':>9}"
    )
    print(header)
    print("-" * 80)
    for _, row in decile_df.iterrows():
        print(
            f"D{int(row['decile']):>2} "
            f"{int(row['n']):>6} "
            f"{int(row['n_egene']):>6} "
            f"{row['egene_rate']*100:>7.1f}% "
            f"{row['ci_lo']*100:>7.1f} "
            f"{row['ci_hi']*100:>7.1f} "
            f"{row['mean_R']:>8.4f} "
            f"{row['median_R']:>9.4f}"
        )
    print("-" * 80)
    print(f"Overall eGene rate: {overall_rate*100:.1f}%")

    # -----------------------------------------------------------------------
    # Logistic regression
    # -----------------------------------------------------------------------
    print("\nFitting logistic regression ...")

    reg_df = df.dropna(subset=["PearsonR", "MeanObs", "StdObs", "log10_gene_length"]).copy()

    reg_df["PearsonR2"] = reg_df["PearsonR"] ** 2
    reg_df["MeanObs_z"] = zscore(reg_df["MeanObs"])
    reg_df["StdObs_z"] = zscore(reg_df["StdObs"])
    reg_df["log10_gl_z"] = zscore(reg_df["log10_gene_length"])

    X = sm.add_constant(
        reg_df[["PearsonR", "PearsonR2", "MeanObs_z", "StdObs_z", "log10_gl_z"]]
    )
    y = reg_df["is_egene"]

    logit_model = sm.Logit(y, X)
    logit_result = logit_model.fit(maxiter=200, disp=False)

    print("\nLogistic Regression Results")
    print("-" * 60)
    print(logit_result.summary2())

    print("\nKey coefficients:")
    for coef_name in ["PearsonR", "PearsonR2"]:
        idx = list(logit_result.params.index).index(coef_name)
        coef = logit_result.params[coef_name]
        pval = logit_result.pvalues[coef_name]
        ci_lo, ci_hi = logit_result.conf_int().loc[coef_name]
        print(f"  {coef_name:15s}: coef={coef:+.4f}, p={pval:.3e}, 95%CI [{ci_lo:.4f}, {ci_hi:.4f}]")

    r2_pval = logit_result.pvalues.get("PearsonR2", np.nan)
    print(f"\n  PearsonR² term p-value: {r2_pval:.3e}")

    # -----------------------------------------------------------------------
    # Class A (D10) vs Class B (D1) eGene table
    # -----------------------------------------------------------------------
    print("\nBuilding Class A / Class B comparison table ...")

    class_a = df[(df["decile"] == 9) & (df["is_egene"] == 1)].copy()  # D10
    class_b = df[(df["decile"] == 0) & (df["is_egene"] == 1)].copy()  # D1
    class_a["class"] = "A (D10)"
    class_b["class"] = "B (D1)"

    print(f"  Class A (D10 eGenes): n={len(class_a)}")
    print(f"  Class B (D1  eGenes): n={len(class_b)}")

    def class_stats(sub: pd.DataFrame, label: str) -> dict:
        row = {"class": label, "n": len(sub)}
        for col in ["PearsonR", "PrediXcanR", "MeanObs", "StdObs", "PearsonRfineTuned"]:
            if col in sub.columns:
                vals = sub[col].dropna()
                row[f"{col}_mean"] = vals.mean()
                row[f"{col}_median"] = vals.median()
                row[f"{col}_std"] = vals.std()
        return row

    stats_a = class_stats(class_a, "Class A (D10)")
    stats_b = class_stats(class_b, "Class B (D1)")

    # Mann-Whitney test on PrediXcan R
    a_pred = class_a["PrediXcanR"].dropna().values
    b_pred = class_b["PrediXcanR"].dropna().values
    if len(a_pred) > 0 and len(b_pred) > 0:
        mw_stat, mw_pval = stats.mannwhitneyu(a_pred, b_pred, alternative="two-sided")
        print(f"\n  PrediXcan R Mann-Whitney: U={mw_stat:.1f}, p={mw_pval:.3e}")
        print(f"  Class A median PrediXcanR: {np.median(a_pred):.3f}")
        print(f"  Class B median PrediXcanR: {np.median(b_pred):.3f}")

    table_df = pd.DataFrame([stats_a, stats_b])
    table_df.to_csv(CLASS_TABLE_OUT, index=False)
    print(f"\nSaved class comparison table: {CLASS_TABLE_OUT}")

    # -----------------------------------------------------------------------
    # Save decile summary (add decile labels back to df for downstream)
    # -----------------------------------------------------------------------
    save_cols = ["decile", "n", "n_egene", "egene_rate", "ci_lo", "ci_hi", "mean_R", "median_R"]
    decile_df[save_cols].to_csv(DECILE_OUT, index=False)
    print(f"Saved decile summary: {DECILE_OUT}")

    # Also save decile column back to gene-level data for downstream scripts
    df_with_decile = df.copy()
    df_with_decile.to_csv(INPUT_FILE, index=False)
    print(f"Updated gene_level_results.csv with decile column.")

    print("\nDone.")


if __name__ == "__main__":
    main()
