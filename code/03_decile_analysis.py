"""
03_decile_analysis.py
---------------------
Decile analysis of Enformer per-gene |Pearson R| vs. GTEx Brain Cortex eGene status.

Key analyses:
  - Assign absR_decile (D1–D10) by |PearsonR| via pd.qcut
  - Bootstrap 95% CIs on per-decile eGene rates (2000 resamples)
  - Per-decile direction accuracy (fraction of eGenes with R > 0)
  - 4-model AIC comparison to formally test magnitude vs. direction:
        M1: is_egene ~ R + R2 + covariates  (U-shape model, AIC ~8378)
        M2: is_egene ~ |R| + covariates     (magnitude model, AIC ~8274)
        M3: is_egene ~ |R| + R + covariates
        M4: is_egene ~ |R| + R + R2 + covariates
  - Save summary tables and AIC comparison

Outputs:
  data/processed/decile_summary.csv        (per |R| decile: n, eGene%, CI)
  data/processed/direction_by_decile.csv   (per |R| decile: direction accuracy)
  data/processed/aic_model_comparison.csv  (4-model AIC comparison)
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
DIRECTION_OUT = os.path.join(PROC_DIR, "direction_by_decile.csv")
AIC_OUT = os.path.join(PROC_DIR, "aic_model_comparison.csv")

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
    Assign decile labels 0-9 (D1=0, D10=9) using pd.qcut.
    Returns integer decile index (0-based).
    """
    labels = list(range(N_DECILES))
    decile, _ = pd.qcut(series, q=N_DECILES, labels=labels, retbins=True, duplicates="drop")
    return decile.astype(int)


def bootstrap_ci(data: np.ndarray, n_boot: int = 2000, seed: int = 42,
                 alpha: float = 0.05) -> tuple:
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
    print("03_decile_analysis.py  --  |R| Decile and AIC analysis")
    print("=" * 60)

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} genes, {int(df['is_egene'].sum())} eGenes")

    df = df.dropna(subset=["PearsonR"]).reset_index(drop=True)
    print(f"After dropping missing PearsonR: {len(df)} genes")

    # Compute |R| and assign |R|-based deciles (the primary analysis variable).
    # absR_decile: 0 = lowest |R| (D1), 9 = highest |R| (D10)
    df["abs_R"] = df["PearsonR"].abs()
    df["absR_decile"] = assign_deciles(df["abs_R"])

    # Binary direction indicator for downstream analysis
    df["R_positive"] = (df["PearsonR"] > 0).astype(int)

    # -----------------------------------------------------------------------
    # Per |R| decile: eGene rate + bootstrap CIs + direction accuracy
    # -----------------------------------------------------------------------
    print("\nComputing per-|R|-decile statistics ...")
    decile_records = []
    direction_records = []

    for d in range(N_DECILES):
        sub = df[df["absR_decile"] == d]
        n = len(sub)
        n_egene = int(sub["is_egene"].sum())
        egene_rate = n_egene / n if n > 0 else np.nan
        mean_absR = sub["abs_R"].mean()
        median_absR = sub["abs_R"].median()
        ci_lo, ci_hi = bootstrap_ci(sub["is_egene"].values, n_boot=N_BOOT, seed=BOOT_SEED + d)

        # Direction accuracy among eGenes in this decile
        sub_eg = sub[sub["is_egene"] == 1]
        n_eg = len(sub_eg)
        if n_eg > 0:
            n_rpos = int(sub_eg["R_positive"].sum())
            dir_acc = n_rpos / n_eg
            binom_result = stats.binomtest(n_rpos, n_eg, p=0.5, alternative="greater")
            dir_pval = binom_result.pvalue
        else:
            dir_acc = np.nan
            dir_pval = np.nan

        decile_records.append({
            "absR_decile": d + 1,
            "absR_decile_idx": d,
            "n": n,
            "n_egene": n_egene,
            "egene_rate": egene_rate,
            "ci_lo": ci_lo,
            "ci_hi": ci_hi,
            "mean_absR": mean_absR,
            "median_absR": median_absR,
        })
        direction_records.append({
            "absR_decile": d + 1,
            "n_egene": n_eg,
            "n_R_positive": int(sub_eg["R_positive"].sum()) if n_eg > 0 else 0,
            "direction_accuracy": dir_acc,
            "direction_pval": dir_pval,
        })

    decile_df = pd.DataFrame(decile_records)
    direction_df = pd.DataFrame(direction_records)
    overall_rate = df["is_egene"].mean()

    print("\n|R| Decile Summary")
    print("-" * 80)
    print(f"{'D':>3} {'n':>6} {'n_eG':>6} {'eGene%':>8} "
          f"{'CI_lo':>7} {'CI_hi':>7} {'mean|R|':>8} {'DirAcc%':>8}")
    print("-" * 80)
    for i, row in decile_df.iterrows():
        d_idx = int(row["absR_decile"]) - 1
        dr = direction_records[d_idx]
        dir_str = f"{dr['direction_accuracy']*100:.1f}" if not np.isnan(dr['direction_accuracy']) else "  n/a"
        print(
            f"D{int(row['absR_decile']):>2} "
            f"{int(row['n']):>6} "
            f"{int(row['n_egene']):>6} "
            f"{row['egene_rate']*100:>7.1f}% "
            f"{row['ci_lo']*100:>7.1f} "
            f"{row['ci_hi']*100:>7.1f} "
            f"{row['mean_absR']:>8.4f} "
            f"{dir_str:>8}%"
        )
    print("-" * 80)
    print(f"Overall eGene rate: {overall_rate*100:.1f}%")
    print(f"D1 eGene rate: {decile_df.iloc[0]['egene_rate']*100:.1f}%")
    print(f"D10 eGene rate: {decile_df.iloc[-1]['egene_rate']*100:.1f}%")

    # Monotonicity test
    spearman_r, spearman_p = stats.spearmanr(
        decile_df["absR_decile"], decile_df["egene_rate"]
    )
    print(f"\nMonotonicity (Spearman rho): {spearman_r:.4f}, p={spearman_p:.3e}")

    # -----------------------------------------------------------------------
    # 2x2 Factorial: |R| stratum x sign(R) -> eGene rate
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("2x2 Factorial: |R| stratum x sign(R) -> eGene rate")
    print("=" * 60)

    high_absR = df[df["absR_decile"] == 9]
    low_absR  = df[df["absR_decile"] == 0]

    cells = [
        ("High|R|, R>0", high_absR[high_absR["PearsonR"] > 0]),
        ("High|R|, R<0", high_absR[high_absR["PearsonR"] < 0]),
        ("Low|R|,  R>0", low_absR[low_absR["PearsonR"] > 0]),
        ("Low|R|,  R<0", low_absR[low_absR["PearsonR"] < 0]),
    ]
    for label, sub in cells:
        n_sub = len(sub)
        n_eg = int(sub["is_egene"].sum())
        rate = n_eg / n_sub * 100 if n_sub > 0 else np.nan
        print(f"  {label}: n={n_sub}, eGenes={n_eg}, rate={rate:.1f}%")

    # Chi-squared test within high |R| stratum (key test)
    hi_pos = high_absR[high_absR["PearsonR"] > 0]
    hi_neg = high_absR[high_absR["PearsonR"] < 0]
    ct = np.array([
        [int((hi_pos["is_egene"] == 1).sum()), int((hi_pos["is_egene"] == 0).sum())],
        [int((hi_neg["is_egene"] == 1).sum()), int((hi_neg["is_egene"] == 0).sum())],
    ])
    chi2, chi2_p, _, _ = stats.chi2_contingency(ct)
    print(f"\n  Chi-squared (High|R|: R>0 vs R<0): chi2={chi2:.3f}, p={chi2_p:.3f}")

    # Overall direction for eGenes
    egenes = df[df["is_egene"] == 1]
    n_rpos_eg = int(egenes["R_positive"].sum())
    dir_overall = n_rpos_eg / len(egenes)
    binom_res = stats.binomtest(n_rpos_eg, len(egenes), p=0.5, alternative="two-sided")
    print(f"\n  eGene overall direction accuracy (R>0): {dir_overall*100:.1f}%")
    print(f"  Binomial test (vs 50%): p={binom_res.pvalue:.3e}")
    non_egenes = df[df["is_egene"] == 0]
    print(f"  Non-eGene R>0 fraction: {non_egenes['R_positive'].mean()*100:.1f}%")

    # -----------------------------------------------------------------------
    # 4-model AIC comparison
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("4-Model AIC Comparison")
    print("=" * 60)

    reg_df = df.dropna(subset=["PearsonR", "MeanObs", "StdObs", "log10_gene_length"]).copy()
    reg_df["PearsonR2"] = reg_df["PearsonR"] ** 2
    reg_df["MeanObs_z"] = zscore(reg_df["MeanObs"])
    reg_df["StdObs_z"] = zscore(reg_df["StdObs"])
    reg_df["log10_gl_z"] = zscore(reg_df["log10_gene_length"])
    y_reg = reg_df["is_egene"]
    covs = ["MeanObs_z", "StdObs_z", "log10_gl_z"]

    # M1: R + R2 + covariates (U-shape model)
    X_m1 = sm.add_constant(reg_df[["PearsonR", "PearsonR2"] + covs])
    res_m1 = sm.Logit(y_reg, X_m1).fit(maxiter=200, disp=False)

    # M2: |R| + covariates (magnitude model)
    X_m2 = sm.add_constant(reg_df[["abs_R"] + covs])
    res_m2 = sm.Logit(y_reg, X_m2).fit(maxiter=200, disp=False)

    # M3: |R| + R + covariates
    X_m3 = sm.add_constant(reg_df[["abs_R", "PearsonR"] + covs])
    res_m3 = sm.Logit(y_reg, X_m3).fit(maxiter=200, disp=False)

    # M4: |R| + R + R2 + covariates (full model)
    X_m4 = sm.add_constant(reg_df[["abs_R", "PearsonR", "PearsonR2"] + covs])
    res_m4 = sm.Logit(y_reg, X_m4).fit(maxiter=200, disp=False)

    aic_rows = [
        {"model": "M1: R + R2 + covariates (U-shape)",    "n_params": int(res_m1.df_model + 1), "AIC": round(res_m1.aic, 1), "BIC": round(res_m1.bic, 1)},
        {"model": "M2: |R| + covariates (magnitude)",      "n_params": int(res_m2.df_model + 1), "AIC": round(res_m2.aic, 1), "BIC": round(res_m2.bic, 1)},
        {"model": "M3: |R| + R + covariates",              "n_params": int(res_m3.df_model + 1), "AIC": round(res_m3.aic, 1), "BIC": round(res_m3.bic, 1)},
        {"model": "M4: |R| + R + R2 + covariates (full)", "n_params": int(res_m4.df_model + 1), "AIC": round(res_m4.aic, 1), "BIC": round(res_m4.bic, 1)},
    ]
    aic_df = pd.DataFrame(aic_rows)
    best_aic = aic_df["AIC"].min()
    aic_df["delta_AIC"] = (aic_df["AIC"] - best_aic).round(1)

    print("\n  AIC Model Comparison")
    print("  " + "-" * 72)
    print(f"  {'Model':<44} {'k':>4} {'AIC':>8} {'BIC':>8} {'dAIC':>7}")
    print("  " + "-" * 72)
    for _, row in aic_df.iterrows():
        print(f"  {row['model']:<44} {int(row['n_params']):>4} "
              f"{row['AIC']:>8.1f} {row['BIC']:>8.1f} {row['delta_AIC']:>7.1f}")
    print("  " + "-" * 72)
    best_model = aic_df.loc[aic_df["AIC"].idxmin(), "model"]
    print(f"\n  Best model: {best_model}")
    print(f"  M1 AIC = {res_m1.aic:.1f}")
    print(f"  M2 AIC = {res_m2.aic:.1f}")
    print(f"  Delta AIC (M2 - M1) = {res_m2.aic - res_m1.aic:.1f}  "
          f"(negative = M2 better; magnitude model preferred)")

    r_coef_m3 = res_m3.params.get("PearsonR", np.nan)
    r_pval_m3 = res_m3.pvalues.get("PearsonR", np.nan)
    r2_pval_m4 = res_m4.pvalues.get("PearsonR2", np.nan)
    print(f"\n  M3 -- R coefficient after |R| controlled: "
          f"coef={r_coef_m3:+.4f}, p={r_pval_m3:.3e}")
    print(f"  M4 -- R2 coefficient after |R|+R: p={r2_pval_m4:.3e}")

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------
    decile_df.to_csv(DECILE_OUT, index=False)
    print(f"\nSaved |R| decile summary: {DECILE_OUT}")

    direction_df.to_csv(DIRECTION_OUT, index=False)
    print(f"Saved direction accuracy: {DIRECTION_OUT}")

    aic_df.to_csv(AIC_OUT, index=False)
    print(f"Saved AIC comparison: {AIC_OUT}")

    df.to_csv(INPUT_FILE, index=False)
    print(f"Updated gene_level_results.csv with abs_R, absR_decile, R_positive.")

    print("\nDone.")


if __name__ == "__main__":
    main()
