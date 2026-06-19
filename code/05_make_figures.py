"""
05_make_figures.py
------------------
Generate all publication-quality figures for the eQTL / Enformer U-shape paper.

Figures:
  fig1_egene_decile_ushape.png/.pdf
    Bar chart of per-decile eGene rates with bootstrap CI error bars,
    overall mean dashed line, cubic spline overlay, D1=red, D10=blue, rest=gray.

  fig2_spline_logistic.png/.pdf
    Panel A: Logistic regression predicted probability vs. PearsonR
    Panel B: Scatter of decile mean_R vs. egene_rate with quadratic fit

  fig3_predixcan_classA_B.png/.pdf
    Panel A: Violin/box comparison of PrediXcan R, Class A vs. B
    Panel B: Fine-tuning response (PearsonRfineTuned - PearsonR) bar chart
    Panel C: Expression stats (MeanObs, StdObs) comparison

Outputs: manuscript/figures/
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.interpolate import make_interp_spline
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
PROC_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
FIG_DIR = os.path.join(PROJECT_ROOT, "manuscript", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

GENE_FILE = os.path.join(PROC_DIR, "gene_level_results.csv")
DECILE_FILE = os.path.join(PROC_DIR, "decile_summary.csv")
CLASSA_FILE = os.path.join(PROC_DIR, "classA_D10_eGenes.csv")
CLASSB_FILE = os.path.join(PROC_DIR, "classB_D1_eGenes.csv")

# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "lines.linewidth": 1.5,
    "patch.linewidth": 0.8,
})

# Color scheme
COLOR_D1 = "#D62728"    # red — D1 (worst Enformer performance)
COLOR_D10 = "#1F77B4"   # blue — D10 (best Enformer performance)
COLOR_GRAY = "#AAAAAA"  # intermediate deciles
COLOR_SPLINE = "#2CA02C"  # green spline overlay
COLOR_OVERALL = "#555555"


def assign_deciles(series: pd.Series, n: int = 10) -> pd.Series:
    labels = list(range(n))
    decile, _ = pd.qcut(series, q=n, labels=labels, retbins=True, duplicates="drop")
    return decile.astype(int)


def save_fig(fig: plt.Figure, name: str) -> None:
    for ext in ("png", "pdf"):
        path = os.path.join(FIG_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_data():
    df = pd.read_csv(GENE_FILE)
    df = df.dropna(subset=["PearsonR"]).reset_index(drop=True)
    if "decile" not in df.columns:
        df["decile"] = assign_deciles(df["PearsonR"])

    decile_df = pd.read_csv(DECILE_FILE)

    class_a = pd.read_csv(CLASSA_FILE) if os.path.exists(CLASSA_FILE) else pd.DataFrame()
    class_b = pd.read_csv(CLASSB_FILE) if os.path.exists(CLASSB_FILE) else pd.DataFrame()

    return df, decile_df, class_a, class_b


# ---------------------------------------------------------------------------
# FIG 1: U-shaped bar chart
# ---------------------------------------------------------------------------
def make_fig1(decile_df: pd.DataFrame, overall_rate: float) -> None:
    print("\nGenerating fig1_egene_decile_ushape ...")

    fig, ax = plt.subplots(figsize=(8, 5))

    deciles = decile_df["decile"].values          # 1-based
    rates = decile_df["egene_rate"].values * 100  # percent
    ci_lo = decile_df["ci_lo"].values * 100
    ci_hi = decile_df["ci_hi"].values * 100

    # Error bars: asymmetric (rate - lo, hi - rate)
    yerr_lo = rates - ci_lo
    yerr_hi = ci_hi - rates

    # Bar colors
    bar_colors = []
    for d in deciles:
        if d == 1:
            bar_colors.append(COLOR_D1)
        elif d == 10:
            bar_colors.append(COLOR_D10)
        else:
            bar_colors.append(COLOR_GRAY)

    x = np.arange(len(deciles))
    bars = ax.bar(
        x, rates, color=bar_colors, width=0.7, zorder=3,
        edgecolor="white", linewidth=0.5,
    )
    ax.errorbar(
        x, rates,
        yerr=[yerr_lo, yerr_hi],
        fmt="none", color="black", capsize=4, linewidth=1.2, zorder=4,
    )

    # Overall mean dashed line
    ax.axhline(overall_rate * 100, color=COLOR_OVERALL, linestyle="--",
               linewidth=1.5, zorder=2, label=f"Overall mean ({overall_rate*100:.1f}%)")

    # Smooth spline overlay (fit on 1..10 then interpolate)
    x_pts = np.array(deciles, dtype=float)
    y_pts = rates
    x_smooth = np.linspace(x_pts.min(), x_pts.max(), 300)
    if len(x_pts) >= 4:
        k = min(3, len(x_pts) - 1)
        spl = make_interp_spline(x_pts, y_pts, k=k)
        y_smooth = spl(x_smooth)
        # Map decile 1..10 to bar positions 0..9
        ax.plot(x_smooth - 1, y_smooth, color=COLOR_SPLINE, linewidth=2,
                zorder=5, label="Spline", alpha=0.85)

    # Axes
    ax.set_xticks(x)
    ax.set_xticklabels([f"D{d}" for d in deciles])
    ax.set_xlabel("Enformer Performance Decile (PearsonR)", labelpad=8)
    ax.set_ylabel("Brain Cortex eGene Rate (%)")
    ax.set_title("eGene Rate by Enformer Performance Decile", pad=12)
    ax.set_ylim(0, max(rates) * 1.3)
    ax.yaxis.grid(True, linestyle=":", alpha=0.5, zorder=1)
    ax.set_axisbelow(True)

    # Legend
    patch_d1 = mpatches.Patch(color=COLOR_D1, label="D1 (lowest R)")
    patch_d10 = mpatches.Patch(color=COLOR_D10, label="D10 (highest R)")
    patch_mid = mpatches.Patch(color=COLOR_GRAY, label="D2–D9")
    line_overall = plt.Line2D([0], [0], color=COLOR_OVERALL, linestyle="--",
                               linewidth=1.5, label=f"Overall mean ({overall_rate*100:.1f}%)")
    line_spline = plt.Line2D([0], [0], color=COLOR_SPLINE, linewidth=2, label="Spline")
    ax.legend(handles=[patch_d1, patch_d10, patch_mid, line_overall, line_spline],
              loc="upper center", frameon=False, ncol=3, fontsize=9)

    # Annotate D1 and D10 bars
    d1_idx, d10_idx = 0, len(deciles) - 1
    ax.text(d1_idx, rates[d1_idx] + yerr_hi[d1_idx] + 1.0,
            f"{rates[d1_idx]:.1f}%", ha="center", va="bottom", color=COLOR_D1,
            fontsize=9, fontweight="bold")
    ax.text(d10_idx, rates[d10_idx] + yerr_hi[d10_idx] + 1.0,
            f"{rates[d10_idx]:.1f}%", ha="center", va="bottom", color=COLOR_D10,
            fontsize=9, fontweight="bold")

    fig.tight_layout()
    save_fig(fig, "fig1_egene_decile_ushape")
    plt.close(fig)


# ---------------------------------------------------------------------------
# FIG 2: Logistic regression predicted probability + quadratic decile fit
# ---------------------------------------------------------------------------
def make_fig2(df: pd.DataFrame, decile_df: pd.DataFrame) -> None:
    print("\nGenerating fig2_spline_logistic ...")

    # Fit logistic regression for panel A
    reg = df.dropna(subset=["PearsonR", "MeanObs", "StdObs", "log10_gene_length"]).copy()
    reg["PearsonR2"] = reg["PearsonR"] ** 2
    reg["MeanObs_z"] = (reg["MeanObs"] - reg["MeanObs"].mean()) / reg["MeanObs"].std()
    reg["StdObs_z"] = (reg["StdObs"] - reg["StdObs"].mean()) / reg["StdObs"].std()
    reg["log10_gl_z"] = (
        (reg["log10_gene_length"] - reg["log10_gene_length"].mean())
        / reg["log10_gene_length"].std()
    )

    X = sm.add_constant(reg[["PearsonR", "PearsonR2", "MeanObs_z", "StdObs_z", "log10_gl_z"]])
    y = reg["is_egene"]
    logit_result = sm.Logit(y, X).fit(maxiter=200, disp=False)

    # Prediction grid (vary PearsonR, hold covariates at zero/mean)
    r_range = np.linspace(reg["PearsonR"].min(), reg["PearsonR"].max(), 400)
    X_pred = pd.DataFrame({
        "const": 1.0,
        "PearsonR": r_range,
        "PearsonR2": r_range ** 2,
        "MeanObs_z": 0.0,
        "StdObs_z": 0.0,
        "log10_gl_z": 0.0,
    })
    y_pred = logit_result.predict(X_pred)

    # Quadratic fit for panel B
    mean_r = decile_df["mean_R"].values
    egene_rate = decile_df["egene_rate"].values * 100
    quad_coef = np.polyfit(mean_r, egene_rate, 2)
    quad_fit = np.polyval(quad_coef, np.sort(mean_r))

    # -------- Plot --------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Logistic predicted probability
    ax = axes[0]
    # Rug plot of individual genes (subsample for clarity)
    sample_1 = reg[reg["is_egene"] == 1]["PearsonR"].values
    sample_0 = reg[reg["is_egene"] == 0]["PearsonR"].values
    rng = np.random.default_rng(42)
    idx_sub = rng.choice(len(sample_0), size=min(3000, len(sample_0)), replace=False)
    ax.scatter(sample_0[idx_sub], rng.uniform(-0.02, 0.02, size=len(idx_sub)),
               s=1, alpha=0.08, color=COLOR_GRAY, label="Non-eGene", rasterized=True)
    idx_sub1 = rng.choice(len(sample_1), size=min(1000, len(sample_1)), replace=False)
    ax.scatter(sample_1[idx_sub1], rng.uniform(0.98, 1.02, size=len(idx_sub1)),
               s=1, alpha=0.15, color=COLOR_D10, label="eGene", rasterized=True)
    # Logistic curve
    ax.plot(r_range, y_pred, color="#D62728", linewidth=2.5, zorder=5,
            label="Logistic fit")

    r2_pval = logit_result.pvalues.get("PearsonR2", np.nan)
    ax.set_xlabel("Enformer Pearson R")
    ax.set_ylabel("Predicted eGene Probability")
    ax.set_title("A   Logistic Regression Predicted Probability", loc="left")
    ax.text(0.03, 0.93, f"R² term p = {r2_pval:.2e}", transform=ax.transAxes,
            fontsize=9, color="#333333")
    ax.set_ylim(-0.08, 1.08)
    ax.legend(loc="center right", frameon=False, markerscale=6, fontsize=9)

    # Panel B: Decile mean_R vs. egene_rate with quadratic fit
    ax = axes[1]
    colors_dec = [COLOR_D1 if i == 0 else (COLOR_D10 if i == 9 else COLOR_GRAY)
                  for i in range(len(mean_r))]
    ax.scatter(mean_r, egene_rate, c=colors_dec, s=70, zorder=4, edgecolors="white",
               linewidth=0.5)

    # Quadratic fit line
    r_sort = np.sort(mean_r)
    ax.plot(r_sort, np.polyval(quad_coef, r_sort), "--", color=COLOR_SPLINE,
            linewidth=2, label="Quadratic fit")

    # Label D1 and D10
    ax.annotate("D1", (mean_r[0], egene_rate[0]), textcoords="offset points",
                xytext=(-18, 5), color=COLOR_D1, fontsize=9, fontweight="bold")
    ax.annotate("D10", (mean_r[-1], egene_rate[-1]), textcoords="offset points",
                xytext=(5, 5), color=COLOR_D10, fontsize=9, fontweight="bold")

    ax.set_xlabel("Decile Mean Pearson R")
    ax.set_ylabel("eGene Rate (%)")
    ax.set_title("B   Decile Mean R vs. eGene Rate", loc="left")
    ax.legend(frameon=False, fontsize=9)

    for ax in axes:
        ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
        ax.set_axisbelow(True)

    fig.tight_layout(pad=2.0)
    save_fig(fig, "fig2_spline_logistic")
    plt.close(fig)


# ---------------------------------------------------------------------------
# FIG 3: PrediXcan R, fine-tuning, expression stats comparison
# ---------------------------------------------------------------------------
def make_fig3(class_a: pd.DataFrame, class_b: pd.DataFrame) -> None:
    print("\nGenerating fig3_predixcan_classA_B ...")

    if class_a.empty or class_b.empty:
        print("  Class A or B data missing — skipping fig3")
        return

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # ---- Panel A: Violin + box for PrediXcanR ----
    ax = axes[0]
    a_pred = class_a["PrediXcanR"].dropna().values
    b_pred = class_b["PrediXcanR"].dropna().values

    parts = ax.violinplot(
        [a_pred, b_pred],
        positions=[1, 2],
        showmedians=False,
        showextrema=False,
    )
    colors_vio = [COLOR_D10, COLOR_D1]
    for i, (pc, c) in enumerate(zip(parts["bodies"], colors_vio)):
        pc.set_facecolor(c)
        pc.set_alpha(0.55)
        pc.set_edgecolor("white")

    # Box overlay
    bp = ax.boxplot(
        [a_pred, b_pred],
        positions=[1, 2],
        widths=0.12,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(color="gray"),
        capprops=dict(color="gray"),
        flierprops=dict(marker=".", markersize=3, alpha=0.4, color="gray"),
    )
    for patch, c in zip(bp["boxes"], colors_vio):
        patch.set_facecolor("white")
        patch.set_edgecolor(c)
        patch.set_linewidth(1.5)

    # Mann-Whitney p
    if len(a_pred) > 1 and len(b_pred) > 1:
        _, pval = stats.mannwhitneyu(a_pred, b_pred, alternative="two-sided")
        y_top = max(np.percentile(a_pred, 99), np.percentile(b_pred, 99))
        ax.plot([1, 2], [y_top * 1.05, y_top * 1.05], "k-", linewidth=0.8)
        sig = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
        ax.text(1.5, y_top * 1.07, f"{sig}\np={pval:.2e}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Class A\n(D10 eGenes)", "Class B\n(D1 eGenes)"])
    ax.set_ylabel("PrediXcan R")
    ax.set_title("A   PrediXcan R Comparison", loc="left")
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    # ---- Panel B: Fine-tuning delta bar chart ----
    ax = axes[1]

    delta_data = {}
    for label, sub, color in [("Class A\n(D10)", class_a, COLOR_D10),
                                ("Class B\n(D1)", class_b, COLOR_D1)]:
        if "PearsonRfineTuned" in sub.columns:
            delta = (sub["PearsonRfineTuned"] - sub["PearsonR"]).dropna()
            delta_data[label] = (delta.mean(), delta.sem(), color)

    if delta_data:
        xlabels = list(delta_data.keys())
        means = [v[0] for v in delta_data.values()]
        sems = [v[1] for v in delta_data.values()]
        colors_bar = [v[2] for v in delta_data.values()]

        x = np.arange(len(xlabels))
        bars = ax.bar(x, means, yerr=sems, color=colors_bar, width=0.5,
                      capsize=5, error_kw={"linewidth": 1.2},
                      edgecolor="white", zorder=3)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="-")
        ax.set_xticks(x)
        ax.set_xticklabels(xlabels)
        ax.set_ylabel("ΔPearsonR (fine-tuned − base)")
        ax.set_title("B   Fine-Tuning Response", loc="left")

        # Mann-Whitney on delta if both available
        if len(delta_data) == 2:
            keys = list(delta_data.keys())
            d_a = (class_a["PearsonRfineTuned"] - class_a["PearsonR"]).dropna().values
            d_b = (class_b["PearsonRfineTuned"] - class_b["PearsonR"]).dropna().values
            if len(d_a) > 1 and len(d_b) > 1:
                _, pval_ft = stats.mannwhitneyu(d_a, d_b, alternative="two-sided")
                ax.text(0.5, 0.92, f"p={pval_ft:.2e}", transform=ax.transAxes,
                        ha="center", fontsize=9, color="#333333")
    else:
        ax.text(0.5, 0.5, "PearsonRfineTuned\nnot available",
                ha="center", va="center", transform=ax.transAxes, fontsize=10,
                color="#888888")
        ax.set_title("B   Fine-Tuning Response", loc="left")

    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    # ---- Panel C: MeanObs / StdObs comparison ----
    ax = axes[2]

    stats_rows = []
    for col, label_col in [("MeanObs", "Mean Expression"), ("StdObs", "Std Expression")]:
        if col not in class_a.columns:
            continue
        a_vals = class_a[col].dropna().values
        b_vals = class_b[col].dropna().values
        _, pval_mw = stats.mannwhitneyu(a_vals, b_vals, alternative="two-sided") \
            if len(a_vals) > 1 and len(b_vals) > 1 else (np.nan, np.nan)
        stats_rows.append({
            "metric": label_col,
            "A_mean": np.mean(a_vals), "A_sem": stats.sem(a_vals),
            "B_mean": np.mean(b_vals), "B_sem": stats.sem(b_vals),
            "pval": pval_mw,
        })

    if stats_rows:
        width = 0.3
        x = np.arange(len(stats_rows))
        ax.bar(x - width / 2,
               [r["A_mean"] for r in stats_rows],
               width, yerr=[r["A_sem"] for r in stats_rows],
               color=COLOR_D10, label="Class A (D10)", capsize=4,
               edgecolor="white", zorder=3)
        ax.bar(x + width / 2,
               [r["B_mean"] for r in stats_rows],
               width, yerr=[r["B_sem"] for r in stats_rows],
               color=COLOR_D1, label="Class B (D1)", capsize=4,
               edgecolor="white", zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([r["metric"] for r in stats_rows])
        ax.set_ylabel("Value")
        ax.set_title("C   Expression Statistics", loc="left")
        ax.legend(frameon=False, fontsize=9)

        # Significance annotations
        for i, row in enumerate(stats_rows):
            if not np.isnan(row["pval"]):
                sig = "***" if row["pval"] < 0.001 else \
                      ("**" if row["pval"] < 0.01 else ("*" if row["pval"] < 0.05 else "ns"))
                y_top = max(row["A_mean"] + row["A_sem"], row["B_mean"] + row["B_sem"]) * 1.08
                ax.text(i, y_top, sig, ha="center", va="bottom", fontsize=10)
    else:
        ax.text(0.5, 0.5, "MeanObs/StdObs\nnot available",
                ha="center", va="center", transform=ax.transAxes, fontsize=10,
                color="#888888")
        ax.set_title("C   Expression Statistics", loc="left")

    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    fig.tight_layout(pad=2.0)
    save_fig(fig, "fig3_predixcan_classA_B")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 60)
    print("05_make_figures.py  —  Generating publication figures")
    print("=" * 60)

    df, decile_df, class_a, class_b = load_data()
    overall_rate = df["is_egene"].mean()

    print(f"Gene data: {len(df)} genes, {df['is_egene'].sum()} eGenes")
    print(f"Overall eGene rate: {overall_rate*100:.1f}%")
    print(f"Decile summary: {len(decile_df)} rows")
    print(f"Class A: {len(class_a)} genes | Class B: {len(class_b)} genes")

    make_fig1(decile_df, overall_rate)
    make_fig2(df, decile_df)
    make_fig3(class_a, class_b)

    print(f"\nAll figures saved to: {FIG_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()
