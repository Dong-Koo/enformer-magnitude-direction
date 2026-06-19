"""
05_make_figures.py
------------------
Generate all publication-quality figures for the magnitude-direction dissociation paper.

Figures produced:
  fig1_magnitude_direction.png/.pdf
    Panel A: eGene rate by |R| decile (monotone bar chart with bootstrap CIs)
    Panel B: 2x2 factorial (|R| stratum x sign(R); eGene rates)
    Panel C: Direction accuracy (fraction of eGenes with R>0) by |R| decile

  fig2_finetune_direction.png/.pdf
    Panel A: FinetuneR distribution by base R sign (eGenes only)
    Panel B: PrediXcan R across 4 gene categories
    Panel C: StdObs vs |PearsonR| scatter (StdEnf proxy; labeled if available)

  fig3_model_comparison.png/.pdf
    Panel A: AIC comparison — 4-model horizontal bar chart
    Panel B: |R| vs PearsonR scatter colored by eGene status

Prerequisites:
  Run 03_decile_analysis.py and 04_predixcan_analysis.py first to generate:
    data/processed/decile_summary.csv
    data/processed/direction_by_decile.csv
    data/processed/aic_model_comparison.csv
    data/processed/gene_level_results.csv  (with abs_R, absR_decile, R_positive)
    data/processed/high_absR_Rpos_egenes.csv
    data/processed/high_absR_Rneg_egenes.csv
    data/processed/low_absR_Rpos_egenes.csv
    data/processed/low_absR_Rneg_egenes.csv

Outputs: manuscript/figures/
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
PROC_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
FIG_DIR = os.path.join(PROJECT_ROOT, "manuscript", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

GENE_FILE    = os.path.join(PROC_DIR, "gene_level_results.csv")
DECILE_FILE  = os.path.join(PROC_DIR, "decile_summary.csv")
DIR_FILE     = os.path.join(PROC_DIR, "direction_by_decile.csv")
AIC_FILE     = os.path.join(PROC_DIR, "aic_model_comparison.csv")
HI_POS_FILE  = os.path.join(PROC_DIR, "high_absR_Rpos_egenes.csv")
HI_NEG_FILE  = os.path.join(PROC_DIR, "high_absR_Rneg_egenes.csv")

# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
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

COLOR_HIGH  = "#1F77B4"   # blue — high |R|
COLOR_LOW   = "#D62728"   # red  — low |R|
COLOR_GRAY  = "#AAAAAA"
COLOR_RPOS  = "#2CA02C"   # green — R > 0 (correct direction)
COLOR_RNEG  = "#FF7F0E"   # orange — R < 0 (wrong direction)
COLOR_MEAN  = "#555555"


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

    if "abs_R" not in df.columns:
        df["abs_R"] = df["PearsonR"].abs()
    if "absR_decile" not in df.columns:
        df["absR_decile"] = assign_deciles(df["abs_R"])
    if "R_positive" not in df.columns:
        df["R_positive"] = (df["PearsonR"] > 0).astype(int)

    decile_df  = pd.read_csv(DECILE_FILE) if os.path.exists(DECILE_FILE) else None
    dir_df     = pd.read_csv(DIR_FILE)    if os.path.exists(DIR_FILE)    else None
    aic_df     = pd.read_csv(AIC_FILE)    if os.path.exists(AIC_FILE)    else None
    hi_pos_eg  = pd.read_csv(HI_POS_FILE) if os.path.exists(HI_POS_FILE) else pd.DataFrame()
    hi_neg_eg  = pd.read_csv(HI_NEG_FILE) if os.path.exists(HI_NEG_FILE) else pd.DataFrame()

    return df, decile_df, dir_df, aic_df, hi_pos_eg, hi_neg_eg


# ---------------------------------------------------------------------------
# FIG 1: Magnitude-direction dissociation
# Panel A: |R| decile bar chart (monotone eGene rate)
# Panel B: 2x2 factorial (|R| stratum x sign, eGene rate)
# Panel C: Direction accuracy (fraction R>0 among eGenes) by |R| decile
# ---------------------------------------------------------------------------
def make_fig1(df: pd.DataFrame, decile_df: pd.DataFrame, dir_df: pd.DataFrame) -> None:
    print("\nGenerating fig1_magnitude_direction ...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # ---- Panel A: |R| decile monotone bar chart ----
    ax = axes[0]
    deciles = decile_df["absR_decile"].values
    rates   = decile_df["egene_rate"].values * 100
    ci_lo   = decile_df["ci_lo"].values * 100
    ci_hi   = decile_df["ci_hi"].values * 100
    yerr_lo = rates - ci_lo
    yerr_hi = ci_hi - rates

    bar_colors = []
    for d in deciles:
        if d == 1:
            bar_colors.append(COLOR_LOW)
        elif d == 10:
            bar_colors.append(COLOR_HIGH)
        else:
            bar_colors.append(COLOR_GRAY)

    x = np.arange(len(deciles))
    ax.bar(x, rates, color=bar_colors, width=0.7, zorder=3,
           edgecolor="white", linewidth=0.5)
    ax.errorbar(x, rates, yerr=[yerr_lo, yerr_hi],
                fmt="none", color="black", capsize=4, linewidth=1.2, zorder=4)

    overall_rate = df["is_egene"].mean()
    ax.axhline(overall_rate * 100, color=COLOR_MEAN, linestyle="--",
               linewidth=1.5, zorder=2, label=f"Overall ({overall_rate*100:.1f}%)")

    ax.set_xticks(x)
    ax.set_xticklabels([f"D{d}" for d in deciles], fontsize=8)
    ax.set_xlabel("|Enformer Pearson R| Decile", labelpad=8)
    ax.set_ylabel("Brain Cortex eGene Rate (%)")
    ax.set_title("A   eGene rate by |R| decile", loc="left")
    ax.set_ylim(0, max(rates) * 1.25)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    # Annotate D1 and D10
    ax.text(0, rates[0] + yerr_hi[0] + 1.5, f"{rates[0]:.1f}%",
            ha="center", va="bottom", color=COLOR_LOW, fontsize=9, fontweight="bold")
    ax.text(len(deciles)-1, rates[-1] + yerr_hi[-1] + 1.5, f"{rates[-1]:.1f}%",
            ha="center", va="bottom", color=COLOR_HIGH, fontsize=9, fontweight="bold")

    patch_lo  = mpatches.Patch(color=COLOR_LOW,  label="D1 (lowest |R|)")
    patch_hi  = mpatches.Patch(color=COLOR_HIGH, label="D10 (highest |R|)")
    patch_mid = mpatches.Patch(color=COLOR_GRAY, label="D2-D9")
    line_ov   = plt.Line2D([0], [0], color=COLOR_MEAN, linestyle="--",
                            linewidth=1.5, label=f"Overall ({overall_rate*100:.1f}%)")
    ax.legend(handles=[patch_lo, patch_hi, patch_mid, line_ov],
              loc="upper left", frameon=False, fontsize=8)

    # ---- Panel B: 2x2 Factorial ----
    ax = axes[1]

    high_absR = df[df["absR_decile"] == 9]
    low_absR  = df[df["absR_decile"] == 0]

    groups = [
        ("High|R|\nR>0", high_absR[high_absR["PearsonR"] > 0], COLOR_RPOS, "//"),
        ("High|R|\nR<0", high_absR[high_absR["PearsonR"] < 0], COLOR_RNEG, "//"),
        ("Low|R|\nR>0",  low_absR[low_absR["PearsonR"] > 0],   COLOR_RPOS, ""),
        ("Low|R|\nR<0",  low_absR[low_absR["PearsonR"] < 0],   COLOR_RNEG, ""),
    ]

    x_pos = np.array([0, 1, 2.8, 3.8])
    fact_rates = []
    for i, (label, sub, color, hatch) in enumerate(groups):
        n_sub = len(sub)
        n_eg  = int(sub["is_egene"].sum())
        rate  = n_eg / n_sub * 100 if n_sub > 0 else 0
        fact_rates.append(rate)
        bar = ax.bar(x_pos[i], rate, width=0.75, color=color, zorder=3,
                     edgecolor="white", hatch=hatch, alpha=0.85)
        ax.text(x_pos[i], rate + 0.5, f"{rate:.1f}%",
                ha="center", va="bottom", fontsize=9, fontweight="bold")

    # Chi-squared annotation
    hi_pos_sub = high_absR[high_absR["PearsonR"] > 0]
    hi_neg_sub = high_absR[high_absR["PearsonR"] < 0]
    ct = np.array([
        [int((hi_pos_sub["is_egene"] == 1).sum()), int((hi_pos_sub["is_egene"] == 0).sum())],
        [int((hi_neg_sub["is_egene"] == 1).sum()), int((hi_neg_sub["is_egene"] == 0).sum())],
    ])
    _, chi2_p, _, _ = stats.chi2_contingency(ct)
    y_top = max(fact_rates[:2]) * 1.12
    ax.plot([x_pos[0], x_pos[1]], [y_top, y_top], "k-", linewidth=0.8)
    ax.text(x_pos[0:2].mean(), y_top * 1.01,
            f"p={chi2_p:.3f}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([g[0] for g in groups], fontsize=9)
    ax.set_ylabel("eGene Rate (%)")
    ax.set_title("B   2x2 factorial: |R| stratum x sign(R)", loc="left")
    ax.set_ylim(0, max(fact_rates) * 1.25)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    patch_rpos = mpatches.Patch(color=COLOR_RPOS, label="R>0")
    patch_rneg = mpatches.Patch(color=COLOR_RNEG, label="R<0")
    ax.legend(handles=[patch_rpos, patch_rneg], loc="upper right",
              frameon=False, fontsize=9)

    # ---- Panel C: Direction accuracy by |R| decile ----
    ax = axes[2]

    if dir_df is not None:
        dir_acc = dir_df["direction_accuracy"].values * 100
        dir_deciles = dir_df["absR_decile"].values

        bar_colors_c = []
        for d in dir_deciles:
            if d == 1:
                bar_colors_c.append(COLOR_LOW)
            elif d == 10:
                bar_colors_c.append(COLOR_HIGH)
            else:
                bar_colors_c.append(COLOR_GRAY)

        x_c = np.arange(len(dir_deciles))
        ax.bar(x_c, dir_acc, color=bar_colors_c, width=0.7, zorder=3,
               edgecolor="white", linewidth=0.5)
        ax.axhline(50, color="black", linestyle="--", linewidth=1.2, zorder=2,
                   label="Chance (50%)")

        ax.set_xticks(x_c)
        ax.set_xticklabels([f"D{d}" for d in dir_deciles], fontsize=8)
        ax.set_xlabel("|Enformer Pearson R| Decile", labelpad=8)
        ax.set_ylabel("eGenes with R>0 (%)")
        ax.set_title("C   Direction accuracy by |R| decile", loc="left")
        ax.set_ylim(0, 100)
        ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=9)

        # Annotate D10
        ax.text(len(dir_deciles)-1, dir_acc[-1] + 1.5, f"{dir_acc[-1]:.1f}%",
                ha="center", va="bottom", color=COLOR_HIGH, fontsize=9, fontweight="bold")
    else:
        ax.text(0.5, 0.5, "Run 03_decile_analysis.py\nto generate direction_by_decile.csv",
                ha="center", va="center", transform=ax.transAxes, fontsize=10,
                color="#888888")
        ax.set_title("C   Direction accuracy by |R| decile", loc="left")

    fig.tight_layout(pad=2.0)
    save_fig(fig, "fig1_magnitude_direction")
    plt.close(fig)


# ---------------------------------------------------------------------------
# FIG 2: Fine-tuning stability and PrediXcan comparison
# Panel A: FinetuneR distribution by base R sign (eGenes only)
# Panel B: PrediXcan R across 4 gene categories
# Panel C: StdObs vs |R| scatter (proxy for prediction spread analysis)
# ---------------------------------------------------------------------------
def make_fig2(df: pd.DataFrame) -> None:
    print("\nGenerating fig2_finetune_direction ...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    egenes = df[df["is_egene"] == 1]
    eg_rpos = egenes[egenes["PearsonR"] > 0]
    eg_rneg = egenes[egenes["PearsonR"] < 0]

    # ---- Panel A: FinetuneR distribution by base R sign ----
    ax = axes[0]

    if "PearsonRfineTuned" in df.columns:
        ft_rpos = eg_rpos["PearsonRfineTuned"].dropna().values
        ft_rneg = eg_rneg["PearsonRfineTuned"].dropna().values

        # Violin + box plot
        parts = ax.violinplot(
            [ft_rpos, ft_rneg],
            positions=[1, 2],
            showmedians=False,
            showextrema=False,
        )
        for pc, c in zip(parts["bodies"], [COLOR_RPOS, COLOR_RNEG]):
            pc.set_facecolor(c)
            pc.set_alpha(0.55)
            pc.set_edgecolor("white")

        bp = ax.boxplot(
            [ft_rpos, ft_rneg],
            positions=[1, 2],
            widths=0.12,
            patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(color="gray"),
            capprops=dict(color="gray"),
            flierprops=dict(marker=".", markersize=3, alpha=0.4, color="gray"),
        )
        for patch, c in zip(bp["boxes"], [COLOR_RPOS, COLOR_RNEG]):
            patch.set_facecolor("white")
            patch.set_edgecolor(c)
            patch.set_linewidth(1.5)

        ax.axhline(0, color="gray", linewidth=0.8, linestyle="-")

        # Mann-Whitney p
        if len(ft_rpos) > 1 and len(ft_rneg) > 1:
            _, pval_ft = stats.mannwhitneyu(ft_rpos, ft_rneg, alternative="two-sided")
            y_top = max(np.percentile(ft_rpos, 97), np.percentile(ft_rneg, 97))
            ax.plot([1, 2], [y_top * 1.05, y_top * 1.05], "k-", linewidth=0.8)
            sig = "***" if pval_ft < 0.001 else ("**" if pval_ft < 0.01 else "ns")
            ax.text(1.5, y_top * 1.07,
                    f"{sig}\np={pval_ft:.2e}", ha="center", va="bottom", fontsize=9)

        mean_rpos = np.mean(ft_rpos) if len(ft_rpos) > 0 else 0
        mean_rneg = np.mean(ft_rneg) if len(ft_rneg) > 0 else 0
        ax.text(1, ax.get_ylim()[0] * 0.95, f"mean={mean_rpos:+.3f}",
                ha="center", va="top", color=COLOR_RPOS, fontsize=8)
        ax.text(2, ax.get_ylim()[0] * 0.95, f"mean={mean_rneg:+.3f}",
                ha="center", va="top", color=COLOR_RNEG, fontsize=8)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["eGene\nR>0", "eGene\nR<0"])
        ax.set_ylabel("Fine-tuned Pearson R")
        ax.set_title("A   Fine-tuning stability by base R sign", loc="left")
    else:
        ax.text(0.5, 0.5, "PearsonRfineTuned\nnot available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888888")
        ax.set_title("A   Fine-tuning stability by base R sign", loc="left")

    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    # ---- Panel B: PrediXcan R across 4 categories ----
    ax = axes[1]

    if "PrediXcanR" in df.columns:
        non_eg = df[df["is_egene"] == 0]
        absR_median = df["abs_R"].median()

        cat1 = non_eg[non_eg["abs_R"] <= absR_median]["PrediXcanR"].dropna().values
        cat2 = non_eg[non_eg["abs_R"] > absR_median]["PrediXcanR"].dropna().values
        cat3 = df[(df["is_egene"] == 1) & (df["PearsonR"] < 0)]["PrediXcanR"].dropna().values
        cat4 = df[(df["is_egene"] == 1) & (df["PearsonR"] > 0)]["PrediXcanR"].dropna().values

        cats = [cat1, cat2, cat3, cat4]
        labels_cat = [
            "Non-eGene\nlow|R|",
            "Non-eGene\nhigh|R|",
            "eGene\nR<0",
            "eGene\nR>0",
        ]
        colors_cat = [COLOR_GRAY, "#7F7F7F", COLOR_RNEG, COLOR_RPOS]

        means  = [np.mean(c) if len(c) > 0 else 0 for c in cats]
        sems   = [stats.sem(c) if len(c) > 1 else 0 for c in cats]

        x_b = np.arange(len(labels_cat))
        bars = ax.bar(x_b, means, yerr=sems, color=colors_cat, width=0.6,
                      capsize=5, error_kw={"linewidth": 1.2},
                      edgecolor="white", zorder=3)

        for i, (m, label) in enumerate(zip(means, labels_cat)):
            ax.text(i, m + sems[i] + 0.005, f"{m:.3f}",
                    ha="center", va="bottom", fontsize=8.5)

        # Bracket: eGene R<0 significantly different from Non-eGene low|R|
        if len(cat3) > 1 and len(cat1) > 1:
            _, p_mw = stats.mannwhitneyu(cat3, cat1, alternative="two-sided")
            y_br = max(means[0] + sems[0], means[2] + sems[2]) + 0.02
            ax.plot([0, 2], [y_br, y_br], "k-", linewidth=0.8)
            sig = "***" if p_mw < 0.001 else ("**" if p_mw < 0.01 else "*" if p_mw < 0.05 else "ns")
            ax.text(1, y_br + 0.003, sig, ha="center", va="bottom", fontsize=10)

        ax.set_xticks(x_b)
        ax.set_xticklabels(labels_cat, fontsize=9)
        ax.set_ylabel("Mean PrediXcan R")
        ax.set_title("B   PrediXcan R across 4 gene categories", loc="left")
        ax.set_ylim(0, max(means) * 1.4)
        ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
        ax.set_axisbelow(True)
    else:
        ax.text(0.5, 0.5, "PrediXcanR\nnot available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888888")
        ax.set_title("B   PrediXcan R across 4 gene categories", loc="left")

    # ---- Panel C: StdObs vs |R| scatter ----
    ax = axes[2]

    if "StdObs" in df.columns:
        sample = df.sample(min(3000, len(df)), random_state=42)
        eg_mask = sample["is_egene"] == 1
        ax.scatter(sample.loc[~eg_mask, "abs_R"], sample.loc[~eg_mask, "StdObs"],
                   s=4, alpha=0.15, color=COLOR_GRAY, rasterized=True, label="Non-eGene")
        ax.scatter(sample.loc[eg_mask, "abs_R"], sample.loc[eg_mask, "StdObs"],
                   s=6, alpha=0.3, color=COLOR_HIGH, rasterized=True, label="eGene")

        # Spearman correlation
        spear_r, spear_p = stats.spearmanr(df["abs_R"], df["StdObs"])
        ax.text(0.05, 0.93, f"Spearman rho = {spear_r:.3f}\np = {spear_p:.2e}",
                transform=ax.transAxes, fontsize=9, color="#333333",
                verticalalignment="top")

        ax.set_xlabel("|Enformer Pearson R|")
        ax.set_ylabel("Observed Expression Std Dev (StdObs)")
        ax.set_title("C   Prediction magnitude vs. observed spread", loc="left")
        ax.legend(loc="upper right", frameon=False, fontsize=9, markerscale=3)
        ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
        ax.set_axisbelow(True)
    else:
        ax.text(0.5, 0.5, "StdObs\nnot available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888888")
        ax.set_title("C   Prediction magnitude vs. observed spread", loc="left")

    fig.tight_layout(pad=2.0)
    save_fig(fig, "fig2_finetune_direction")
    plt.close(fig)


# ---------------------------------------------------------------------------
# FIG 3: Formal model comparison
# Panel A: AIC comparison — 4-model horizontal bar chart
# Panel B: |R| vs PearsonR scatter colored by eGene status
# ---------------------------------------------------------------------------
def make_fig3(df: pd.DataFrame, aic_df: pd.DataFrame) -> None:
    print("\nGenerating fig3_model_comparison ...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ---- Panel A: AIC horizontal bar chart ----
    ax = axes[0]

    if aic_df is not None:
        model_labels = aic_df["model"].values
        aic_vals     = aic_df["AIC"].values
        delta_aic    = aic_df["delta_AIC"].values

        # Short labels for display
        short_labels = [
            "M1: R + R² (U-shape)",
            "M2: |R| (magnitude)",
            "M3: |R| + R",
            "M4: |R| + R + R²",
        ]

        # Color: best model highlighted
        best_idx = np.argmin(aic_vals)
        bar_colors = [COLOR_HIGH if i == best_idx else COLOR_GRAY
                      for i in range(len(aic_vals))]

        y = np.arange(len(short_labels))
        bars = ax.barh(y, aic_vals - aic_vals.min(),
                       color=bar_colors, height=0.55, zorder=3,
                       edgecolor="white", linewidth=0.5)

        # Annotate with actual AIC and delta
        for i, (val, daic) in enumerate(zip(aic_vals, delta_aic)):
            label_str = f"AIC={val:.0f}" + (f" (+{daic:.0f})" if daic > 0 else " (best)")
            ax.text(val - aic_vals.min() + 5, i, label_str,
                    va="center", fontsize=9)

        ax.set_yticks(y)
        ax.set_yticklabels(short_labels, fontsize=10)
        ax.set_xlabel("ΔAIC (relative to best model)")
        ax.set_title("A   Model comparison: AIC", loc="left")
        ax.set_xlim(0, (aic_vals.max() - aic_vals.min()) * 1.5)
        ax.xaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
        ax.set_axisbelow(True)

        # Annotate ΔAIC between M1 and M2
        m1_aic = aic_vals[0]
        m2_aic = aic_vals[1]
        ax.text(0.98, 0.05,
                f"ΔAIC(M2 vs M1) = {m1_aic - m2_aic:.0f}\n|R| model preferred",
                ha="right", va="bottom", transform=ax.transAxes,
                fontsize=9, color=COLOR_HIGH, style="italic")
    else:
        ax.text(0.5, 0.5, "Run 03_decile_analysis.py\nto generate aic_model_comparison.csv",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888888")
        ax.set_title("A   Model comparison: AIC", loc="left")

    # ---- Panel B: |R| vs PearsonR scatter colored by eGene status ----
    ax = axes[1]

    sample = df.sample(min(4000, len(df)), random_state=42)
    eg_mask = sample["is_egene"] == 1

    ax.scatter(sample.loc[~eg_mask, "PearsonR"], sample.loc[~eg_mask, "abs_R"],
               s=4, alpha=0.12, color=COLOR_GRAY, rasterized=True, label="Non-eGene")
    ax.scatter(sample.loc[eg_mask, "PearsonR"], sample.loc[eg_mask, "abs_R"],
               s=5, alpha=0.3, color=COLOR_HIGH, rasterized=True, label="eGene")

    ax.set_xlabel("Enformer Pearson R (direction)")
    ax.set_ylabel("|Enformer Pearson R| (magnitude)")
    ax.set_title("B   Magnitude vs. direction scatter", loc="left")
    ax.axvline(0, color="gray", linewidth=0.8, linestyle="--")
    ax.legend(loc="upper left", frameon=False, fontsize=9, markerscale=3)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.xaxis.grid(True, linestyle=":", alpha=0.4, zorder=1)
    ax.set_axisbelow(True)

    ax.text(0.55, 0.95,
            "eGenes cluster at high |R|\nregardless of R sign",
            transform=ax.transAxes, fontsize=9, color="#333333",
            va="top", style="italic")

    fig.tight_layout(pad=2.0)
    save_fig(fig, "fig3_model_comparison")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 60)
    print("05_make_figures.py  --  Generating publication figures")
    print("=" * 60)

    df, decile_df, dir_df, aic_df, hi_pos_eg, hi_neg_eg = load_data()
    print(f"Gene data: {len(df)} genes, {int(df['is_egene'].sum())} eGenes")

    if decile_df is not None:
        print(f"Decile summary: {len(decile_df)} rows")
    else:
        print("WARNING: decile_summary.csv not found -- run 03_decile_analysis.py first")

    if aic_df is not None:
        print(f"AIC comparison: {len(aic_df)} models")
    else:
        print("WARNING: aic_model_comparison.csv not found -- run 03_decile_analysis.py first")

    if dir_df is not None:
        print(f"Direction by decile: {len(dir_df)} rows")
    else:
        print("WARNING: direction_by_decile.csv not found -- run 03_decile_analysis.py first")

    make_fig1(df, decile_df, dir_df)
    make_fig2(df)
    make_fig3(df, aic_df)

    print(f"\nAll figures saved to: {FIG_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()
