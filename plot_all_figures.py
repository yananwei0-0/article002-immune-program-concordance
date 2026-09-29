#!/usr/bin/env python3
"""Rebuild all publication figures for article 002 from one workbook.

The script reads the locked result registries in
``data/Article2_Consolidated_Tables.xlsx`` and recreates four main figures and
six supplementary figures.  It does not recompute upstream gene scores or
statistical models; those registered numerical results are treated as inputs.
"""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parent
DEFAULT_WORKBOOK = ROOT / "data" / "Article2_Consolidated_Tables.xlsx"
DEFAULT_OUTDIR = ROOT / "generated"

PROGRAM_ORDER = [
    "MHC-I",
    "MHC-II",
    "IFN-gamma",
    "Cytolytic",
    "Checkpoint",
    "Myeloid",
    "TGF/EMT",
    "Proteasome",
]
CANCER_ORDER = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"]
CONTRAST_ORDER = ["malignant_minus_immune", "malignant_minus_stromal", "immune_minus_stromal"]
CONTRAST_TITLES = {
    "malignant_minus_immune": "Malignant − immune",
    "malignant_minus_stromal": "Malignant − stromal",
    "immune_minus_stromal": "Immune − stromal",
}
DISPLAY_PROGRAM = {"IFN-gamma": "IFN-γ"}

BLUE = "#0b78b4"
ORANGE = "#d95f02"
GREEN = "#179c77"
CHARCOAL = "#242a31"
MIDGREY = "#747d86"
LIGHTGREY = "#e4e7ea"
GRID = "#dde3e8"
WHITE = "#ffffff"

SEQUENTIAL = mpl.colors.LinearSegmentedColormap.from_list(
    "article002_sequential", ["#f3f5f6", "#8fc2d6", "#087c9d"]
)
SEQUENTIAL.set_bad(LIGHTGREY)
DIVERGING = mpl.colors.LinearSegmentedColormap.from_list(
    "article002_diverging", ["#b7352d", "#f7f7f7", "#2f73b3"]
)
DIVERGING.set_bad(LIGHTGREY)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.4,
            "axes.titlesize": 9.0,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.1,
            "ytick.labelsize": 7.1,
            "axes.linewidth": 0.8,
            "axes.edgecolor": CHARCOAL,
            "axes.labelcolor": CHARCOAL,
            "xtick.color": CHARCOAL,
            "ytick.color": CHARCOAL,
            "text.color": CHARCOAL,
            "svg.fonttype": "none",
            "savefig.facecolor": WHITE,
            "figure.facecolor": WHITE,
        }
    )


def display_program(label: str) -> str:
    return DISPLAY_PROGRAM.get(label, label)


def load_tables(workbook: Path) -> dict[str, pd.DataFrame]:
    required = [
        "S3_Primary",
        "S4a_FeatureSens",
        "S4b_Threshold1",
        "S5_ResidualContext",
        "S6a_TCGA_Full",
        "S6b_TCGA_Summary",
        "S6c_TCGA_Anchors",
        "S6d_TCGA_Subtypes",
        "S7_APOLLO",
        "S8b_scRNA_Source",
        "S9_Spatial_BRCA",
        "Derived_Heterogeneity",
        "Source_Alias",
    ]
    tables = pd.read_excel(workbook, sheet_name=required)
    validate_tables(tables)
    return tables


def validate_tables(t: dict[str, pd.DataFrame]) -> None:
    primary = t["S3_Primary"]
    evaluable = primary[primary["status"].eq("evaluable")]
    layer_counts = evaluable["comparison_layer"].value_counts().to_dict()
    assert len(primary) == 160, f"Expected 160 CPTAC rows, found {len(primary)}"
    assert len(evaluable) == 129, f"Expected 129 evaluable CPTAC rows, found {len(evaluable)}"
    assert layer_counts.get("proteomics") == 77
    assert layer_counts.get("phosphoproteomics") == 52
    assert len(t["S6a_TCGA_Full"]) == 5440
    assert int((t["S6a_TCGA_Full"]["q_value_bh"] < 0.05).sum()) == 4850
    anchors = set(
        zip(
            t["S6c_TCGA_Anchors"]["signature_label"],
            t["S6c_TCGA_Anchors"]["immune_signature_68"],
        )
    )
    assert len(anchors) == 17
    apollo_counts = t["S7_APOLLO"]["evidence_status"].value_counts().to_dict()
    assert apollo_counts == {"supported": 11, "not_evaluable": 7, "null_result": 5, "opposite": 1}
    scrna = t["S8b_scRNA_Source"]
    assert len(scrna) == 336
    assert int(scrna["status"].eq("evaluable").sum()) == 256
    assert int(scrna["within_source_support_window"].fillna(False).astype(bool).sum()) == 116
    spatial = t["S9_Spatial_BRCA"]
    assert len(spatial) == 24
    assert set(spatial["paired_sample_n"].astype(int)) == {4}


def panel_label(ax: mpl.axes.Axes, label: str, x: float = -0.12, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        clip_on=False,
    )


def clean_axis(ax: mpl.axes.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.7)
        ax.set_axisbelow(True)


def draw_heatmap(
    ax: mpl.axes.Axes,
    values: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    *,
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
    annotations: np.ndarray | None = None,
    annotate_threshold: float = 0.55,
    x_rotation: float = 0,
    x_fontsize: float = 7.0,
    y_fontsize: float = 7.0,
) -> mpl.image.AxesImage:
    masked = np.ma.masked_invalid(np.asarray(values, dtype=float))
    image = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(col_labels)), labels=col_labels)
    ax.set_yticks(np.arange(len(row_labels)), labels=row_labels)
    ax.tick_params(length=0, pad=3, labelsize=y_fontsize)
    ax.tick_params(axis="x", labelsize=x_fontsize)
    if x_rotation:
        for tick in ax.get_xticklabels():
            tick.set_rotation(x_rotation)
            tick.set_ha("right")
            tick.set_rotation_mode("anchor")
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color=WHITE, linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if annotations is not None:
        annotations = np.asarray(annotations, dtype=object)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                text = annotations[i, j]
                if text is None or text == "":
                    continue
                value = values[i, j]
                color = CHARCOAL
                if np.isfinite(value):
                    scaled = abs(float(norm(value)) - 0.5) * 2 if isinstance(norm, mpl.colors.TwoSlopeNorm) else float(norm(value))
                    if scaled >= annotate_threshold:
                        color = WHITE
                else:
                    color = MIDGREY
                ax.text(j, i, str(text), ha="center", va="center", fontsize=6.4, color=color)
    return image


def matrix_from_rows(
    df: pd.DataFrame,
    row_order: list[str],
    col_order: list[str],
    row_field: str,
    col_field: str,
    value_field: str,
) -> np.ndarray:
    return (
        df.pivot(index=row_field, columns=col_field, values=value_field)
        .reindex(index=row_order, columns=col_order)
        .to_numpy(dtype=float)
        .copy()
    )


def figure_1(t: dict[str, pd.DataFrame]) -> plt.Figure:
    primary = t["S3_Primary"].copy()
    primary["cancer_upper"] = primary["cancer"].str.upper()
    fig = plt.figure(figsize=(7.31, 6.47))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.033], hspace=0.14, wspace=0.05,
                          left=0.12, right=0.91, top=0.94, bottom=0.12)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    image = None
    for idx, (layer, title) in enumerate(
        [("proteomics", "RNA–protein concordance"), ("phosphoproteomics", "RNA–phosphoprotein concordance")]
    ):
        d = primary[primary["comparison_layer"].eq(layer)].copy()
        values = matrix_from_rows(d, PROGRAM_ORDER, CANCER_ORDER, "program_label", "cancer_upper", "spearman_rho")
        status = d.pivot(index="program_label", columns="cancer_upper", values="status").reindex(index=PROGRAM_ORDER, columns=CANCER_ORDER)
        qvals = d.pivot(index="program_label", columns="cancer_upper", values="q_value_bh").reindex(index=PROGRAM_ORDER, columns=CANCER_ORDER)
        annotations = np.empty(values.shape, dtype=object)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                if status.iloc[i, j] != "evaluable" or not np.isfinite(values[i, j]):
                    values[i, j] = np.nan
                    annotations[i, j] = "NE"
                else:
                    marker = "†" if float(qvals.iloc[i, j]) >= 0.05 else ""
                    annotations[i, j] = f"{values[i, j]:.2f}{marker}"
        ax = axes[idx]
        image = draw_heatmap(
            ax,
            values,
            [display_program(x) for x in PROGRAM_ORDER],
            CANCER_ORDER,
            cmap=SEQUENTIAL,
            norm=mpl.colors.Normalize(0, 1),
            annotations=annotations,
            annotate_threshold=0.58,
        )
        ax.set_title(title, loc="left", fontweight="bold", pad=7)
        panel_label(ax, chr(ord("a") + idx), x=-0.095, y=1.03)
        if idx == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("CPTAC cancer cohort", labelpad=7)
    cax = fig.add_subplot(gs[:, 1])
    cb = fig.colorbar(image, cax=cax, ticks=[0, 0.25, 0.5, 0.75, 1])
    cb.set_label("Spearman ρ", labelpad=7)
    cb.outline.set_linewidth(0.8)
    fig.text(0.12, 0.055, "Cell text is Spearman ρ. Grey cells were not evaluable (NE); † indicates q ≥ 0.05.", fontsize=6.8)
    return fig


def _layer_style(layer: str) -> tuple[str, str, float]:
    if layer == "proteomics":
        return BLUE, "o", -0.12
    return ORANGE, "s", 0.12


def figure_2(t: dict[str, pd.DataFrame]) -> plt.Figure:
    primary = t["S3_Primary"]
    hetero = t["Derived_Heterogeneity"]
    feature = t["S4a_FeatureSens"]
    fig = plt.figure(figsize=(7.38, 6.56))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0], width_ratios=[2.7, 1.0],
                          hspace=0.26, wspace=0.16, left=0.12, right=0.97, top=0.94, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)
    ax_c = fig.add_subplot(gs[1, :])

    ybase = np.arange(len(PROGRAM_ORDER))
    for layer in ["proteomics", "phosphoproteomics"]:
        color, marker, offset = _layer_style(layer)
        for yi, program in enumerate(PROGRAM_ORDER):
            d = primary[
                primary["comparison_layer"].eq(layer)
                & primary["program_label"].eq(program)
                & primary["status"].eq("evaluable")
            ]
            if d.empty:
                continue
            vals = d["spearman_rho"].astype(float).to_numpy()
            ypos = yi + offset
            q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
            ax_a.plot([vals.min(), vals.max()], [ypos, ypos], color=color, lw=0.8, alpha=0.45)
            ax_a.plot([q1, q3], [ypos, ypos], color=color, lw=3.0, solid_capstyle="round")
            supported = d["q_value_bh"].astype(float).to_numpy() < 0.05
            ax_a.scatter(vals[supported], np.full(supported.sum(), ypos), s=13, c=color, marker=marker, alpha=0.8, linewidths=0.25)
            if (~supported).any():
                ax_a.scatter(vals[~supported], np.full((~supported).sum(), ypos), s=18, facecolors=WHITE,
                             edgecolors=color, marker=marker, linewidths=0.9, zorder=4)
            ax_a.scatter([med], [ypos], s=30, c=color, marker="D", edgecolors=WHITE, linewidths=0.6, zorder=5)
    ax_a.set_yticks(ybase, labels=[display_program(x) for x in PROGRAM_ORDER])
    ax_a.invert_yaxis()
    ax_a.set_xlim(-0.03, 1.02)
    ax_a.set_xlabel("Cancer-specific Spearman ρ")
    ax_a.set_title("Cancer-specific estimates, median and IQR", loc="left", fontweight="bold")
    clean_axis(ax_a, "x")
    legend = [
        Line2D([0], [0], color=BLUE, marker="o", lw=2, label="RNA–protein"),
        Line2D([0], [0], color=ORANGE, marker="s", lw=2, label="RNA–phosphoprotein"),
        Line2D([0], [0], marker="o", markerfacecolor=WHITE, markeredgecolor=MIDGREY, lw=0, label="q ≥ 0.05"),
    ]
    ax_a.legend(handles=legend, loc="upper left", frameon=False, fontsize=6.8)
    panel_label(ax_a, "a")

    for layer in ["proteomics", "phosphoproteomics"]:
        color, marker, offset = _layer_style(layer)
        for yi, program in enumerate(PROGRAM_ORDER):
            d = hetero[hetero["comparison_layer"].eq(layer) & hetero["program_label"].eq(program)]
            if d.empty or int(d.iloc[0]["evaluable_cancer_n"]) == 0:
                if layer == "phosphoproteomics":
                    ax_b.text(0.04, yi + offset, "NE", color=MIDGREY, va="center", fontsize=6.5)
                continue
            ax_b.scatter(float(d.iloc[0]["i2_fraction"]), yi + offset, c=color, marker=marker, s=25, zorder=3)
    ax_b.set_xlim(0, 1.02)
    ax_b.set_xlabel("I² (descriptive)")
    ax_b.set_title("Across-cancer\nheterogeneity", loc="left", fontweight="bold", fontsize=8.0, linespacing=1.0)
    ax_b.tick_params(axis="y", left=False, labelleft=False)
    clean_axis(ax_b, "x")
    panel_label(ax_b, "b", x=-0.12)

    keys = ["cancer", "program_label", "comparison_layer"]
    p = primary[primary["status"].eq("evaluable")][keys + ["spearman_rho"]].rename(columns={"spearman_rho": "primary_rho"})
    f = feature[feature["status"].eq("evaluable")][keys + ["spearman_rho"]].rename(columns={"spearman_rho": "feature_rho"})
    joined = p.merge(f, on=keys, validate="one_to_one")
    for layer in ["proteomics", "phosphoproteomics"]:
        color, marker, _ = _layer_style(layer)
        d = joined[joined["comparison_layer"].eq(layer)]
        ax_c.scatter(d["primary_rho"], d["feature_rho"], c=color, marker=marker, s=15, alpha=0.78, linewidths=0.2)
    ax_c.plot([-0.05, 1], [-0.05, 1], ls="--", color=MIDGREY, lw=0.9)
    joined["abs_delta"] = (joined["feature_rho"] - joined["primary_rho"]).abs()
    row = joined.loc[joined["abs_delta"].idxmax()]
    ax_c.annotate(
        f"Largest |Δρ|={row['abs_delta']:.3f}\n{str(row['cancer']).upper()}, {display_program(str(row['program_label']))}",
        xy=(row["primary_rho"], row["feature_rho"]),
        xytext=(0.03, 0.88),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "-", "color": MIDGREY, "lw": 0.8},
        fontsize=6.6,
        va="top",
    )
    phospho = joined[joined["comparison_layer"].eq("phosphoproteomics")]
    ax_c.text(
        0.98,
        0.04,
        f"{len(joined)} jointly evaluable rows\nSign concordance: {(np.sign(joined.primary_rho)==np.sign(joined.feature_rho)).sum()}/{len(joined)}\n"
        f"Phosphoprotein median |Δρ|={np.median(np.abs(phospho.feature_rho-phospho.primary_rho)):.3f}",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color=MIDGREY,
        bbox={"boxstyle": "square,pad=0.25", "facecolor": WHITE, "edgecolor": GRID, "linewidth": 0.7},
    )
    ax_c.set_xlim(-0.05, 1.0)
    ax_c.set_ylim(-0.05, 1.0)
    ax_c.set_xlabel("Primary gene-equal Spearman ρ")
    ax_c.set_ylabel("Feature-row sensitivity Spearman ρ")
    ax_c.set_title("Scoring sensitivity", loc="left", fontweight="bold")
    clean_axis(ax_c, "both")
    panel_label(ax_c, "c", x=-0.085)
    return fig


def figure_3(t: dict[str, pd.DataFrame]) -> plt.Figure:
    summary = t["S6b_TCGA_Summary"].set_index("signature_label")
    apollo = t["S7_APOLLO"]
    fig = plt.figure(figsize=(7.39, 3.94))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.5, 0.045], wspace=0.30,
                          left=0.10, right=0.93, top=0.86, bottom=0.23)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    vals = summary.reindex(PROGRAM_ORDER)["selected_anchor_median_rho"].astype(float).to_numpy()
    y = np.arange(len(PROGRAM_ORDER))
    for i, value in enumerate(vals):
        color = MIDGREY if PROGRAM_ORDER[i] == "Proteasome" else GREEN
        ax_a.plot([0, value], [i, i], color=color, lw=1.3, alpha=0.65)
        ax_a.scatter(value, i, s=30, c=color, marker="s" if PROGRAM_ORDER[i] == "Proteasome" else "o", zorder=3)
        ax_a.text(value - 0.02 if value > 0.5 else value + 0.025, i, f"{value:.3f}",
                  ha="right" if value > 0.5 else "left", va="center", fontsize=7.0)
    ax_a.set_yticks(y, labels=[display_program(x) for x in PROGRAM_ORDER])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 1.0)
    ax_a.set_xlabel("Median within-cancer Spearman ρ")
    ax_a.set_title("TCGA RNA locked pairs", loc="left", fontweight="bold")
    clean_axis(ax_a, "x")
    panel_label(ax_a, "a", x=-0.23)

    columns = [
        ("APOLLO-LUAD", "proteomics", "LUAD\nProtein"),
        ("APOLLO-LUAD", "phosphoproteomics", "LUAD\nPhosphoprotein"),
        ("APOLLO-OV", "proteomics", "OV\nProtein"),
    ]
    values = np.full((len(PROGRAM_ORDER), len(columns)), np.nan)
    annotations = np.empty(values.shape, dtype=object)
    for i, program in enumerate(PROGRAM_ORDER):
        for j, (cohort, layer, _) in enumerate(columns):
            d = apollo[
                apollo["cohort_label"].eq(cohort)
                & apollo["comparison_layer"].eq(layer)
                & apollo["program_label"].eq(program)
            ]
            if d.empty or d.iloc[0]["status"] != "evaluable":
                annotations[i, j] = "NE"
                continue
            row = d.iloc[0]
            value = float(row["spearman_rho"])
            values[i, j] = value
            marker = "*" if row["evidence_status"] == "supported" else "‡" if row["evidence_status"] == "opposite" else ""
            annotations[i, j] = f"{value:.2f}{marker}"
    image = draw_heatmap(
        ax_b,
        values,
        [display_program(x) for x in PROGRAM_ORDER],
        [x[2] for x in columns],
        cmap=DIVERGING,
        norm=mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
        annotations=annotations,
        annotate_threshold=0.52,
    )
    ax_b.set_xlabel("Cohort and omic layer", labelpad=8)
    ax_b.set_title("APOLLO cross-omic concordance", loc="left", fontweight="bold")
    panel_label(ax_b, "b", x=-0.13)
    cb = fig.colorbar(image, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label("Spearman ρ", labelpad=7)
    fig.text(
        0.50,
        0.08,
        "* FDR-supported positive; ‡ FDR-supported opposite direction; NE, not evaluable.",
        ha="center",
        fontsize=6.7,
    )
    return fig


def figure_4(t: dict[str, pd.DataFrame]) -> plt.Figure:
    scrna = t["S8b_scRNA_Source"]
    spatial = t["S9_Spatial_BRCA"]
    fig = plt.figure(figsize=(7.64, 6.07))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 0.90], hspace=0.45,
                             left=0.10, right=0.92, top=0.93, bottom=0.10)
    top = outer[0].subgridspec(1, 3, wspace=0.28)
    top_axes = [fig.add_subplot(top[0, i]) for i in range(3)]
    y = np.arange(len(PROGRAM_ORDER))
    for j, contrast in enumerate(CONTRAST_ORDER):
        ax = top_axes[j]
        neg_counts, pos_counts = [], []
        for program in PROGRAM_ORDER:
            d = scrna[
                scrna["program_label"].eq(program)
                & scrna["contrast"].eq(contrast)
                & scrna["within_source_support_window"].fillna(False).astype(bool)
            ]
            neg_counts.append(int((d["median_delta"] < 0).sum()))
            pos_counts.append(int((d["median_delta"] > 0).sum()))
        neg = np.asarray(neg_counts)
        pos = np.asarray(pos_counts)
        ax.barh(y, -neg, color=ORANGE, height=0.55, label="Negative median difference")
        ax.barh(y, pos, color=BLUE, height=0.55, label="Positive median difference")
        for yi, n in enumerate(neg):
            if n:
                ax.text(-n - 0.18, yi, str(n), ha="right", va="center", fontsize=6.6)
        for yi, n in enumerate(pos):
            if n:
                ax.text(n + 0.18, yi, str(n), ha="left", va="center", fontsize=6.6)
        ax.axvline(0, color=CHARCOAL, lw=0.8)
        ax.set_xlim(-8.4, 8.4)
        ax.set_xticks([-6, -3, 0, 3, 6], labels=["6", "3", "0", "3", "6"])
        ax.set_yticks(y, labels=[display_program(x) for x in PROGRAM_ORDER] if j == 0 else [])
        ax.invert_yaxis()
        ax.set_title(CONTRAST_TITLES[contrast], fontweight="bold", pad=6)
        clean_axis(ax, "x")
    panel_label(top_axes[0], "a", x=-0.20)
    fig.text(0.51, 0.51, "FDR-supported source–cancer rows (count)", ha="center", fontsize=8.0)
    handles = [Rectangle((0, 0), 1, 1, color=ORANGE), Rectangle((0, 0), 1, 1, color=BLUE)]
    fig.legend(handles, ["Negative median difference", "Positive median difference"], loc="upper center",
               bbox_to_anchor=(0.60, 0.465), ncol=2, frameon=False, fontsize=6.8)

    ax_b = fig.add_subplot(outer[1])
    # Reserve a clean band between panels for the shared legend and title.
    ax_b.set_position([0.10, 0.14, 0.82, 0.25])
    values = matrix_from_rows(spatial, PROGRAM_ORDER, CONTRAST_ORDER, "signature_label", "contrast", "median_delta")
    annotations = np.vectorize(lambda x: f"{x:.2f}")(values)
    image = draw_heatmap(
        ax_b,
        values,
        [display_program(x) for x in PROGRAM_ORDER],
        [CONTRAST_TITLES[x] for x in CONTRAST_ORDER],
        cmap=DIVERGING,
        norm=mpl.colors.TwoSlopeNorm(vmin=-0.6, vcenter=0, vmax=0.6),
        annotations=annotations,
        annotate_threshold=0.62,
    )
    ax_b.set_xlabel("Compartment contrast", labelpad=6)
    ax_b.set_title("BRCA Visium sample-level median contrasts (paired n=4)", loc="left", fontweight="bold")
    panel_label(ax_b, "b", x=-0.095)
    cax = fig.add_axes([0.94, 0.14, 0.018, 0.25])
    cb = fig.colorbar(image, cax=cax, ticks=[-0.6, -0.3, 0, 0.3, 0.6])
    cb.set_label("Median score difference", labelpad=7)
    fig.text(0.10, 0.025, "Single-cell counts are source-bounded summaries, not pooled independent tests. Spatial contrasts are descriptive; no spot-level P values were used.", fontsize=6.2)
    return fig


def figure_s1(t: dict[str, pd.DataFrame]) -> plt.Figure:
    primary = t["S3_Primary"].copy()
    primary["cancer_upper"] = primary["cancer"].str.upper()
    fig = plt.figure(figsize=(7.31, 6.32))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.035], hspace=0.20, wspace=0.05,
                          left=0.12, right=0.91, top=0.94, bottom=0.12)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    image = None
    for idx, (layer, title) in enumerate(
        [("proteomics", "RNA–protein complete pairs"), ("phosphoproteomics", "RNA–phosphoprotein complete pairs")]
    ):
        d = primary[primary["comparison_layer"].eq(layer)]
        values = matrix_from_rows(d, PROGRAM_ORDER, CANCER_ORDER, "program_label", "cancer_upper", "complete_pair_n")
        status = d.pivot(index="program_label", columns="cancer_upper", values="status").reindex(index=PROGRAM_ORDER, columns=CANCER_ORDER)
        annotations = np.vectorize(lambda x: str(int(x)))(values)
        ax = axes[idx]
        image = draw_heatmap(
            ax,
            values,
            [display_program(x) for x in PROGRAM_ORDER],
            CANCER_ORDER,
            cmap=SEQUENTIAL,
            norm=mpl.colors.Normalize(0, 125),
            annotations=annotations,
            annotate_threshold=0.58,
        )
        for i in range(len(PROGRAM_ORDER)):
            for j in range(len(CANCER_ORDER)):
                if status.iloc[i, j] != "evaluable":
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=ORANGE, linewidth=1.1))
        ax.set_title(title, loc="left", fontweight="bold", pad=7)
        panel_label(ax, chr(ord("a") + idx), x=-0.095)
        if idx == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("CPTAC cancer cohort", labelpad=7)
    cax = fig.add_subplot(gs[:, 1])
    cb = fig.colorbar(image, cax=cax, ticks=[0, 20, 40, 60, 80, 100, 120])
    cb.set_label("Complete matched pairs (n)", labelpad=7)
    fig.text(0.12, 0.05, "Outlined cells did not meet the primary n≥20/distinct-value gate.", fontsize=6.8)
    return fig


def figure_s2(t: dict[str, pd.DataFrame]) -> plt.Figure:
    residual = t["S5_ResidualContext"].copy()
    supported = residual[
        residual["linear_status"].eq("evaluable")
        & residual["linear_q_value_bh"].notna()
        & residual["linear_q_value_bh"].lt(0.05)
    ].copy()
    context_order = [
        "estimate_tumor_purity",
        "estimate_immune_score",
        "cibersort_t_cell",
        "cibersort_myeloid",
    ]
    context_label = {
        "estimate_tumor_purity": "Purity",
        "estimate_immune_score": "ESTIMATE immune",
        "cibersort_t_cell": "CIBERSORT T-cell",
        "cibersort_myeloid": "CIBERSORT myeloid",
    }
    categories = []
    for layer, prefix in [("proteomics", "P"), ("phosphoproteomics", "pP")]:
        for context in context_order:
            categories.append((layer, context, f"{prefix}: {context_label[context]}"))

    fig = plt.figure(figsize=(7.31, 4.42))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.035], wspace=0.04,
                          left=0.15, right=0.91, top=0.86, bottom=0.25)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    max_count = 1
    for yi, program in enumerate(PROGRAM_ORDER):
        for xi, (layer, context, _) in enumerate(categories):
            d = supported[
                supported["program_label"].eq(program)
                & supported["comparison_layer"].eq(layer)
                & supported["context_variable"].eq(context)
            ]
            count = len(d)
            if count == 0:
                continue
            max_count = max(max_count, count)
            med = float(d["linear_rho"].median())
            ax.scatter(xi, yi, s=32 + count * 28, c=[med], cmap=DIVERGING,
                       norm=mpl.colors.TwoSlopeNorm(vmin=-0.6, vcenter=0, vmax=0.6),
                       edgecolors=CHARCOAL, linewidths=0.35, zorder=3)
            ax.text(xi, yi, str(count), ha="center", va="center", fontsize=6.2,
                    color=WHITE if abs(med) > 0.28 else CHARCOAL, zorder=4)
    ax.set_xlim(-0.6, len(categories) - 0.4)
    ax.set_ylim(len(PROGRAM_ORDER) - 0.5, -0.5)
    ax.set_xticks(np.arange(len(categories)), labels=[x[2] for x in categories], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(PROGRAM_ORDER)), labels=[display_program(x) for x in PROGRAM_ORDER])
    ax.set_xlabel("Layer and residual-context variable", labelpad=8)
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.TwoSlopeNorm(vmin=-0.6, vcenter=0, vmax=0.6), cmap=DIVERGING)
    cb = fig.colorbar(sm, cax=cax, ticks=[-0.4, -0.2, 0, 0.2, 0.4])
    cb.set_label("Median residual-context ρ", labelpad=7)
    fig.text(0.15, 0.91, "Circle area and numeral show the number of linear-residual rows with q<0.05.", fontsize=7.0)
    return fig


def figure_s3(t: dict[str, pd.DataFrame]) -> plt.Figure:
    full = t["S6a_TCGA_Full"]
    anchors = t["S6c_TCGA_Anchors"]
    metric_order = list(dict.fromkeys(full["correlated_program_or_gene"].tolist()))
    assert len(metric_order) == 68
    median = (
        full.groupby(["program_label", "correlated_program_or_gene"], sort=False)["spearman_rho"]
        .median()
        .unstack("correlated_program_or_gene")
        .reindex(index=PROGRAM_ORDER, columns=metric_order)
    )
    anchor_pairs = set(zip(anchors["signature_label"], anchors["immune_signature_68"]))
    fig = plt.figure(figsize=(10.0, 8.2))
    panel_positions = [[0.200, 0.145, 0.220, 0.775], [0.670, 0.145, 0.220, 0.775]]
    image = None
    for panel_index, position in enumerate(panel_positions):
        start = panel_index * 34
        stop = start + 34
        metrics = metric_order[start:stop]
        values = median.loc[:, metrics].T.to_numpy(dtype=float)
        ax = fig.add_axes(position)
        image = draw_heatmap(
            ax,
            values,
            metrics,
            [display_program(x) for x in PROGRAM_ORDER],
            cmap=DIVERGING,
            norm=mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
            annotations=None,
            x_rotation=34,
            x_fontsize=8.0,
            y_fontsize=7.2,
        )
        for row_index, metric in enumerate(metrics):
            for col_index, program in enumerate(PROGRAM_ORDER):
                if (program, metric) in anchor_pairs:
                    ax.add_patch(Rectangle((col_index - 0.5, row_index - 0.5), 1, 1, fill=False,
                                           edgecolor=CHARCOAL, linewidth=1.15))
        ax.set_title(f"Immune metrics {start + 1}–{stop}", loc="left", pad=2, fontweight="bold")
        panel_label(ax, chr(ord("a") + panel_index), x=-0.77, y=1.025)
    cax = fig.add_axes([0.925, 0.31, 0.016, 0.42])
    cb = fig.colorbar(image, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label("Median within-cancer Spearman ρ", labelpad=8)
    return fig


def figure_s4(t: dict[str, pd.DataFrame]) -> plt.Figure:
    subtypes = t["S6d_TCGA_Subtypes"]
    subtype_order = [f"Immune C{i}" for i in range(1, 7)]
    values = matrix_from_rows(subtypes, PROGRAM_ORDER, subtype_order, "signature_label", "immune_subtype_short", "median_score")
    annotations = np.vectorize(lambda x: f"{x:.2f}")(values)
    ns = (
        subtypes[["immune_subtype_short", "n"]]
        .drop_duplicates("immune_subtype_short")
        .set_index("immune_subtype_short")
        .reindex(subtype_order)["n"]
        .astype(int)
    )
    col_labels = [f"{s}\n(n={ns[s]})" + ("†" if s == "Immune C5" else "") for s in subtype_order]
    vmax = max(abs(float(np.nanmin(values))), abs(float(np.nanmax(values))))
    fig = plt.figure(figsize=(7.31, 4.37))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.035], wspace=0.05,
                          left=0.15, right=0.91, top=0.87, bottom=0.21)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    image = draw_heatmap(
        ax,
        values,
        [display_program(x) for x in PROGRAM_ORDER],
        col_labels,
        cmap=DIVERGING,
        norm=mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
        annotations=annotations,
        annotate_threshold=0.60,
    )
    ax.set_xlabel("TCGA pan-cancer immune subtype", labelpad=8)
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("Median program score", labelpad=7)
    fig.text(0.15, 0.08, "Values are median standardized RNA program scores. † Immune C5 (n=4) is descriptive only.", fontsize=6.8)
    return fig


def figure_s5(t: dict[str, pd.DataFrame]) -> plt.Figure:
    scrna = t["S8b_scRNA_Source"].copy()
    aliases = dict(zip(t["Source_Alias"]["source_id"], t["Source_Alias"]["source_alias"]))
    scrna["source_alias"] = scrna["source_id"].map(aliases)
    scrna["row_key"] = scrna["cancer_code"].str.upper() + " / " + scrna["source_alias"]
    cancer_rank = {c: i for i, c in enumerate(CANCER_ORDER)}
    rows = (
        scrna[["cancer_code", "source_alias", "row_key"]]
        .drop_duplicates()
        .assign(cancer_rank=lambda x: x["cancer_code"].str.upper().map(cancer_rank))
        .sort_values(["cancer_rank", "source_alias"])["row_key"]
        .tolist()
    )
    vmax = max(2.5, float(scrna["median_delta"].abs().max(skipna=True)))
    fig = plt.figure(figsize=(10.31, 5.92))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.04], wspace=0.12,
                          left=0.09, right=0.93, top=0.86, bottom=0.18)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    image = None
    for j, contrast in enumerate(CONTRAST_ORDER):
        d = scrna[scrna["contrast"].eq(contrast)]
        values = matrix_from_rows(d, rows, PROGRAM_ORDER, "row_key", "program_label", "median_delta")
        status = d.pivot(index="row_key", columns="program_label", values="status").reindex(index=rows, columns=PROGRAM_ORDER)
        support = d.pivot(index="row_key", columns="program_label", values="within_source_support_window").reindex(index=rows, columns=PROGRAM_ORDER)
        values[status.to_numpy() != "evaluable"] = np.nan
        ax = axes[j]
        image = draw_heatmap(
            ax,
            values,
            rows if j == 0 else [""] * len(rows),
            [display_program(x) for x in PROGRAM_ORDER],
            cmap=DIVERGING,
            norm=mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
            annotations=None,
            x_rotation=35,
            x_fontsize=6.4,
            y_fontsize=6.4,
        )
        for i in range(len(rows)):
            for k in range(len(PROGRAM_ORDER)):
                if bool(support.iloc[i, k]) and np.isfinite(values[i, k]):
                    ax.scatter(k, i, s=7, c=CHARCOAL, marker="o", linewidths=0)
        ax.set_title(CONTRAST_TITLES[contrast], fontweight="bold", pad=6)
        ax.set_xlabel("Program", labelpad=8)
        if j == 0:
            ax.set_ylabel("Cancer / source alias")
    cax = fig.add_subplot(gs[0, 3])
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("Paired-donor median difference", labelpad=7)
    fig.text(0.09, 0.06, "Black dots mark source–cancer rows meeting paired n≥6 and within-family q<0.05.", fontsize=6.8)
    return fig


def figure_s6(t: dict[str, pd.DataFrame]) -> plt.Figure:
    primary = t["S3_Primary"]
    threshold = t["S4b_Threshold1"]
    keys = ["cancer", "program_label", "comparison_layer"]
    p = primary[primary["status"].eq("evaluable")][keys + ["spearman_rho", "max_abs_leave_one_out_rho_delta"]].rename(columns={"spearman_rho": "primary_rho"})
    q = threshold[threshold["status"].eq("evaluable")][keys + ["spearman_rho"]].rename(columns={"spearman_rho": "threshold_rho"})
    joined = p.merge(q, on=keys, validate="one_to_one")
    fig = plt.figure(figsize=(7.32, 3.57))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.28,
                          left=0.10, right=0.97, top=0.86, bottom=0.19)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    for layer in ["proteomics", "phosphoproteomics"]:
        color, marker, _ = _layer_style(layer)
        d = joined[joined["comparison_layer"].eq(layer)]
        ax_a.scatter(d["primary_rho"], d["threshold_rho"], c=color, marker=marker, s=15, alpha=0.82, linewidths=0.2)
    ax_a.plot([-0.05, 1], [-0.05, 1], ls="--", color=MIDGREY, lw=0.8)
    ax_a.set_xlim(-0.05, 1.0)
    ax_a.set_ylim(-0.05, 1.0)
    ax_a.set_xlabel("Primary 3-gene threshold ρ")
    ax_a.set_ylabel("1-gene threshold ρ")
    ax_a.text(0.04, 0.94, f"Jointly evaluable: {len(joined)}\nThreshold-1 evaluable: {int(threshold.status.eq('evaluable').sum())}/160",
              transform=ax_a.transAxes, va="top", fontsize=6.8)
    clean_axis(ax_a, "both")
    panel_label(ax_a, "a", x=-0.20)

    ybase = np.arange(len(PROGRAM_ORDER))
    rng = np.random.default_rng(20260929)
    for layer in ["proteomics", "phosphoproteomics"]:
        color, marker, offset = _layer_style(layer)
        for yi, program in enumerate(PROGRAM_ORDER):
            d = p[p["comparison_layer"].eq(layer) & p["program_label"].eq(program)]
            if d.empty:
                continue
            jitter = rng.uniform(-0.045, 0.045, len(d))
            ax_b.scatter(d["max_abs_leave_one_out_rho_delta"], yi + offset + jitter,
                         c=color, marker=marker, s=14, alpha=0.82, linewidths=0.2)
    maximum = float(p["max_abs_leave_one_out_rho_delta"].max())
    ax_b.axvline(maximum, color=MIDGREY, ls="--", lw=0.8)
    ax_b.text(maximum - 0.002, len(PROGRAM_ORDER) - 0.32, f"max={maximum:.3f}", ha="right", va="bottom", fontsize=6.4)
    ax_b.set_yticks(ybase, labels=[display_program(x) for x in PROGRAM_ORDER])
    ax_b.invert_yaxis()
    ax_b.set_xlim(-0.003, max(0.125, maximum * 1.04))
    ax_b.set_xlabel("Maximum |leave-one-case-out Δρ|")
    clean_axis(ax_b, "x")
    panel_label(ax_b, "b", x=-0.14)
    return fig


FIGURES = {
    "Fig1_CPTAC_primary_concordance": figure_1,
    "Fig2_heterogeneity_and_scoring_robustness": figure_2,
    "Fig3_TCGA_and_APOLLO_triangulation": figure_3,
    "Fig4_single_cell_and_spatial_compartments": figure_4,
    "FigS1_complete_pair_counts": figure_s1,
    "FigS2_residual_context_support": figure_s2,
    "FigS3_TCGA_full_correlation_atlas": figure_s3,
    "FigS4_TCGA_immune_subtypes": figure_s4,
    "FigS5_source_resolved_single_cell": figure_s5,
    "FigS6_threshold_and_influence_sensitivities": figure_s6,
}


def save_figure(fig: plt.Figure, stem: str, outdir: Path, formats: Iterable[str]) -> list[Path]:
    saved: list[Path] = []
    formats = set(formats)
    if "svg" in formats:
        path = outdir / "SVG_editable" / f"{stem}.svg"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, format="svg", facecolor=WHITE)
        saved.append(path)
    if "png" in formats:
        path = outdir / "PNG_300dpi" / f"{stem}.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="article002_png_") as tmpdir:
            tmp_path = Path(tmpdir) / path.name
            fig.savefig(tmp_path, format="png", dpi=300, facecolor=WHITE)
            with Image.open(tmp_path) as check:
                check.load()
                if check.info.get("dpi") and round(check.info["dpi"][0]) != 300:
                    raise RuntimeError(f"Unexpected PNG resolution for {stem}: {check.info.get('dpi')}")
            shutil.copyfile(tmp_path, path)
        saved.append(path)
    if "tiff" in formats:
        path = outdir / "TIFF_600dpi" / f"{stem}.tiff"
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="article002_tiff_") as tmpdir:
            raw_path = Path(tmpdir) / f"{stem}.raw.tiff"
            rgb_path = Path(tmpdir) / f"{stem}.tiff"
            fig.savefig(raw_path, format="tiff", dpi=600, facecolor=WHITE, pil_kwargs={"compression": "tiff_lzw"})
            with Image.open(raw_path) as source:
                source.load()
                rgb = source.convert("RGB").copy()
            rgb.save(rgb_path, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
            with Image.open(rgb_path) as check:
                check.load()
                if check.mode != "RGB" or round(check.info.get("dpi", (0, 0))[0]) != 600:
                    raise RuntimeError(
                        f"Invalid TIFF export for {stem}: mode={check.mode}, dpi={check.info.get('dpi')}"
                    )
            shutil.copyfile(rgb_path, path)
        saved.append(path)
    plt.close(fig)
    return saved


def write_manifest(workbook: Path, outdir: Path, figures: list[str], files: list[Path], t: dict[str, pd.DataFrame]) -> None:
    def portable_path(path: Path) -> str:
        resolved = path.resolve()
        try:
            return str(resolved.relative_to(ROOT))
        except ValueError:
            return str(resolved)

    manifest = {
        "workbook": portable_path(workbook),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "package_versions": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
            "Pillow": Image.__version__,
        },
        "validated_counts": {
            "cptac_registered": len(t["S3_Primary"]),
            "cptac_evaluable": int(t["S3_Primary"]["status"].eq("evaluable").sum()),
            "tcga_registry": len(t["S6a_TCGA_Full"]),
            "tcga_q_lt_0_05": int(t["S6a_TCGA_Full"]["q_value_bh"].lt(0.05).sum()),
            "tcga_fixed_anchors": len(t["S6c_TCGA_Anchors"]),
            "single_cell_registry": len(t["S8b_scRNA_Source"]),
            "spatial_contrasts": len(t["S9_Spatial_BRCA"]),
        },
        "figures": figures,
        "files": [portable_path(p) for p in files],
    }
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


def validate_saved_files(files: list[Path]) -> None:
    """Fail the run if any declared deliverable is absent, truncated, or mis-specified."""
    for path in files:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
        if path.suffix.lower() == ".png":
            with Image.open(path) as image:
                image.load()
                if round(image.info.get("dpi", (0, 0))[0]) != 300:
                    raise RuntimeError(f"Unexpected PNG resolution: {path}")
        elif path.suffix.lower() in {".tif", ".tiff"}:
            with Image.open(path) as image:
                image.load()
                if image.mode != "RGB" or round(image.info.get("dpi", (0, 0))[0]) != 600:
                    raise RuntimeError(f"Unexpected TIFF properties: {path}")


def cleanup_temporary_files(outdir: Path) -> None:
    for pattern in ("*.tmp.*", "*.raw.*", "*.rgb.*"):
        for path in outdir.rglob(pattern):
            path.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK, help="Consolidated result workbook")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Output directory")
    parser.add_argument("--only", nargs="*", choices=list(FIGURES), default=None, help="Optional subset of figure stems")
    parser.add_argument("--formats", nargs="+", choices=["png", "svg", "tiff"], default=["png", "svg", "tiff"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()
    tables = load_tables(args.workbook)
    selected = args.only or list(FIGURES)
    saved: list[Path] = []
    for stem in selected:
        print(f"Rendering {stem} ...", flush=True)
        fig = FIGURES[stem](tables)
        saved.extend(save_figure(fig, stem, args.outdir, args.formats))
    cleanup_temporary_files(args.outdir)
    validate_saved_files(saved)
    write_manifest(args.workbook, args.outdir, selected, saved, tables)
    print(f"Generated {len(selected)} figures and {len(saved)} files in {args.outdir}")


if __name__ == "__main__":
    main()
