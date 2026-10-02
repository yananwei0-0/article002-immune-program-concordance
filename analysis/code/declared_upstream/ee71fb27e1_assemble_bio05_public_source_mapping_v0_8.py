#!/usr/bin/env python3
from __future__ import annotations

import os as _workspace_os


def _workspace_path(value: str) -> str:
    """Resolve portable workspace markers without embedding a local user path."""
    mapping = {
        "${DATA_WORKSPACE}": "DATA_WORKSPACE",
        "${RESEARCH_WORKSPACE}": "RESEARCH_WORKSPACE",
        "${HOME}": "HOME",
    }
    resolved = value
    for marker, variable in mapping.items():
        if marker not in resolved:
            continue
        root = _workspace_os.environ.get(variable, "").rstrip("/")
        if not root:
            raise RuntimeError(f"Set {variable} before running this script")
        resolved = resolved.replace(marker, root)
    return resolved


import importlib.util
import platform
import zlib
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


PROJECT_ROOT = Path(
    _workspace_path('${DATA_WORKSPACE}/14_pan_cancer_proteogenomic_immune_evasion')
)
VERSION = "v0_8"


def load_v07_module():
    script_path = PROJECT_ROOT / "scripts/run_bio05_scrna_source_mapping_v0_7.py"
    spec = importlib.util.spec_from_file_location("bio05_scrna_v07", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load base script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V07 = load_v07_module()


BASE_DIR = PROJECT_ROOT / "04_processed/discovery_v0_7_public_validation"
HNSCC_GBM_DIR = PROJECT_ROOT / "04_processed/discovery_v0_8_public_validation/hnscc_gbm_incremental"
THREECA_DIR = PROJECT_ROOT / "04_processed/discovery_v0_8_public_validation/threeca_incremental"
OUT_DIR = PROJECT_ROOT / "04_processed/discovery_v0_8_public_validation"
FIG_DIR = PROJECT_ROOT / "06_figures/discovery_v0_8_public_validation"
REPORT_PATH = PROJECT_ROOT / "05_reports/public_source_mapping_v0_8_report.md"


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def concat_frames(paths: list[Path]) -> pd.DataFrame:
    frames = [read_csv(path) for path in paths]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).drop_duplicates()


def build_pairwise_inference(score_df: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("immune", "malignant"),
        ("immune", "stromal"),
        ("malignant", "stromal"),
    ]
    rows = []
    group_cols = ["source_id", "dataset_title", "cancer_code", "signature", "signature_label"]
    for keys, group in score_df.groupby(group_cols, sort=True):
        pivot = group.pivot_table(
            index="donor_id",
            columns="compartment",
            values="signature_score",
            aggfunc="mean",
        )
        for case, reference in pairs:
            if case not in pivot.columns or reference not in pivot.columns:
                continue
            delta = (pivot[case] - pivot[reference]).dropna().to_numpy(dtype=float)
            if delta.size < 2:
                continue
            nonzero = delta[delta != 0]
            n_positive = int(np.sum(delta > 0))
            sign_p = (
                float(
                    stats.binomtest(
                        int(np.sum(nonzero > 0)),
                        int(nonzero.size),
                        p=0.5,
                        alternative="two-sided",
                    ).pvalue
                )
                if nonzero.size
                else 1.0
            )
            wilcoxon_p = (
                float(stats.wilcoxon(delta, alternative="two-sided", method="auto").pvalue)
                if nonzero.size
                else 1.0
            )
            seed_text = "|".join(map(str, keys)) + f"|{case}|{reference}"
            rng = np.random.default_rng(zlib.crc32(seed_text.encode("utf-8")))
            draws = rng.choice(delta, size=(10_000, delta.size), replace=True).mean(axis=1)
            sd = float(np.std(delta, ddof=1))
            rows.append(
                {
                    "source_id": keys[0],
                    "dataset_title": keys[1],
                    "cancer_code": keys[2],
                    "signature": keys[3],
                    "signature_label": keys[4],
                    "contrast": f"{case}_minus_{reference}",
                    "paired_donor_n": int(delta.size),
                    "mean_delta": float(np.mean(delta)),
                    "mean_delta_bootstrap_ci95_low": float(np.quantile(draws, 0.025)),
                    "mean_delta_bootstrap_ci95_high": float(np.quantile(draws, 0.975)),
                    "median_delta": float(np.median(delta)),
                    "paired_standardized_mean_difference_dz": (
                        float(np.mean(delta) / sd) if sd > 0 else np.nan
                    ),
                    "positive_donor_n": n_positive,
                    "positive_fraction": float(n_positive / delta.size),
                    "two_sided_sign_p": sign_p,
                    "two_sided_wilcoxon_p": wilcoxon_p,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["sign_fdr_global"] = stats.false_discovery_control(
        out["two_sided_sign_p"].to_numpy(dtype=float), method="bh"
    )
    out["wilcoxon_fdr_global"] = stats.false_discovery_control(
        out["two_sided_wilcoxon_p"].to_numpy(dtype=float), method="bh"
    )
    out["sign_fdr_within_source_cancer"] = np.nan
    out["wilcoxon_fdr_within_source_cancer"] = np.nan
    for _, idx in out.groupby(["source_id", "cancer_code"]).groups.items():
        loc = list(idx)
        out.loc[loc, "sign_fdr_within_source_cancer"] = stats.false_discovery_control(
            out.loc[loc, "two_sided_sign_p"].to_numpy(dtype=float), method="bh"
        )
        out.loc[loc, "wilcoxon_fdr_within_source_cancer"] = stats.false_discovery_control(
            out.loc[loc, "two_sided_wilcoxon_p"].to_numpy(dtype=float), method="bh"
        )
    return out


def plot_public_heatmap(summary_df: pd.DataFrame, out_path: Path) -> None:
    plot_df = summary_df.copy()
    plot_df["value"] = plot_df["median_immune"] - plot_df["median_malignant"]
    plot_df["cohort_label"] = plot_df.apply(
        lambda row: (
            f"{str(row['cancer_code']).upper()}\n"
            + (
                "3CA-Peng2019"
                if row["source_id"] == "3ca-peng2019-pdac"
                else "3CA-Regner2021"
                if row["source_id"] == "3ca-regner2021-ucec"
                else f"CXG-{str(row['source_id'])[:4]}"
            )
        ),
        axis=1,
    )
    cohort_order = (
        plot_df[["cancer_code", "source_id", "cohort_label"]]
        .drop_duplicates()
        .sort_values(["cancer_code", "source_id"])["cohort_label"]
        .tolist()
    )
    signature_order = [
        "Checkpoint",
        "Cytolytic",
        "IFN-gamma",
        "MHC-I",
        "MHC-II",
        "Myeloid",
        "Proteasome",
        "TGF/EMT",
    ]
    matrix = plot_df.pivot_table(
        index="signature_label",
        columns="cohort_label",
        values="value",
        aggfunc="median",
    ).reindex(index=signature_order, columns=cohort_order)
    finite = matrix.to_numpy(dtype=float)
    vmax = float(np.nanquantile(np.abs(finite), 0.98))
    vmax = max(vmax, 0.25)

    fig_width = max(12.0, 0.82 * len(cohort_order))
    fig, ax = plt.subplots(figsize=(fig_width, 5.8), constrained_layout=True)
    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(cohort_order)))
    ax.set_xticklabels(cohort_order, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(signature_order)))
    ax.set_yticklabels(signature_order)
    ax.set_title("Immune-minus-malignant pseudobulk signature medians across public cohorts")
    ax.set_xlabel("Cancer and source cohort")
    ax.set_ylabel("Signature")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Median signature-score difference")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_ifngamma_forest(inference_df: pd.DataFrame, out_path: Path) -> None:
    forest = inference_df.loc[
        (inference_df["signature_label"] == "IFN-gamma")
        & (inference_df["contrast"] == "immune_minus_stromal")
    ].copy()
    if forest.empty:
        return
    forest["cohort_label"] = forest.apply(
        lambda row: (
            f"{str(row['cancer_code']).upper()} | "
            + (
                "3CA-Peng2019"
                if row["source_id"] == "3ca-peng2019-pdac"
                else "3CA-Regner2021"
                if row["source_id"] == "3ca-regner2021-ucec"
                else f"CXG-{str(row['source_id'])[:4]}"
            )
            + f" | n={int(row['paired_donor_n'])}"
        ),
        axis=1,
    )
    forest = forest.sort_values(["mean_delta", "cancer_code"]).reset_index(drop=True)
    y = np.arange(forest.shape[0])
    colors = np.where(forest["wilcoxon_fdr_global"] < 0.05, "#b63a3a", "#547a91")
    xerr = np.vstack(
        [
            forest["mean_delta"] - forest["mean_delta_bootstrap_ci95_low"],
            forest["mean_delta_bootstrap_ci95_high"] - forest["mean_delta"],
        ]
    )
    fig_height = max(5.0, 0.42 * forest.shape[0] + 1.8)
    fig, ax = plt.subplots(figsize=(8.8, fig_height), constrained_layout=True)
    for idx in range(forest.shape[0]):
        ax.errorbar(
            forest.loc[idx, "mean_delta"],
            y[idx],
            xerr=xerr[:, idx : idx + 1],
            fmt="o",
            color=colors[idx],
            ecolor=colors[idx],
            capsize=3,
            markersize=6,
        )
    ax.axvline(0.0, color="#333333", linewidth=1, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(forest["cohort_label"])
    ax.set_xlabel("Mean immune-minus-stromal IFN-gamma score difference (bootstrap 95% CI)")
    ax.set_title("Donor-paired IFN-gamma compartment contrast across public cohorts")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_report(
    schema_df: pd.DataFrame,
    group_df: pd.DataFrame,
    pair_df: pd.DataFrame,
    inference_df: pd.DataFrame,
) -> None:
    source_ids = sorted(schema_df["source_id"].dropna().unique().tolist()) if not schema_df.empty else []
    cancer_codes = sorted(group_df["cancer_code"].dropna().unique().tolist()) if not group_df.empty else []
    cellxgene_count = sum(1 for x in source_ids if not str(x).startswith("3ca-"))
    threeca_count = sum(1 for x in source_ids if str(x).startswith("3ca-"))
    lines = [
        "# BIO-05 Public Source Mapping v0.8",
        "",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        f"- Platform: {platform.platform()}",
        f"- Official public datasets with analyzable pseudobulk groups: {len(source_ids)}",
        f"- Source routes: {cellxgene_count} CELLxGENE H5ADs + {threeca_count} 3CA tarball matrix packages",
        f"- Cancer codes covered: {', '.join(cancer_codes) if cancer_codes else 'none'}",
        "",
        "## Schema Audit",
        "",
    ]
    for source_id, sub in schema_df.groupby("source_id", dropna=False):
        dataset_title = sub["dataset_title"].iloc[0]
        lines.append(f"- `{source_id}` `{dataset_title}`")
        for item in sub.itertuples():
            field_name = "" if pd.isna(item.field_name) else str(item.field_name)
            top_values = "" if pd.isna(item.top_values) else str(item.top_values)
            lines.append(
                f"  {item.field_role}: `{field_name}`; nunique={item.nunique}; top={top_values[:160]}"
            )

    lines.extend(["", "## Pairwise Contrast Snapshot", ""])
    if pair_df.empty:
        lines.append("- No pairwise donor-aware contrast was available.")
    else:
        top = pair_df.sort_values(["signature_label", "contrast", "paired_donor_n"], ascending=[True, True, False])
        for row in top.head(36).itertuples():
            lines.append(
                f"- `{row.cancer_code}` `{row.signature_label}` `{row.contrast}`: n={row.paired_donor_n}, median_delta={row.median_delta:.3f}, positive_fraction={row.positive_fraction:.3f}"
            )
    lines.extend(["", "## Donor-Paired Inference", ""])
    if inference_df.empty:
        lines.append("- No donor-paired inferential contrast was available.")
    else:
        ranked = inference_df.sort_values(
            ["wilcoxon_fdr_global", "paired_donor_n"],
            ascending=[True, False],
        )
        for row in ranked.head(24).itertuples():
            lines.append(
                f"- `{row.cancer_code}` `{row.signature_label}` `{row.contrast}`: "
                f"n={row.paired_donor_n}, mean_delta={row.mean_delta:.3f} "
                f"(bootstrap 95% CI {row.mean_delta_bootstrap_ci95_low:.3f} to "
                f"{row.mean_delta_bootstrap_ci95_high:.3f}), dz={row.paired_standardized_mean_difference_dz:.3f}, "
                f"Wilcoxon P={row.two_sided_wilcoxon_p:.4g}, global FDR={row.wilcoxon_fdr_global:.4g}"
            )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)

    schema_df = concat_frames(
        [
            BASE_DIR / "scrna_schema_audit_v0_7.csv",
            HNSCC_GBM_DIR / "scrna_schema_audit_v0_7.csv",
            THREECA_DIR / "threeca_schema_audit_v0_8.csv",
        ]
    )
    group_df = concat_frames(
        [
            BASE_DIR / "scrna_group_counts_v0_7.csv",
            HNSCC_GBM_DIR / "scrna_group_counts_v0_7.csv",
            THREECA_DIR / "threeca_group_counts_v0_8.csv",
        ]
    )
    score_df = concat_frames(
        [
            BASE_DIR / "scrna_signature_scores_v0_7.csv",
            HNSCC_GBM_DIR / "scrna_signature_scores_v0_7.csv",
            THREECA_DIR / "threeca_signature_scores_v0_8.csv",
        ]
    )
    pair_df = concat_frames(
        [
            BASE_DIR / "scrna_pairwise_contrasts_v0_7.csv",
            HNSCC_GBM_DIR / "scrna_pairwise_contrasts_v0_7.csv",
            THREECA_DIR / "threeca_pairwise_contrasts_v0_8.csv",
        ]
    )
    summary_df = concat_frames(
        [
            BASE_DIR / "scrna_dataset_cancer_signature_summary_v0_7.csv",
            HNSCC_GBM_DIR / "scrna_dataset_cancer_signature_summary_v0_7.csv",
            THREECA_DIR / "threeca_dataset_cancer_signature_summary_v0_8.csv",
        ]
    )
    inference_df = build_pairwise_inference(score_df)

    schema_df.to_csv(OUT_DIR / f"public_source_schema_audit_{VERSION}.csv", index=False)
    group_df.to_csv(OUT_DIR / f"public_source_group_counts_{VERSION}.csv", index=False)
    score_df.to_csv(OUT_DIR / f"public_source_signature_scores_{VERSION}.csv", index=False)
    pair_df.to_csv(OUT_DIR / f"public_source_pairwise_contrasts_{VERSION}.csv", index=False)
    summary_df.to_csv(OUT_DIR / f"public_source_dataset_cancer_signature_summary_{VERSION}.csv", index=False)
    inference_df.to_csv(OUT_DIR / f"public_source_pairwise_inference_{VERSION}.csv", index=False)

    plot_public_heatmap(
        summary_df,
        FIG_DIR / f"fig_public_immune_minus_malignant_heatmap_{VERSION}.png",
    )
    plot_ifngamma_forest(
        inference_df,
        FIG_DIR / f"fig_public_ifngamma_immune_minus_stromal_forest_{VERSION}.png",
    )
    write_report(schema_df, group_df, pair_df, inference_df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
