#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import math
import os
import platform
from datetime import datetime
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


VERSION = "v0_7"
MIN_SIGNATURE_GENES = 3
DATA_ROOT = Path(os.environ["DATA_WORKSPACE"]) / "14_pan_cancer_proteogenomic_immune_evasion"

SIGNATURE_LABELS = {
    "antigen_presentation_mhc_i": "MHC-I",
    "antigen_presentation_mhc_ii": "MHC-II",
    "checkpoint_exhaustion_context": "Checkpoint",
    "cytolytic_t_cell_context": "Cytolytic",
    "ifn_gamma_response": "IFN-gamma",
    "myeloid_inflammatory_context": "Myeloid",
    "proteasome_antigen_processing": "Proteasome",
    "tgfb_emt_exclusion_context": "TGF/EMT",
}

CANCER_SCORE_COLS = [
    "Cancer Basal SC",
    "Cancer Cycling",
    "Cancer Her2 SC",
    "Cancer LumA SC",
    "Cancer LumB SC",
]

IMMUNE_SCORE_COLS = [
    "B cells Memory",
    "B cells Naive",
    "Cycling T-cells",
    "Cycling_Myeloid",
    "DCs",
    "Macrophage",
    "Monocyte",
    "NK cells",
    "NKT cells",
    "Plasmablasts",
    "T cells CD4+",
    "T cells CD8+",
]

STROMAL_SCORE_COLS = [
    "CAFs MSC iCAF-like",
    "CAFs myCAF-like",
    "Endothelial ACKR1",
    "Endothelial CXCL12",
    "Endothelial Lymphatic LYVE1",
    "Endothelial RGS5",
    "Myoepithelial",
    "PVL Differentiated",
    "PVL Immature",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="BIO-05 BRCA spatial triangulation.")
    parser.add_argument(
        "--queue",
        default=expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/manifests/BIO-05/bio05_brca_spatial_queue_20260730.tsv"),
    )
    parser.add_argument(
        "--signature-table",
        default=str(DATA_ROOT / "04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv"),
    )
    parser.add_argument(
        "--out-dir",
        default=str(DATA_ROOT / "04_processed/discovery_v0_7_public_validation"),
    )
    parser.add_argument(
        "--figure-dir",
        default=str(DATA_ROOT / "06_figures/discovery_v0_7_public_validation"),
    )
    parser.add_argument(
        "--report",
        default=str(DATA_ROOT / "05_reports/public_brca_spatial_triangulation_v0_7_report.md"),
    )
    return parser.parse_args()


def safe_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def normalize_text(value: Any) -> str:
    return safe_str(value).strip().lower()


def load_signature_genes(path: Path) -> dict[str, list[str]]:
    table = pd.read_csv(path)
    return {
        signature: sorted({safe_str(g).strip().upper() for g in group["gene_symbol"] if safe_str(g).strip()})
        for signature, group in table.groupby("signature", dropna=False)
    }


def file_is_complete(path: Path, expected_size: str) -> bool:
    if not path.exists():
        return False
    expected_text = safe_str(expected_size).strip()
    if not expected_text:
        return path.stat().st_size > 0
    expected_value = int(float(expected_text))
    return path.stat().st_size >= expected_value


def load_queue(path: Path) -> pd.DataFrame:
    queue = pd.read_csv(path, sep="\t", dtype=str)
    queue["remote_path"] = queue["remote_path"].map(Path)
    queue["observed_size_bytes"] = queue["remote_path"].map(lambda p: p.stat().st_size if p.exists() else -1)
    queue["is_complete"] = [
        file_is_complete(path, expected)
        for path, expected in zip(queue["remote_path"], queue["expected_size"])
    ]
    return queue


def pick_gene_symbols(adata: ad.AnnData) -> pd.Index:
    if "feature_name" in adata.var.columns:
        return adata.var["feature_name"].astype(str).str.upper()
    return adata.var_names.astype(str).str.upper()


def cpm_log1p(matrix: sparse.csr_matrix) -> sparse.csr_matrix:
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    totals[totals == 0] = np.nan
    scale = np.divide(1e6, totals, out=np.zeros_like(totals), where=np.isfinite(totals))
    scaled = sparse.diags(scale) @ matrix
    scaled.data = np.log1p(scaled.data)
    return scaled.tocsr()


def zscore_dense(values: np.ndarray) -> np.ndarray:
    out = values.astype(float).copy()
    means = np.nanmean(out, axis=0)
    stds = np.nanstd(out, axis=0, ddof=0)
    for idx in range(out.shape[1]):
        std = stds[idx]
        if not np.isfinite(std) or std == 0:
            out[:, idx] = np.nan
        else:
            out[:, idx] = (out[:, idx] - means[idx]) / std
    return out


def assign_spot_compartment(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for label, cols in [
        ("cancer_score", CANCER_SCORE_COLS),
        ("immune_score", IMMUNE_SCORE_COLS),
        ("stromal_score", STROMAL_SCORE_COLS),
    ]:
        present = [col for col in cols if col in out.columns]
        if not present:
            out[label] = 0.0
        else:
            out[label] = out[present].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=1)

    compartments: list[str] = []
    rules: list[str] = []
    for row in out.itertuples():
        classification = normalize_text(getattr(row, "Classification", ""))
        scores = {
            "malignant": float(getattr(row, "cancer_score")),
            "immune": float(getattr(row, "immune_score")),
            "stromal": float(getattr(row, "stromal_score")),
        }
        best_compartment = max(scores, key=scores.get)
        best_value = scores[best_compartment]
        if best_value > 0:
            compartments.append(best_compartment)
            rules.append("spot_score_argmax")
            continue
        if classification == "stroma":
            compartments.append("stromal")
            rules.append("classification_stroma")
            continue
        if "lymphocytes" in classification:
            compartments.append("immune")
            rules.append("classification_lymphocytes")
            continue
        if "invasive cancer" in classification:
            compartments.append("malignant")
            rules.append("classification_invasive_cancer")
            continue
        compartments.append("other")
        rules.append("unclassified")
    out["spot_compartment"] = compartments
    out["compartment_rule"] = rules
    return out


def process_spatial_dataset(
    row: pd.Series,
    signature_genes: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = Path(row["remote_path"])
    adata = ad.read_h5ad(path)
    obs = adata.obs.copy()
    obs["_row_index"] = np.arange(adata.n_obs)
    obs["spot_id"] = obs.index.astype(str)
    obs["sample_id"] = safe_str(obs["donor_id"].iloc[0]) if "donor_id" in obs.columns else safe_str(row["title"])
    obs["source_id"] = safe_str(row["source_id"])
    obs["dataset_title"] = safe_str(row["title"])

    if "in_tissue" in obs.columns:
        in_tissue = pd.to_numeric(obs["in_tissue"], errors="coerce").fillna(0).astype(int)
        obs = obs.loc[in_tissue.eq(1)].copy()
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    obs = assign_spot_compartment(obs)
    obs = obs[obs["spot_compartment"].isin(["malignant", "immune", "stromal"])].copy()
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    matrix = adata.X
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(np.asarray(matrix))
    else:
        matrix = matrix.tocsr()
    matrix = matrix[obs["_row_index"].to_numpy(), :]
    matrix = cpm_log1p(matrix)

    gene_symbols = pick_gene_symbols(adata)
    gene_to_idx: dict[str, list[int]] = {}
    for idx, gene in enumerate(gene_symbols):
        gene_to_idx.setdefault(gene, []).append(idx)

    score_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for signature, genes in signature_genes.items():
        recovered_idx: list[int] = []
        recovered_genes: list[str] = []
        for gene in genes:
            if gene in gene_to_idx:
                recovered_idx.extend(gene_to_idx[gene])
                recovered_genes.append(gene)
        if len(recovered_genes) < MIN_SIGNATURE_GENES:
            continue
        sub = matrix[:, recovered_idx]
        dense = sub.toarray() if sparse.issparse(sub) else np.asarray(sub)
        z = zscore_dense(dense)
        scores = np.nanmean(z, axis=1)
        score_df = obs[["spot_id", "sample_id", "spot_compartment", "compartment_rule", "Classification", "annotation"]].copy()
        score_df["signature"] = signature
        score_df["signature_label"] = SIGNATURE_LABELS.get(signature, signature)
        score_df["signature_score"] = scores
        score_df["n_genes_used"] = len(recovered_genes)
        score_df["source_id"] = safe_str(row["source_id"])
        score_df["dataset_title"] = safe_str(row["title"])
        score_rows.extend(score_df.to_dict(orient="records"))

        med = score_df.groupby("spot_compartment")["signature_score"].median().to_dict()
        summary_rows.append(
            {
                "source_id": safe_str(row["source_id"]),
                "dataset_title": safe_str(row["title"]),
                "sample_id": safe_str(score_df["sample_id"].iloc[0]),
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "n_spots": int(score_df["spot_id"].nunique()),
                "n_genes_used": len(recovered_genes),
                "median_malignant": med.get("malignant"),
                "median_immune": med.get("immune"),
                "median_stromal": med.get("stromal"),
            }
        )

    score_long = pd.DataFrame(score_rows)
    summary_df = pd.DataFrame(summary_rows)
    counts = (
        obs.groupby(["sample_id", "spot_compartment", "compartment_rule"], dropna=False)
        .size()
        .reset_index(name="n_spots")
    )
    counts["source_id"] = safe_str(row["source_id"])
    counts["dataset_title"] = safe_str(row["title"])
    return counts, score_long, summary_df


def summarize_pairwise(summary_df: pd.DataFrame) -> pd.DataFrame:
    if summary_df.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for signature_label, sub in summary_df.groupby("signature_label", dropna=False):
        for left, right in [("malignant", "immune"), ("malignant", "stromal"), ("immune", "stromal")]:
            left_col = f"median_{left}"
            right_col = f"median_{right}"
            delta = pd.to_numeric(sub[left_col] - sub[right_col], errors="coerce").dropna()
            if delta.empty:
                continue
            rows.append(
                {
                    "signature_label": signature_label,
                    "contrast": f"{left}_minus_{right}",
                    "paired_sample_n": int(delta.shape[0]),
                    "median_delta": float(delta.median()),
                    "mean_delta": float(delta.mean()),
                    "positive_fraction": float((delta > 0).mean()),
                }
            )
    return pd.DataFrame(rows)


def summarize_pathology(score_long: pd.DataFrame) -> pd.DataFrame:
    if score_long.empty:
        return pd.DataFrame()
    unique_spots = score_long[["sample_id", "Classification", "spot_compartment", "spot_id"]].drop_duplicates().copy()
    out = (
        unique_spots.groupby(["sample_id", "Classification", "spot_compartment"], dropna=False)
        .size()
        .reset_index(name="n_spots")
        .sort_values(["sample_id", "Classification", "spot_compartment"])
    )
    return out


def plot_heatmap(summary_df: pd.DataFrame, out_path: Path) -> None:
    if summary_df.empty:
        return
    plot_df = summary_df.copy()
    frames = []
    for compartment in ["malignant", "immune", "stromal"]:
        temp = plot_df[["sample_id", "signature_label", f"median_{compartment}"]].copy()
        temp["sample_compartment"] = temp["sample_id"] + "|" + compartment
        temp["median_score"] = temp[f"median_{compartment}"]
        frames.append(temp[["signature_label", "sample_compartment", "median_score"]])
    heat_df = pd.concat(frames, ignore_index=True)
    heat = heat_df.pivot(index="signature_label", columns="sample_compartment", values="median_score")
    if heat.empty:
        return
    fig, ax = plt.subplots(figsize=(max(8, 0.8 * heat.shape[1]), 5))
    im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="coolwarm")
    ax.set_xticks(np.arange(heat.shape[1]))
    ax.set_xticklabels(heat.columns, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(heat.shape[0]))
    ax.set_yticklabels(heat.index, fontsize=9)
    ax.set_title("BRCA spatial compartment medians")
    fig.colorbar(im, ax=ax, shrink=0.8, label="median signature z-score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def render_report(
    queue_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    pairwise_df: pd.DataFrame,
    report_path: Path,
) -> None:
    lines = [
        "# BIO-05 BRCA Spatial Triangulation v0.7",
        "",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        f"- Platform: {platform.platform()}",
        f"- Completed spatial files audited: {int(queue_df['is_complete'].sum())}",
        f"- Samples with analyzable spatial compartments: {counts_df['sample_id'].nunique() if not counts_df.empty else 0}",
        "",
        "## Completed Files",
        "",
    ]
    for row in queue_df[queue_df["is_complete"]].itertuples():
        lines.append(
            f"- `{row.source_id}`: `{row.title}` -> `{row.remote_path}` ({int(row.observed_size_bytes)} bytes)"
        )
    lines.extend(["", "## Pairwise Contrast Snapshot", ""])
    if pairwise_df.empty:
        lines.append("- No analyzable pairwise spatial contrasts yet.")
    else:
        for row in pairwise_df.sort_values(["signature_label", "contrast"]).itertuples():
            lines.append(
                f"- `{row.signature_label}` `{row.contrast}`: n={row.paired_sample_n}, median_delta={row.median_delta:.3f}, positive_fraction={row.positive_fraction:.3f}"
            )
    report_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = Path(args.figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_path = Path(args.report)

    queue_df = load_queue(Path(args.queue))
    signature_genes = load_signature_genes(Path(args.signature_table))

    count_frames: list[pd.DataFrame] = []
    score_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    for row in queue_df[queue_df["is_complete"]].itertuples(index=False):
        counts_df, score_long, summary_df = process_spatial_dataset(pd.Series(row._asdict()), signature_genes)
        if not counts_df.empty:
            count_frames.append(counts_df)
        if not score_long.empty:
            score_frames.append(score_long)
        if not summary_df.empty:
            summary_frames.append(summary_df)

    counts_df = pd.concat(count_frames, ignore_index=True) if count_frames else pd.DataFrame()
    score_long = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    summary_df = pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame()
    pairwise_df = summarize_pairwise(summary_df)
    pathology_df = summarize_pathology(score_long)

    queue_df.to_csv(out_dir / f"spatial_brca_queue_audit_{VERSION}.csv", index=False)
    counts_df.to_csv(out_dir / f"spatial_brca_spot_compartment_counts_{VERSION}.csv", index=False)
    score_long.to_csv(out_dir / f"spatial_brca_signature_scores_{VERSION}.csv", index=False)
    summary_df.to_csv(out_dir / f"spatial_brca_sample_signature_summary_{VERSION}.csv", index=False)
    pairwise_df.to_csv(out_dir / f"spatial_brca_pairwise_contrasts_{VERSION}.csv", index=False)
    pathology_df.to_csv(out_dir / f"spatial_brca_pathology_alignment_{VERSION}.csv", index=False)
    plot_heatmap(summary_df, figure_dir / f"fig_brca_spatial_compartment_heatmap_{VERSION}.png")
    render_report(queue_df, counts_df, pairwise_df, report_path)


if __name__ == "__main__":
    main()
