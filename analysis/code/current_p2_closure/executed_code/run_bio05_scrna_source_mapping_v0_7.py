#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import json
import math
import platform
import re
import time
from collections import Counter
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


PROJECT_ROOT = Path(
    expand_legacy_paths("${SCI27_DATA_ROOT}/14_pan_cancer_proteogenomic_immune_evasion")
)
VERSION = "v0_7"
MIN_GROUP_CELLS = 20
MIN_SIGNATURE_GENES = 3

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

DATASET_DEFAULT_CANCER = {
    "3f7c572c-cd73-4b51-a313-207c7f20f188": None,
    "dea97145-f712-431c-a223-6b5f565f362a": "brca",
    "1df8c90d-d299-4b2e-a54d-a5a80f36e780": "ccrcc",
    "3f50314f-bdc9-40c6-8e4a-b0901ebfbe4c": "ccrcc",
    "4d82bd8e-8827-410b-8fc7-685b7dd92585": "coad",
    "0bebef1a-4607-4584-9070-dacf89a0d635": "luad",
    "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": "hnscc",
    "999f2a15-3d7e-440b-96ae-2c806799c08c": "gbm",
}

IMMUNE_KEYWORDS = {
    "immune cell",
    "b cell",
    "plasma cell",
    "t cell",
    "nk",
    "natural killer",
    "macrophage",
    "monocyte",
    "microgl",
    "dendritic",
    "myeloid",
    "mast cell",
    "neutrophil",
    "granulocyte",
    "lymphocyte",
    "megakaryocyte",
}

STROMAL_KEYWORDS = {
    "fibroblast",
    "myofibroblast",
    "endothelial",
    "pericyte",
    "smooth muscle",
    "mural cell",
    "vascular associated smooth muscle cell",
    "mesangial",
    "podocyte",
    "schwann",
    "stromal",
    "astrocyte",
    "oligodendrocyte",
    "oligodendrocyte precursor",
    "radial glial",
    "neuron",
    "support cell",
}

MALIGNANT_KEYWORDS = {
    "malignant",
    "abnormal cell",
    "tumor",
    "tumour",
    "cancer epithelial",
    "neoplastic",
}

DATASET_SPECIFIC_MALIGNANT = {
    "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": {"squamous epithelial cell"},
    "0bebef1a-4607-4584-9070-dacf89a0d635": {"epithelial cell of lung"},
}

DATASET_SPECIFIC_OTHER = {
    "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": {"salivary gland glandular cell"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="BIO-05 public scRNA source mapping with donor-aware pseudobulk."
    )
    parser.add_argument(
        "--queue",
        default=expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/manifests/BIO-05/bio05_public_validation_queue_20260730.tsv"),
    )
    parser.add_argument(
        "--signature-table",
        default="04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="04_processed/discovery_v0_7_public_validation",
    )
    parser.add_argument(
        "--figure-dir",
        default="06_figures/discovery_v0_7_public_validation",
    )
    parser.add_argument(
        "--report",
        default="05_reports/public_scrna_source_mapping_v0_7_report.md",
    )
    parser.add_argument("--min-group-cells", type=int, default=MIN_GROUP_CELLS)
    parser.add_argument("--min-signature-genes", type=int, default=MIN_SIGNATURE_GENES)
    return parser.parse_args()


def safe_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def normalize_text(value: Any) -> str:
    return re.sub(r"\s+", " ", safe_str(value).strip().lower())


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
    if expected_text:
        try:
            expected_value = int(float(expected_text))
        except ValueError:
            expected_value = None
        if expected_value is not None and expected_value > 0:
            return path.stat().st_size >= expected_value
    return path.stat().st_size > 0


def load_queue(path: Path) -> pd.DataFrame:
    queue = pd.read_csv(path, sep="\t", dtype=str)
    queue["remote_path"] = queue["remote_path"].map(Path)
    return queue


def pick_gene_symbols(adata: ad.AnnData) -> pd.Index:
    if adata.raw is not None:
        raw_var = adata.raw.var.copy()
        if "feature_name" in raw_var.columns:
            return raw_var["feature_name"].astype(str).str.upper()
        return raw_var.index.astype(str).str.upper()
    if "feature_name" in adata.var.columns:
        return adata.var["feature_name"].astype(str).str.upper()
    return adata.var_names.astype(str).str.upper()


def candidate_obs_columns(adata: ad.AnnData, pattern: str) -> list[str]:
    cols = []
    for col in adata.obs.columns:
        if re.search(pattern, col, flags=re.I):
            cols.append(col)
    return cols


def best_obs_column(adata: ad.AnnData, preferred: list[str], pattern: str) -> str | None:
    for col in preferred:
        if col in adata.obs.columns:
            return col
    matches = candidate_obs_columns(adata, pattern)
    return matches[0] if matches else None


def infer_cancer_code(disease: str, tissue: str, default_code: str | None) -> str | None:
    disease_text = normalize_text(disease)
    tissue_text = normalize_text(tissue)
    if disease_text in {"normal", "healthy"}:
        return None
    text = f"{disease_text} | {tissue_text}"
    if "clear cell renal" in text or "renal cell carcinoma" in text or "kidney cancer" in text:
        return "ccrcc"
    if "glioblastoma" in text:
        return "gbm"
    if "head and neck" in text or "oropharynx squamous" in text or "oral cavity squamous" in text:
        return "hnscc"
    if "lung adenocarcinoma" in text:
        return "luad"
    if "squamous cell lung" in text:
        return "lscc"
    if "squamous cell carcinoma" in text and "lung" in text:
        return "lscc"
    if "breast" in text:
        return "brca"
    if "ovarian" in text or "serous adenocarcinoma" in text:
        return "ov"
    if "colon adenocarcinoma" in text or "colorectal" in text or "cecum adenocarcinoma" in text or "rectosigmoid" in text:
        return "coad"
    return default_code


def classify_compartment(
    row: pd.Series,
    source_id: str,
) -> tuple[str, str]:
    cell_type = normalize_text(row.get("cell_type"))
    celltype_major = normalize_text(row.get("celltype_major"))
    author_cell_type = normalize_text(row.get("author_cell_type"))
    author_cluster = normalize_text(row.get("author_cluster_label"))
    normal_cell_call = normalize_text(row.get("normal_cell_call"))
    combined = " | ".join([cell_type, celltype_major, author_cell_type, author_cluster])

    if normal_cell_call == "cancer":
        return "malignant", "explicit_normal_cell_call_cancer"
    if normal_cell_call == "normal" and "epithelial" in combined:
        return "other", "explicit_normal_cell_call_normal_epithelial"

    if celltype_major == "cancer epithelial":
        return "malignant", "explicit_celltype_major_cancer_epithelial"
    if celltype_major in {"t-cells", "myeloid", "plasmablasts", "b-cells"}:
        return "immune", "explicit_celltype_major_immune"
    if celltype_major in {"cafs", "endothelial", "pvl"}:
        return "stromal", "explicit_celltype_major_stromal"
    if celltype_major == "normal epithelial":
        return "other", "explicit_celltype_major_normal_epithelial"

    other_terms = {normalize_text(x) for x in DATASET_SPECIFIC_OTHER.get(source_id, set())}
    if cell_type in other_terms:
        return "other", "dataset_specific_other"

    malignant_terms = {normalize_text(x) for x in DATASET_SPECIFIC_MALIGNANT.get(source_id, set())}
    if cell_type in malignant_terms:
        return "malignant", "dataset_specific_malignant"
    if any(term in combined for term in MALIGNANT_KEYWORDS):
        return "malignant", "keyword_malignant"
    if any(term in combined for term in IMMUNE_KEYWORDS):
        return "immune", "keyword_immune"
    if any(term in combined for term in STROMAL_KEYWORDS):
        return "stromal", "keyword_stromal"
    if "epithelial" in cell_type and source_id in {
        "dea97145-f712-431c-a223-6b5f565f362a",
        "4d82bd8e-8827-410b-8fc7-685b7dd92585",
    }:
        return "malignant", "tumor_epithelial_fallback"
    return "other", "unmapped"


def build_group_mean_matrix(matrix: sparse.spmatrix, group_codes: np.ndarray, group_count: int) -> sparse.csr_matrix:
    group_mat = sparse.csr_matrix(
        (np.ones(len(group_codes), dtype=float), (group_codes, np.arange(len(group_codes)))),
        shape=(group_count, len(group_codes)),
    )
    sizes = np.asarray(group_mat.sum(axis=1)).ravel()
    weights = sparse.diags(1.0 / np.maximum(sizes, 1.0))
    return weights @ group_mat @ matrix


def cpm_log1p(matrix: sparse.csr_matrix) -> sparse.csr_matrix:
    libsize = np.asarray(matrix.sum(axis=1)).ravel()
    scaled = sparse.diags(1e6 / np.maximum(libsize, 1.0)) @ matrix
    scaled = scaled.tocsr(copy=True)
    scaled.data = np.log1p(scaled.data)
    return scaled


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for col in out.columns:
        values = pd.to_numeric(out[col], errors="coerce")
        std = float(values.std(ddof=0))
        if not np.isfinite(std) or std == 0:
            out[col] = np.nan
        else:
            out[col] = (values - float(values.mean())) / std
    return out


def summarize_schema(
    source_id: str,
    row: pd.Series,
    adata: ad.AnnData,
    donor_col: str | None,
    cell_type_col: str | None,
    disease_col: str | None,
    tissue_col: str | None,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for label, col in [
        ("donor", donor_col),
        ("cell_type", cell_type_col),
        ("disease", disease_col),
        ("tissue", tissue_col),
    ]:
        if not col:
            records.append(
                {
                    "source_id": source_id,
                    "dataset_title": safe_str(row.get("title")),
                    "field_role": label,
                    "field_name": "",
                    "nunique": 0,
                    "top_values": "",
                }
            )
            continue
        values = adata.obs[col].astype(str)
        top_values = values.value_counts().head(10)
        records.append(
            {
                "source_id": source_id,
                "dataset_title": safe_str(row.get("title")),
                "field_role": label,
                "field_name": col,
                "nunique": int(values.nunique()),
                "top_values": "; ".join(f"{idx}:{int(val)}" for idx, val in top_values.items()),
            }
        )
    return records


def process_scrna_dataset(
    row: pd.Series,
    signature_genes: dict[str, list[str]],
    min_group_cells: int,
    min_signature_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    path = Path(row["remote_path"])
    source_id = safe_str(row["source_id"])
    adata = ad.read_h5ad(path)

    donor_col = best_obs_column(adata, ["donor_id", "patient_id", "participant_id"], r"donor|patient|participant")
    cell_type_col = best_obs_column(adata, ["cell_type", "author_cell_type"], r"cell.?type|cluster|annotation")
    disease_col = best_obs_column(adata, ["disease"], r"disease|cancer|tumou?r")
    tissue_col = best_obs_column(adata, ["tissue"], r"tissue|organ")

    schema_rows = summarize_schema(source_id, row, adata, donor_col, cell_type_col, disease_col, tissue_col)
    if not donor_col or not cell_type_col:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), schema_rows

    obs = adata.obs.copy()
    obs["donor_id_resolved"] = obs[donor_col].astype(str)
    obs["cell_type_resolved"] = obs[cell_type_col].astype(str)
    obs["disease_resolved"] = obs[disease_col].astype(str) if disease_col else ""
    obs["tissue_resolved"] = obs[tissue_col].astype(str) if tissue_col else ""
    default_cancer = DATASET_DEFAULT_CANCER.get(source_id)
    obs["cancer_code"] = [
        infer_cancer_code(disease, tissue, default_cancer)
        for disease, tissue in zip(obs["disease_resolved"], obs["tissue_resolved"])
    ]
    comp = obs.apply(lambda x: classify_compartment(x, source_id), axis=1)
    obs["broad_compartment"] = [x[0] for x in comp]
    obs["compartment_rule"] = [x[1] for x in comp]

    obs = obs[
        obs["cancer_code"].notna()
        & obs["broad_compartment"].isin(["malignant", "immune", "stromal"])
        & obs["donor_id_resolved"].astype(str).ne("")
    ].copy()
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), schema_rows

    group_counts = (
        obs.groupby(["cancer_code", "donor_id_resolved", "broad_compartment", "compartment_rule"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    group_counts["source_id"] = source_id
    group_counts["dataset_title"] = safe_str(row.get("title"))

    keep = group_counts[group_counts["n_cells"] >= min_group_cells].copy()
    if keep.empty:
        return group_counts, pd.DataFrame(), pd.DataFrame(), schema_rows

    obs = obs.merge(
        keep[["cancer_code", "donor_id_resolved", "broad_compartment"]],
        on=["cancer_code", "donor_id_resolved", "broad_compartment"],
        how="inner",
    )
    obs["group_id"] = (
        obs["cancer_code"].astype(str)
        + "__"
        + obs["donor_id_resolved"].astype(str)
        + "__"
        + obs["broad_compartment"].astype(str)
    )

    matrix = adata.raw.X if adata.raw is not None else adata.X
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(np.asarray(matrix))
    else:
        matrix = matrix.tocsr()
    sub = matrix[obs.index.to_numpy(), :]
    groups = obs["group_id"].astype("category")
    pb = build_group_mean_matrix(sub, groups.cat.codes.to_numpy(), len(groups.cat.categories))
    pb = cpm_log1p(pb)

    gene_symbols = pick_gene_symbols(adata)
    expr = pd.DataFrame(pb.toarray(), index=groups.cat.categories, columns=gene_symbols)
    expr = expr.T.groupby(level=0).mean().T

    score_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for cancer_code, cancer_groups in obs.groupby("cancer_code", dropna=False):
        group_ids = sorted(cancer_groups["group_id"].unique())
        cancer_expr = expr.loc[group_ids].copy()
        cancer_z = zscore_columns(cancer_expr)
        for signature, genes in signature_genes.items():
            recovered = [g for g in genes if g in cancer_z.columns]
            if len(recovered) < min_signature_genes:
                continue
            scores = cancer_z[recovered].mean(axis=1)
            score_df = pd.DataFrame(
                {
                    "group_id": scores.index,
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "signature_score": scores.to_numpy(dtype=float),
                    "n_genes_used": len(recovered),
                    "source_id": source_id,
                    "dataset_title": safe_str(row.get("title")),
                    "cancer_code": cancer_code,
                }
            )
            parts = score_df["group_id"].str.split("__", expand=True)
            score_df["donor_id"] = parts[1]
            score_df["compartment"] = parts[2]
            score_rows.extend(score_df.to_dict(orient="records"))

            med = score_df.groupby("compartment")["signature_score"].median().to_dict()
            summary_rows.append(
                {
                    "source_id": source_id,
                    "dataset_title": safe_str(row.get("title")),
                    "cancer_code": cancer_code,
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "n_groups": int(score_df["group_id"].nunique()),
                    "n_donors": int(score_df["donor_id"].nunique()),
                    "n_genes_used": len(recovered),
                    "median_malignant": med.get("malignant"),
                    "median_immune": med.get("immune"),
                    "median_stromal": med.get("stromal"),
                }
            )

    score_long = pd.DataFrame(score_rows)
    dataset_summary = pd.DataFrame(summary_rows)
    if score_long.empty:
        return group_counts, score_long, dataset_summary, schema_rows

    pair_rows: list[dict[str, Any]] = []
    pairs = [("malignant", "immune"), ("malignant", "stromal"), ("immune", "stromal")]
    for (source_id_val, dataset_title, cancer_code, signature_label), sub_scores in score_long.groupby(
        ["source_id", "dataset_title", "cancer_code", "signature_label"],
        dropna=False,
    ):
        wide = sub_scores.pivot_table(
            index="donor_id",
            columns="compartment",
            values="signature_score",
            aggfunc="mean",
        )
        for left, right in pairs:
            if left not in wide.columns or right not in wide.columns:
                continue
            delta = pd.to_numeric(wide[left] - wide[right], errors="coerce").dropna()
            if delta.empty:
                continue
            pair_rows.append(
                {
                    "source_id": source_id_val,
                    "dataset_title": dataset_title,
                    "cancer_code": cancer_code,
                    "signature_label": signature_label,
                    "contrast": f"{left}_minus_{right}",
                    "paired_donor_n": int(delta.shape[0]),
                    "median_delta": float(delta.median()),
                    "mean_delta": float(delta.mean()),
                    "positive_fraction": float((delta > 0).mean()),
                }
            )
    return group_counts, score_long, pd.DataFrame(pair_rows), schema_rows


def plot_compartment_heatmap(summary_df: pd.DataFrame, out_path: Path) -> None:
    if summary_df.empty:
        return
    plot_df = summary_df.copy()
    plot_df["dataset_cancer"] = plot_df["source_id"] + "|" + plot_df["cancer_code"]
    plot_df["immune_minus_malignant"] = plot_df["median_immune"] - plot_df["median_malignant"]
    heat = plot_df.pivot(index="signature_label", columns="dataset_cancer", values="immune_minus_malignant")
    if heat.empty:
        return
    fig_w = max(8, 0.6 * heat.shape[1])
    fig_h = max(5, 0.5 * heat.shape[0])
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat.to_numpy(dtype=float), aspect="auto", cmap="coolwarm")
    ax.set_xticks(range(heat.shape[1]))
    ax.set_xticklabels(heat.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(heat.shape[0]))
    ax.set_yticklabels(heat.index, fontsize=9)
    ax.set_title("Immune minus malignant pseudobulk signature median", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_report(
    path: Path,
    queue: pd.DataFrame,
    schema_df: pd.DataFrame,
    group_counts: pd.DataFrame,
    pairwise_df: pd.DataFrame,
) -> None:
    completed = queue[queue["is_complete"].astype(bool)].copy()
    cancer_coverage = (
        sorted(group_counts["cancer_code"].dropna().unique().tolist()) if not group_counts.empty else []
    )
    lines = [
        "# BIO-05 Public scRNA Source Mapping v0.7",
        "",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        f"- Platform: {platform.platform()}",
        f"- Completed scRNA files audited: {len(completed)}",
        f"- Cancer codes with analyzable pseudobulk groups so far: {', '.join(cancer_coverage) if cancer_coverage else 'none'}",
        "",
        "## Completed Files",
        "",
    ]
    for row in completed.itertuples():
        lines.append(
            f"- `{row.source_id}`: `{row.title}` -> `{row.remote_path}` ({row.observed_size_bytes} bytes)"
        )

    lines.extend(["", "## Schema Audit", ""])
    if schema_df.empty:
        lines.append("- No complete file passed schema audit yet.")
    else:
        for source_id, sub in schema_df.groupby("source_id", dropna=False):
            lines.append(f"- `{source_id}`")
            for item in sub.itertuples():
                lines.append(
                    f"  {item.field_role}: `{item.field_name}`; nunique={item.nunique}; top={item.top_values[:160]}"
                )

    lines.extend(["", "## Pairwise Contrast Snapshot", ""])
    if pairwise_df.empty:
        lines.append("- No pairwise donor-aware contrast is available yet.")
    else:
        top = pairwise_df.sort_values(["signature_label", "contrast", "paired_donor_n"], ascending=[True, True, False])
        for row in top.head(24).itertuples():
            lines.append(
                f"- `{row.cancer_code}` `{row.signature_label}` `{row.contrast}`: n={row.paired_donor_n}, median_delta={row.median_delta:.3f}, positive_fraction={row.positive_fraction:.3f}"
            )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = PROJECT_ROOT
    out_dir = root / args.out_dir
    figure_dir = root / args.figure_dir
    report_path = root / args.report
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    queue = load_queue(Path(args.queue))
    queue = queue[queue["source_type"].astype(str).eq("cellxgene")].copy()
    queue["observed_size_bytes"] = [
        Path(p).stat().st_size if Path(p).exists() else 0 for p in queue["remote_path"]
    ]
    queue["is_complete"] = [
        file_is_complete(Path(p), safe_str(expected))
        for p, expected in zip(queue["remote_path"], queue["expected_size"])
    ]
    queue_out = out_dir / f"scrna_queue_audit_{VERSION}.csv"
    queue.assign(remote_path=queue["remote_path"].astype(str)).to_csv(queue_out, index=False)

    signature_genes = load_signature_genes(root / args.signature_table)

    schema_rows: list[dict[str, Any]] = []
    group_frames: list[pd.DataFrame] = []
    score_frames: list[pd.DataFrame] = []
    pair_frames: list[pd.DataFrame] = []

    for row in queue[queue["is_complete"].astype(bool)].itertuples(index=False):
        row_series = pd.Series(row._asdict())
        group_counts, score_long, pairwise_df, schema_info = process_scrna_dataset(
            row_series,
            signature_genes=signature_genes,
            min_group_cells=args.min_group_cells,
            min_signature_genes=args.min_signature_genes,
        )
        schema_rows.extend(schema_info)
        if not group_counts.empty:
            group_frames.append(group_counts)
        if not score_long.empty:
            score_frames.append(score_long)
        if not pairwise_df.empty:
            pair_frames.append(pairwise_df)

    schema_df = pd.DataFrame(schema_rows)
    group_df = pd.concat(group_frames, ignore_index=True) if group_frames else pd.DataFrame()
    score_df = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    pair_df = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()

    summary_df = pd.DataFrame()
    if not score_df.empty:
        summary_df = (
            score_df.groupby(["source_id", "dataset_title", "cancer_code", "signature", "signature_label", "compartment"], dropna=False)[
                "signature_score"
            ]
            .median()
            .reset_index()
            .pivot_table(
                index=["source_id", "dataset_title", "cancer_code", "signature", "signature_label"],
                columns="compartment",
                values="signature_score",
                aggfunc="first",
            )
            .reset_index()
            .rename_axis(None, axis=1)
        )
        for col in ["malignant", "immune", "stromal"]:
            if col not in summary_df.columns:
                summary_df[col] = np.nan
        summary_df = summary_df.rename(
            columns={
                "malignant": "median_malignant",
                "immune": "median_immune",
                "stromal": "median_stromal",
            }
        )

    schema_df.to_csv(out_dir / f"scrna_schema_audit_{VERSION}.csv", index=False)
    group_df.to_csv(out_dir / f"scrna_group_counts_{VERSION}.csv", index=False)
    score_df.to_csv(out_dir / f"scrna_signature_scores_{VERSION}.csv", index=False)
    pair_df.to_csv(out_dir / f"scrna_pairwise_contrasts_{VERSION}.csv", index=False)
    summary_df.to_csv(out_dir / f"scrna_dataset_cancer_signature_summary_{VERSION}.csv", index=False)

    plot_compartment_heatmap(summary_df, figure_dir / f"fig_scrna_immune_minus_malignant_heatmap_{VERSION}.png")
    write_report(report_path, queue, schema_df, group_df, pair_df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
