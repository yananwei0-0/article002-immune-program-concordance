#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

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


import argparse
import importlib.util
import platform
import tarfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse


PROJECT_ROOT = Path(
    _workspace_path('${DATA_WORKSPACE}/14_pan_cancer_proteogenomic_immune_evasion')
)
VERSION = "v0_8"
MIN_GROUP_CELLS = 20
MIN_SIGNATURE_GENES = 3


@dataclass(frozen=True)
class DatasetSpec:
    source_id: str
    dataset_title: str
    data_tar: Path
    cancer_code: str
    sample_filter_col: str
    sample_filter_value: str
    disease_filter_value: str | None


DATASETS = [
    DatasetSpec(
        source_id="3ca-peng2019-pdac",
        dataset_title="Peng 2019 Pancreas (3CA official PDAC route)",
        data_tar=Path(expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/downloads/BIO-05/3ca/Data_Peng2019_Pancreas.tar.gz")),
        cancer_code="pdac",
        sample_filter_col="cancer_type",
        sample_filter_value="Pancreatic Ductal Adenocarcinoma",
        disease_filter_value=None,
    ),
    DatasetSpec(
        source_id="3ca-regner2021-ucec",
        dataset_title="Regner 2021 Ovarian page endometrial subset (3CA official UCEC route)",
        data_tar=Path(expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/downloads/BIO-05/3ca/Data_Regner2021_Ovarian.tar.gz")),
        cancer_code="ucec",
        sample_filter_col="cancer_type",
        sample_filter_value="Endometrial Cancer",
        disease_filter_value="Endometrial cancer",
    ),
]


COMPARTMENT_MAP = {
    "pdac": {
        "malignant": {"malignant"},
        "immune": {"macrophage", "t_cell", "b_cell", "mast"},
        "stromal": {"fibroblast", "stellate", "endothelial"},
    },
    "ucec": {
        "malignant": {"malignant"},
        "immune": {"macrophage", "t_cell", "b_cell", "mast"},
        "stromal": {"fibroblast", "endothelial", "myocyte"},
    },
}


def load_v07_module():
    script_path = PROJECT_ROOT / "scripts/run_bio05_scrna_source_mapping_v0_7.py"
    spec = importlib.util.spec_from_file_location("bio05_scrna_v07", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load base script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V07 = load_v07_module()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="BIO-05 3CA source mapping for PDAC and UCEC.")
    parser.add_argument(
        "--signature-table",
        default="04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="04_processed/discovery_v0_8_public_validation/threeca_incremental",
    )
    parser.add_argument(
        "--report",
        default="05_reports/public_3ca_source_mapping_v0_8_report.md",
    )
    parser.add_argument("--min-group-cells", type=int, default=MIN_GROUP_CELLS)
    parser.add_argument("--min-signature-genes", type=int, default=MIN_SIGNATURE_GENES)
    return parser.parse_args()


def read_tar_csv(tar_path: Path, suffix: str) -> pd.DataFrame:
    with tarfile.open(tar_path, "r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith(suffix))
        with tf.extractfile(member) as fh:
            return pd.read_csv(fh)


def read_tar_genes(tar_path: Path) -> pd.Index:
    with tarfile.open(tar_path, "r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith("/Genes.txt"))
        with tf.extractfile(member) as fh:
            genes = [line.decode("utf-8", "ignore").strip() for line in fh if line.strip()]
    return pd.Index([g.upper() for g in genes], dtype="object")


def read_tar_matrix(tar_path: Path) -> sparse.csr_matrix:
    with tarfile.open(tar_path, "r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith("/Exp_data_UMIcounts.mtx"))
        with tf.extractfile(member) as fh:
            matrix = spio.mmread(fh)
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(np.asarray(matrix))
    else:
        matrix = matrix.tocsr()
    return matrix


def classify_compartment(cell_type: Any, cancer_code: str) -> tuple[str, str]:
    cell_type_norm = V07.normalize_text(cell_type)
    mapping = COMPARTMENT_MAP[cancer_code]
    if cell_type_norm in mapping["malignant"]:
        return "malignant", "explicit_cell_type_malignant"
    if cell_type_norm in mapping["immune"]:
        return "immune", "explicit_cell_type_immune"
    if cell_type_norm in mapping["stromal"]:
        return "stromal", "explicit_cell_type_stromal"
    return "other", "unmapped"


def summarize_schema(spec: DatasetSpec, samples: pd.DataFrame, cells: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for role, field_name, frame in [
        ("sample", "sample", samples),
        ("patient", "patient", samples if "patient" in samples.columns else pd.DataFrame()),
        ("cell_type", "cell_type", cells),
        ("disease", "disease", cells if "disease" in cells.columns else pd.DataFrame()),
        ("histology", "histology", cells if "histology" in cells.columns else pd.DataFrame()),
    ]:
        if field_name not in frame.columns:
            rows.append(
                {
                    "source_id": spec.source_id,
                    "dataset_title": spec.dataset_title,
                    "field_role": role,
                    "field_name": "",
                    "nunique": 0,
                    "top_values": "",
                }
            )
            continue
        values = frame[field_name].astype(str)
        top_values = values.value_counts().head(10)
        rows.append(
            {
                "source_id": spec.source_id,
                "dataset_title": spec.dataset_title,
                "field_role": role,
                "field_name": field_name,
                "nunique": int(values.nunique()),
                "top_values": "; ".join(f"{idx}:{int(val)}" for idx, val in top_values.items()),
            }
        )
    return rows


def process_dataset(
    spec: DatasetSpec,
    signature_genes: dict[str, list[str]],
    min_group_cells: int,
    min_signature_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    samples = read_tar_csv(spec.data_tar, "/Samples.csv")
    cells = read_tar_csv(spec.data_tar, "/Cells.csv")
    genes = read_tar_genes(spec.data_tar)
    matrix = read_tar_matrix(spec.data_tar)
    if matrix.shape[0] != len(genes):
        raise ValueError(f"{spec.source_id}: genes mismatch {matrix.shape[0]} vs {len(genes)}")
    if matrix.shape[1] != len(cells):
        raise ValueError(f"{spec.source_id}: cells mismatch {matrix.shape[1]} vs {len(cells)}")

    keep_samples = samples.loc[
        samples[spec.sample_filter_col].astype(str).eq(spec.sample_filter_value),
        "sample",
    ].astype(str)
    sample_set = set(keep_samples.tolist())

    obs = cells.copy()
    obs["sample"] = obs["sample"].astype(str)
    obs = obs[obs["sample"].isin(sample_set)].copy()
    if spec.disease_filter_value and "disease" in obs.columns:
        obs = obs[obs["disease"].astype(str).eq(spec.disease_filter_value)].copy()

    sample_to_patient = (
        samples[samples["sample"].astype(str).isin(sample_set)][["sample", "patient"]]
        .astype(str)
        .drop_duplicates()
        if "patient" in samples.columns
        else pd.DataFrame(columns=["sample", "patient"])
    )
    if not sample_to_patient.empty:
        obs = obs.merge(sample_to_patient, on="sample", how="left")
        obs["donor_id_resolved"] = obs["patient"].astype(str)
    else:
        obs["donor_id_resolved"] = obs["sample"].astype(str)

    comp = [classify_compartment(x, spec.cancer_code) for x in obs["cell_type"]]
    obs["broad_compartment"] = [x[0] for x in comp]
    obs["compartment_rule"] = [x[1] for x in comp]
    obs["cancer_code"] = spec.cancer_code

    schema_rows = summarize_schema(spec, samples[samples["sample"].astype(str).isin(sample_set)].copy(), obs)
    obs = obs[
        obs["broad_compartment"].isin(["malignant", "immune", "stromal"])
        & obs["donor_id_resolved"].astype(str).ne("")
    ].copy()
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), schema_rows

    group_counts = (
        obs.groupby(["cancer_code", "donor_id_resolved", "broad_compartment", "compartment_rule"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    group_counts["source_id"] = spec.source_id
    group_counts["dataset_title"] = spec.dataset_title

    keep = group_counts[group_counts["n_cells"] >= min_group_cells].copy()
    if keep.empty:
        return group_counts, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), schema_rows

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

    matrix_cells_by_genes = matrix.T.tocsr()
    sub = matrix_cells_by_genes[obs.index.to_numpy(), :]
    groups = obs["group_id"].astype("category")
    pb = V07.build_group_mean_matrix(sub, groups.cat.codes.to_numpy(), len(groups.cat.categories))
    pb = V07.cpm_log1p(pb)
    expr = pd.DataFrame(pb.toarray(), index=groups.cat.categories, columns=genes)
    expr = expr.T.groupby(level=0).mean().T

    score_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for cancer_code, cancer_groups in obs.groupby("cancer_code", dropna=False):
        group_ids = sorted(cancer_groups["group_id"].unique())
        cancer_expr = expr.loc[group_ids].copy()
        cancer_z = V07.zscore_columns(cancer_expr)
        for signature, geneset in signature_genes.items():
            recovered = [g for g in geneset if g in cancer_z.columns]
            if len(recovered) < min_signature_genes:
                continue
            scores = cancer_z[recovered].mean(axis=1)
            score_df = pd.DataFrame(
                {
                    "group_id": scores.index,
                    "signature": signature,
                    "signature_label": V07.SIGNATURE_LABELS.get(signature, signature),
                    "signature_score": scores.to_numpy(dtype=float),
                    "n_genes_used": len(recovered),
                    "source_id": spec.source_id,
                    "dataset_title": spec.dataset_title,
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
                    "source_id": spec.source_id,
                    "dataset_title": spec.dataset_title,
                    "cancer_code": cancer_code,
                    "signature": signature,
                    "signature_label": V07.SIGNATURE_LABELS.get(signature, signature),
                    "n_groups": int(score_df["group_id"].nunique()),
                    "n_donors": int(score_df["donor_id"].nunique()),
                    "n_genes_used": len(recovered),
                    "median_malignant": med.get("malignant"),
                    "median_immune": med.get("immune"),
                    "median_stromal": med.get("stromal"),
                }
            )

    score_long = pd.DataFrame(score_rows)
    summary_df = pd.DataFrame(summary_rows)
    if score_long.empty:
        return group_counts, score_long, pd.DataFrame(), summary_df, schema_rows

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
    return group_counts, score_long, pd.DataFrame(pair_rows), summary_df, schema_rows


def write_report(
    path: Path,
    specs: list[DatasetSpec],
    schema_df: pd.DataFrame,
    group_df: pd.DataFrame,
    pair_df: pd.DataFrame,
) -> None:
    cancer_codes = sorted(group_df["cancer_code"].dropna().unique().tolist()) if not group_df.empty else []
    lines = [
        "# BIO-05 3CA Source Mapping v0.8",
        "",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        f"- Platform: {platform.platform()}",
        f"- Datasets attempted: {len(specs)}",
        f"- Cancer codes with analyzable pseudobulk groups: {', '.join(cancer_codes) if cancer_codes else 'none'}",
        "",
        "## Completed Files",
        "",
    ]
    for spec in specs:
        lines.append(f"- `{spec.source_id}`: `{spec.dataset_title}` -> `{spec.data_tar}`")
    lines.extend(["", "## Schema Audit", ""])
    for source_id, sub in schema_df.groupby("source_id", dropna=False):
        lines.append(f"- `{source_id}`")
        for item in sub.itertuples():
            lines.append(
                f"  {item.field_role}: `{item.field_name}`; nunique={item.nunique}; top={item.top_values[:160]}"
            )
    lines.extend(["", "## Pairwise Contrast Snapshot", ""])
    if pair_df.empty:
        lines.append("- No pairwise donor-aware contrast was available.")
    else:
        top = pair_df.sort_values(["signature_label", "contrast", "paired_donor_n"], ascending=[True, True, False])
        for row in top.head(24).itertuples():
            lines.append(
                f"- `{row.cancer_code}` `{row.signature_label}` `{row.contrast}`: n={row.paired_donor_n}, median_delta={row.median_delta:.3f}, positive_fraction={row.positive_fraction:.3f}"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = PROJECT_ROOT
    out_dir = root / args.out_dir
    report_path = root / args.report
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    signature_genes = V07.load_signature_genes(root / args.signature_table)

    schema_rows: list[dict[str, Any]] = []
    group_frames: list[pd.DataFrame] = []
    score_frames: list[pd.DataFrame] = []
    pair_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []

    for spec in DATASETS:
        group_counts, score_df, pair_df, summary_df, schema_info = process_dataset(
            spec,
            signature_genes=signature_genes,
            min_group_cells=args.min_group_cells,
            min_signature_genes=args.min_signature_genes,
        )
        schema_rows.extend(schema_info)
        if not group_counts.empty:
            group_frames.append(group_counts)
        if not score_df.empty:
            score_frames.append(score_df)
        if not pair_df.empty:
            pair_frames.append(pair_df)
        if not summary_df.empty:
            summary_frames.append(summary_df)

    schema_df = pd.DataFrame(schema_rows)
    group_df = pd.concat(group_frames, ignore_index=True) if group_frames else pd.DataFrame()
    score_df = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    pair_df = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    summary_df = pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame()

    schema_df.to_csv(out_dir / f"threeca_schema_audit_{VERSION}.csv", index=False)
    group_df.to_csv(out_dir / f"threeca_group_counts_{VERSION}.csv", index=False)
    score_df.to_csv(out_dir / f"threeca_signature_scores_{VERSION}.csv", index=False)
    pair_df.to_csv(out_dir / f"threeca_pairwise_contrasts_{VERSION}.csv", index=False)
    summary_df.to_csv(out_dir / f"threeca_dataset_cancer_signature_summary_{VERSION}.csv", index=False)
    write_report(report_path, DATASETS, schema_df, group_df, pair_df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
