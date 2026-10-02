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


import argparse
import csv
import gzip
import json
import platform
import re
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr

from run_pan_cancer_proteomic_immune_evasion_discovery_v0_1 import (
    SIGNATURES,
    ensg_stable,
    safe_str,
    zscore_rows,
)
from run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap import (
    SIGNATURE_ANCHORS,
    SIGNATURE_LABELS,
    bh_fdr,
    depmap_symbol_entrez_map,
    read_selected_tcga_rows,
)


PROJECT_ROOT = Path(
    _workspace_os.environ.get("ARTICLE002_DATA_ROOT", str(Path.cwd()))
).expanduser().resolve()
VERSION = "v0_5"
MIN_CORRELATION_N = 25
MIN_SUBTYPE_N = 20

TCGA_CANCER_MAP = {
    "brca": "BRCA",
    "ccrcc": "KIRC",
    "coad": "COAD",
    "gbm": "GBM",
    "hnscc": "HNSC",
    "lscc": "LUSC",
    "luad": "LUAD",
    "ov": "OV",
    "pdac": "PAAD",
    "ucec": "UCEC",
}
TCGA_CANCER_ORDER = [TCGA_CANCER_MAP[k] for k in TCGA_CANCER_MAP]
SUBTYPE_ORDER = [
    "Wound Healing (Immune C1)",
    "IFN-gamma Dominant (Immune C2)",
    "Inflammatory (Immune C3)",
    "Lymphocyte Depleted (Immune C4)",
    "Immunologically Quiet (Immune C5)",
    "TGF-beta Dominant (Immune C6)",
]
LOCKED_CLAIM_TIERS = {
    signature: (
        "locked_bounded_comparator_program"
        if signature == "proteasome_antigen_processing"
        else "locked_immune_ecology_program"
    )
    for signature in SIGNATURE_LABELS
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="v0.5 TCGA PanCanAtlas sensitivity for matched RNA-concordant proteogenomic immune ecology."
    )
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--harmonization-dir", default="04_processed/harmonization_v0_1")
    parser.add_argument("--discovery-v0-4-dir", default="04_processed/discovery_v0_4")
    parser.add_argument("--public-manifest", default="04_processed/public_pan_cancer_download_manifest_v0_1.csv")
    parser.add_argument("--out-dir", default="04_processed/discovery_v0_5_tcga_sensitivity")
    parser.add_argument("--figure-dir", default="06_figures/discovery_v0_5_tcga_sensitivity")
    parser.add_argument("--report", default="05_reports/tcga_xena_sensitivity_v0_5_report.md")
    parser.add_argument("--summary", default="05_reports/tcga_xena_sensitivity_v0_5_summary.json")
    parser.add_argument("--min-correlation-n", type=int, default=MIN_CORRELATION_N)
    parser.add_argument("--min-subtype-n", type=int, default=MIN_SUBTYPE_N)
    return parser.parse_args()


def rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def load_tcga_meta(harmonization_dir: Path) -> pd.DataFrame:
    usecols = [
        "layer",
        "sample_id",
        "patient_id",
        "sample_type_code",
        "sample_type_label",
        "cancer_hint",
        "primary_disease_hint",
    ]
    df = pd.read_csv(harmonization_dir / "tcga_sample_map_v0_1.csv", usecols=usecols)
    df = df[
        (df["layer"].astype(str) == "rna_expression")
        & (pd.to_numeric(df["sample_type_code"], errors="coerce") == 1)
        & (df["cancer_hint"].astype(str).isin(TCGA_CANCER_ORDER))
    ].copy()
    df = df.sort_values(["patient_id", "sample_id"]).drop_duplicates("patient_id", keep="first")
    return df.reset_index(drop=True)


def build_signature_row_maps(depmap_feature_map: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    symbol_to_entrez, _entrez_to_symbol = depmap_symbol_entrez_map(depmap_feature_map)
    signature_row_ids: dict[str, set[str]] = defaultdict(set)
    row_id_to_symbol: dict[str, set[str]] = defaultdict(set)
    for signature, spec in SIGNATURES.items():
        for gene in spec["genes"]:
            signature_row_ids[signature].add(gene)
            row_id_to_symbol[gene].add(gene)
            for entrez in symbol_to_entrez.get(gene, set()):
                signature_row_ids[signature].add(entrez)
                row_id_to_symbol[entrez].add(gene)
    return signature_row_ids, row_id_to_symbol


def compute_signature_scores(
    root: Path,
    harmonization_dir: Path,
    out_dir: Path,
    sample_meta: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    depmap_feature_map = pd.read_csv(harmonization_dir / "depmap_gene_feature_map_v0_1.csv")
    signature_row_ids, row_id_to_symbol = build_signature_row_maps(depmap_feature_map)
    sample_ids = sample_meta["sample_id"].astype(str).tolist()
    all_targets = set()
    for row_ids in signature_row_ids.values():
        all_targets.update(row_ids)
    expr_path = root / "03_downloads/tcga_xena_pancanatlas/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    expr, kept_samples = read_selected_tcga_rows(expr_path, all_targets, sample_ids)
    sample_meta = sample_meta[sample_meta["sample_id"].isin(kept_samples)].copy()
    sample_meta = sample_meta.set_index("sample_id").loc[kept_samples].reset_index()
    expr = expr.loc[:, kept_samples]

    recovery_rows: list[dict[str, Any]] = []
    score_rows: list[dict[str, Any]] = []
    for signature, spec in SIGNATURES.items():
        row_ids = sorted(signature_row_ids.get(signature, set()).intersection(set(expr.index.astype(str))))
        mapped_symbols = sorted(
            {
                symbol
                for row_id in row_ids
                for symbol in row_id_to_symbol.get(row_id, set())
                if symbol in spec["genes"]
            }
        )
        recovery_rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "requested_gene_count": len(spec["genes"]),
                "mapped_expression_row_count": len(row_ids),
                "mapped_gene_symbol_count": len(mapped_symbols),
                "mapped_gene_symbols": ";".join(mapped_symbols),
                "mapping_source": "TCGA RNA row IDs matched by gene symbol first and DepMap-derived Entrez IDs as fallback",
            }
        )
        if len(mapped_symbols) < 3 or not row_ids:
            continue
        for cancer, cancer_meta in sample_meta.groupby("cancer_hint", sort=False):
            cancer_samples = cancer_meta["sample_id"].astype(str).tolist()
            cancer_expr = expr.loc[row_ids, cancer_samples]
            cancer_z = pd.DataFrame(
                zscore_rows(cancer_expr.to_numpy(dtype=float)),
                index=cancer_expr.index,
                columns=cancer_expr.columns,
            )
            scores = cancer_z.mean(axis=0, skipna=True)
            for sample_id, value in scores.items():
                if not np.isfinite(value):
                    continue
                score_rows.append(
                    {
                        "sample_id": sample_id,
                        "signature": signature,
                        "signature_label": SIGNATURE_LABELS.get(signature, signature),
                        "rna_signature_score": float(value),
                        "n_requested_genes": len(spec["genes"]),
                        "n_mapped_genes": len(mapped_symbols),
                    }
                )
    scores_long = pd.DataFrame(score_rows).merge(sample_meta, on="sample_id", how="left", validate="many_to_one")
    recovery_df = pd.DataFrame(recovery_rows)
    scores_long.to_csv(out_dir / "tcga_xena_signature_scores_v0_5.csv", index=False)
    recovery_df.to_csv(out_dir / "tcga_xena_signature_recovery_v0_5.csv", index=False)
    return scores_long, recovery_df


def load_immune_signature_matrix(root: Path, sample_ids: list[str]) -> pd.DataFrame:
    immune_path = root / "03_downloads/tcga_xena_pancanatlas/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz"
    with open_text(immune_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        immune_names = []
        for row in reader:
            if row:
                immune_names.append(safe_str(row[0]))
    data, kept_samples = read_selected_tcga_rows(immune_path, set(immune_names), sample_ids)
    if data.empty:
        return pd.DataFrame({"sample_id": kept_samples})
    out = data.T.reset_index().rename(columns={"index": "sample_id"})
    out.columns = [safe_str(c) for c in out.columns]
    return out


def load_subtype(root: Path, sample_ids: list[str]) -> pd.DataFrame:
    path = root / "03_downloads/tcga_xena_pancanatlas/Subtype_Immune_Model_Based.txt.gz"
    subtype = pd.read_csv(path, sep="\t", compression="gzip")
    subtype = subtype.rename(columns={"Subtype_Immune_Model_Based": "immune_subtype"})
    subtype["sample_id"] = subtype["sample"].astype(str)
    subtype = subtype[subtype["sample_id"].isin(sample_ids)].copy()
    subtype["immune_subtype"] = subtype["immune_subtype"].astype(str)
    subtype["immune_subtype_short"] = subtype["immune_subtype"].str.extract(r"(Immune C\d)")
    return subtype[["sample_id", "immune_subtype", "immune_subtype_short"]]


def summarize_68immune_correlations(
    scores_long: pd.DataFrame,
    immune_df: pd.DataFrame,
    out_dir: Path,
    min_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = scores_long.merge(immune_df, on="sample_id", how="left", validate="many_to_one")
    immune_cols = [c for c in immune_df.columns if c != "sample_id"]
    corr_rows: list[dict[str, Any]] = []
    for (signature, signature_label, cancer), group in merged.groupby(
        ["signature", "signature_label", "cancer_hint"], dropna=False
    ):
        if not safe_str(cancer):
            continue
        for immune_sig in immune_cols:
            sub = group[["rna_signature_score", immune_sig]].dropna()
            if len(sub) < min_n:
                continue
            if sub["rna_signature_score"].nunique() < 3 or sub[immune_sig].nunique() < 3:
                continue
            rho, p = spearmanr(sub["rna_signature_score"], sub[immune_sig])
            corr_rows.append(
                {
                    "signature": signature,
                    "signature_label": signature_label,
                    "cancer": cancer,
                    "immune_signature_68": immune_sig,
                    "n": int(len(sub)),
                    "spearman_rho": float(rho),
                    "p_value": float(p),
                    "p_family": "tcga_xena_68immune_within_cancer_correlations_v0_5",
                    "interpretation_boundary": "Independent TCGA RNA sensitivity only; not matched CPTAC validation and not malignant-cell attribution.",
                }
            )
    corr_df = pd.DataFrame(corr_rows)
    if not corr_df.empty:
        corr_df["fdr_bh"] = bh_fdr(pd.to_numeric(corr_df["p_value"], errors="coerce").to_numpy())
    corr_df.to_csv(out_dir / "tcga_xena_68immune_within_cancer_correlations_v0_5.csv", index=False)

    meta_rows: list[dict[str, Any]] = []
    for (signature, signature_label, immune_sig), group in corr_df.groupby(
        ["signature", "signature_label", "immune_signature_68"], dropna=False
    ):
        meta_rows.append(
            {
                "signature": signature,
                "signature_label": signature_label,
                "immune_signature_68": immune_sig,
                "cancer_strata": int(group["cancer"].nunique()),
                "total_n": int(group["n"].sum()),
                "median_within_cancer_rho": float(group["spearman_rho"].median()),
                "positive_cancer_fraction": float((group["spearman_rho"] > 0).mean()),
                "max_abs_within_cancer_rho": float(group["spearman_rho"].abs().max()),
                "significant_fdr_fraction": float((pd.to_numeric(group["fdr_bh"], errors="coerce") < 0.05).mean()),
            }
        )
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(out_dir / "tcga_xena_68immune_meta_summary_v0_5.csv", index=False)
    return corr_df, meta_df


def summarize_subtypes(
    scores_long: pd.DataFrame,
    subtype_df: pd.DataFrame,
    out_dir: Path,
    min_subtype_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = scores_long.merge(subtype_df, on="sample_id", how="left", validate="many_to_one")
    merged = merged.dropna(subset=["immune_subtype"]).copy()

    subtype_rows: list[dict[str, Any]] = []
    test_rows: list[dict[str, Any]] = []
    for (signature, signature_label), group in merged.groupby(["signature", "signature_label"], dropna=False):
        subtype_vectors = []
        subtype_medians = {}
        subtype_ns = {}
        for subtype in SUBTYPE_ORDER:
            values = pd.to_numeric(
                group.loc[group["immune_subtype"] == subtype, "rna_signature_score"], errors="coerce"
            ).dropna()
            subtype_ns[subtype] = int(len(values))
            subtype_medians[subtype] = float(values.median()) if len(values) else np.nan
            if len(values) >= min_subtype_n:
                subtype_vectors.append(values.to_numpy())
            subtype_rows.append(
                {
                    "signature": signature,
                    "signature_label": signature_label,
                    "immune_subtype": subtype,
                    "immune_subtype_short": re.search(r"(Immune C\d)", subtype).group(1) if re.search(r"(Immune C\d)", subtype) else "",
                    "n": int(len(values)),
                    "meets_min_n_for_ordering": int(len(values) >= min_subtype_n),
                    "median_score": float(values.median()) if len(values) else np.nan,
                    "mean_score": float(values.mean()) if len(values) else np.nan,
                }
            )
        valid_medians = {
            k: v
            for k, v in subtype_medians.items()
            if np.isfinite(v) and subtype_ns.get(k, 0) >= min_subtype_n
        }
        if not valid_medians:
            valid_medians = {k: v for k, v in subtype_medians.items() if np.isfinite(v)}
        if len(subtype_vectors) >= 2:
            kw_stat, kw_p = kruskal(*subtype_vectors)
        else:
            kw_stat, kw_p = np.nan, np.nan
        top_subtype = max(valid_medians, key=valid_medians.get) if valid_medians else ""
        bottom_subtype = min(valid_medians, key=valid_medians.get) if valid_medians else ""
        top_val = valid_medians.get(top_subtype, np.nan)
        bottom_val = valid_medians.get(bottom_subtype, np.nan)
        test_rows.append(
            {
                "signature": signature,
                "signature_label": signature_label,
                "subtype_groups_with_n_ge_min": int(sum(n >= min_subtype_n for n in subtype_ns.values())),
                "kruskal_statistic": float(kw_stat) if np.isfinite(kw_stat) else np.nan,
                "kruskal_p_value": float(kw_p) if np.isfinite(kw_p) else np.nan,
                "top_median_subtype": top_subtype,
                "top_median_subtype_short": re.search(r"(Immune C\d)", top_subtype).group(1) if re.search(r"(Immune C\d)", top_subtype) else "",
                "bottom_median_subtype": bottom_subtype,
                "bottom_median_subtype_short": re.search(r"(Immune C\d)", bottom_subtype).group(1) if re.search(r"(Immune C\d)", bottom_subtype) else "",
                "top_minus_bottom_median_score": float(top_val - bottom_val) if np.isfinite(top_val) and np.isfinite(bottom_val) else np.nan,
            }
        )
    subtype_df_out = pd.DataFrame(subtype_rows)
    test_df = pd.DataFrame(test_rows)
    if not test_df.empty:
        test_df["kruskal_fdr_bh"] = bh_fdr(pd.to_numeric(test_df["kruskal_p_value"], errors="coerce").to_numpy())
    subtype_df_out.to_csv(out_dir / "tcga_xena_immune_subtype_summary_v0_5.csv", index=False)
    test_df.to_csv(out_dir / "tcga_xena_immune_subtype_tests_v0_5.csv", index=False)
    return subtype_df_out, test_df


def build_signature_summary(
    discovery_v0_4_dir: Path,
    meta_df: pd.DataFrame,
    subtype_tests: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    del discovery_v0_4_dir  # retained in the public CLI for backward compatibility
    rows: list[dict[str, Any]] = []
    for signature, label in SIGNATURE_LABELS.items():
        meta_sig = meta_df[meta_df["signature"] == signature].copy()
        subtype_row = subtype_tests[subtype_tests["signature"] == signature]
        anchors = SIGNATURE_ANCHORS.get(signature, [])
        anchor_meta = meta_sig[meta_sig["immune_signature_68"].isin(anchors)].copy()
        top_hits = meta_sig.sort_values(["median_within_cancer_rho", "positive_cancer_fraction"], ascending=False).head(3)
        row: dict[str, Any] = {
            "signature": signature,
            "signature_label": label,
            "v0_4_claim_tier": LOCKED_CLAIM_TIERS[signature],
            "selected_anchor_count": int(len(anchor_meta)),
            "selected_anchor_median_rho": float(anchor_meta["median_within_cancer_rho"].median()) if len(anchor_meta) else np.nan,
            "selected_anchor_min_positive_cancer_fraction": float(anchor_meta["positive_cancer_fraction"].min()) if len(anchor_meta) else np.nan,
            "selected_anchor_names": ";".join(anchor_meta["immune_signature_68"].astype(str).tolist()),
            "top_median_68immune_hit_1": "",
            "top_median_68immune_rho_1": np.nan,
            "top_median_68immune_hit_2": "",
            "top_median_68immune_rho_2": np.nan,
            "top_median_68immune_hit_3": "",
            "top_median_68immune_rho_3": np.nan,
            "top_median_subtype": "",
            "bottom_median_subtype": "",
            "top_minus_bottom_median_score": np.nan,
            "subtype_kruskal_fdr_bh": np.nan,
        }
        for i, (_, hit) in enumerate(top_hits.iterrows(), start=1):
            row[f"top_median_68immune_hit_{i}"] = safe_str(hit["immune_signature_68"])
            row[f"top_median_68immune_rho_{i}"] = float(hit["median_within_cancer_rho"])
        if len(subtype_row):
            sr = subtype_row.iloc[0]
            row["top_median_subtype"] = safe_str(sr["top_median_subtype_short"]) or safe_str(sr["top_median_subtype"])
            row["bottom_median_subtype"] = safe_str(sr["bottom_median_subtype_short"]) or safe_str(sr["bottom_median_subtype"])
            row["top_minus_bottom_median_score"] = float(sr["top_minus_bottom_median_score"]) if pd.notna(sr["top_minus_bottom_median_score"]) else np.nan
            row["subtype_kruskal_fdr_bh"] = float(sr["kruskal_fdr_bh"]) if pd.notna(sr["kruskal_fdr_bh"]) else np.nan
        rows.append(row)
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(out_dir / "tcga_xena_signature_summary_v0_5.csv", index=False)
    return summary_df


def write_figures(
    meta_df: pd.DataFrame,
    subtype_summary: pd.DataFrame,
    fig_dir: Path,
) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: list[str] = []

    if not meta_df.empty:
        heat = meta_df.pivot(index="signature_label", columns="immune_signature_68", values="median_within_cancer_rho")
        heat = heat.reindex(index=[SIGNATURE_LABELS[s] for s in SIGNATURE_LABELS])
        heat = heat.reindex(columns=sorted(heat.columns), fill_value=np.nan)
        fig, ax = plt.subplots(figsize=(24, 4.6))
        im = ax.imshow(heat.fillna(0).to_numpy(), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
        ax.set_yticks(range(heat.shape[0]))
        ax.set_yticklabels(heat.index, fontsize=8)
        ax.set_xticks(range(heat.shape[1]))
        ax.set_xticklabels(heat.columns, rotation=90, fontsize=6)
        ax.set_title("TCGA PanCanAtlas 68 immune-signature sensitivity across 10 overlapping cancers", fontsize=11)
        cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
        cbar.set_label("Median within-cancer Spearman rho", fontsize=8)
        fig.tight_layout()
        out = fig_dir / "fig15_tcga_xena_68immune_sensitivity_heatmap_v0_5.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        figure_paths.append(str(out))

    if not subtype_summary.empty:
        subtype_pivot = subtype_summary.pivot(
            index="signature_label", columns="immune_subtype_short", values="median_score"
        )
        subtype_order_short = [re.search(r"(Immune C\d)", s).group(1) for s in SUBTYPE_ORDER]
        subtype_pivot = subtype_pivot.reindex(index=[SIGNATURE_LABELS[s] for s in SIGNATURE_LABELS])
        subtype_pivot = subtype_pivot.reindex(columns=subtype_order_short)
        fig, ax = plt.subplots(figsize=(8.8, 4.8))
        im = ax.imshow(subtype_pivot.fillna(0).to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
        ax.set_yticks(range(subtype_pivot.shape[0]))
        ax.set_yticklabels(subtype_pivot.index, fontsize=8)
        ax.set_xticks(range(subtype_pivot.shape[1]))
        ax.set_xticklabels(subtype_pivot.columns, fontsize=8)
        ax.set_title("TCGA immune-subtype sensitivity of within-cancer RNA signature scores", fontsize=11)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Median within-cancer-standardized score", fontsize=8)
        fig.tight_layout()
        out = fig_dir / "fig16_tcga_xena_immune_subtype_median_heatmap_v0_5.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        figure_paths.append(str(out))

    return figure_paths


def format_float(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return "NA"


def write_report(
    report_path: Path,
    root: Path,
    sample_meta: pd.DataFrame,
    recovery_df: pd.DataFrame,
    signature_summary: pd.DataFrame,
    subtype_tests: pd.DataFrame,
    public_manifest: pd.DataFrame,
    figure_paths: list[str],
) -> None:
    reused_ids = [
        "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
        "Subtype_Immune_Model_Based.txt.gz",
        "TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz",
    ]
    reused = public_manifest[public_manifest["dataset_id"].isin(reused_ids)].copy()
    reused = reused.sort_values("dataset_id")
    lines = [
        "# TCGA PanCanAtlas sensitivity v0.5",
        "",
        "## Scope",
        "",
        "This v0.5 module adds an independent TCGA RNA sensitivity layer to the existing matched CPTAC v0.4 package. It does not claim external validation, malignant-cell attribution, or immune-evasion mechanism. The purpose is to test whether the matched RNA-concordant CPTAC programs map onto canonical TCGA immune subtype and 68-immune-signature structure.",
        "",
        "## Sample Surface",
        "",
        f"- Overlapping TCGA cancers: `{', '.join(TCGA_CANCER_ORDER)}`",
        f"- Unique primary-tumor TCGA patients after one-sample-per-patient filtering: `{int(sample_meta['patient_id'].nunique())}`",
        f"- TCGA samples entering signature scoring: `{int(sample_meta['sample_id'].nunique())}`",
        f"- Signatures with >=3 mapped genes: `{int((recovery_df['mapped_gene_symbol_count'] >= 3).sum())}` of `{int(len(recovery_df))}`",
        "",
        "## Reused Official Assets",
        "",
    ]
    for _, row in reused.iterrows():
        lines.extend(
            [
                f"- `{safe_str(row['dataset_id'])}`",
                f"  - URL: `{safe_str(row['source_reference'])}`",
                f"  - Local path: `{safe_str(row['local_path'])}`",
                f"  - Size bytes: `{int(row['size_bytes'])}`",
                f"  - SHA256: `{safe_str(row['sha256'])}`",
                f"  - Existing-file timestamp: `{safe_str(row['timestamp'])}`",
            ]
        )
    lines.extend(["", "## Main Results", ""])
    for _, row in signature_summary.iterrows():
        sig = safe_str(row["signature_label"])
        anchor_rho = format_float(row["selected_anchor_median_rho"], 3)
        top1 = safe_str(row["top_median_68immune_hit_1"])
        top1_rho = format_float(row["top_median_68immune_rho_1"], 3)
        top2 = safe_str(row["top_median_68immune_hit_2"])
        top2_rho = format_float(row["top_median_68immune_rho_2"], 3)
        subtype_top = safe_str(row["top_median_subtype"])
        subtype_bottom = safe_str(row["bottom_median_subtype"])
        subtype_gap = format_float(row["top_minus_bottom_median_score"], 3)
        subtype_fdr = format_float(row["subtype_kruskal_fdr_bh"], 3)
        lines.append(
            f"- `{sig}`: canonical-anchor median rho `{anchor_rho}`; top 68-immune hits `{top1}` ({top1_rho}) and `{top2}` ({top2_rho}); subtype span `{subtype_top}` to `{subtype_bottom}` with median-gap `{subtype_gap}` and Kruskal FDR `{subtype_fdr}`."
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "- Write this as independent TCGA RNA sensitivity supporting canonical immune-landscape alignment.",
            "- Do not write this as matched CPTAC external validation, tumor-cell-intrinsic immune evasion, or source-origin completion.",
            "- Proteasome remains a boundary/comparator context because the v0.4 common-essential limitation still applies.",
            "",
            "## Figures",
            "",
        ]
    )
    for path in figure_paths:
        lines.append(f"- `{rel(Path(path), root)}`")
    lines.append("")
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = Path(args.project_root)
    harmonization_dir = root / args.harmonization_dir
    discovery_v0_4_dir = root / args.discovery_v0_4_dir
    out_dir = root / args.out_dir
    fig_dir = root / args.figure_dir
    report_path = root / args.report
    summary_path = root / args.summary

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    started = time.time()
    sample_meta = load_tcga_meta(harmonization_dir)
    scores_long, recovery_df = compute_signature_scores(root, harmonization_dir, out_dir, sample_meta)
    scored_samples = sorted(scores_long["sample_id"].astype(str).unique().tolist())
    immune_df = load_immune_signature_matrix(root, scored_samples)
    subtype_df = load_subtype(root, scored_samples)
    corr_df, meta_df = summarize_68immune_correlations(scores_long, immune_df, out_dir, args.min_correlation_n)
    subtype_summary, subtype_tests = summarize_subtypes(scores_long, subtype_df, out_dir, args.min_subtype_n)
    signature_summary = build_signature_summary(discovery_v0_4_dir, meta_df, subtype_tests, out_dir)
    figure_paths = write_figures(meta_df, subtype_summary, fig_dir)

    public_manifest = pd.read_csv(root / args.public_manifest)
    write_report(
        report_path=report_path,
        root=root,
        sample_meta=sample_meta[sample_meta["sample_id"].isin(scored_samples)].copy(),
        recovery_df=recovery_df,
        signature_summary=signature_summary,
        subtype_tests=subtype_tests,
        public_manifest=public_manifest,
        figure_paths=figure_paths,
    )

    summary = {
        "version": VERSION,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "elapsed_seconds": round(time.time() - started, 2),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "tcga_primary_tumor_patients": int(sample_meta["patient_id"].nunique()),
        "tcga_scored_samples": int(len(scored_samples)),
        "overlapping_tcga_cancers": TCGA_CANCER_ORDER,
        "recovery_rows": int(len(recovery_df)),
        "within_cancer_correlation_rows": int(len(corr_df)),
        "within_cancer_correlation_meta_rows": int(len(meta_df)),
        "immune_subtype_rows": int(len(subtype_summary)),
        "immune_subtype_test_rows": int(len(subtype_tests)),
        "figures": [rel(Path(p), root) for p in figure_paths],
        "outputs": {
            "scores": rel(out_dir / "tcga_xena_signature_scores_v0_5.csv", root),
            "recovery": rel(out_dir / "tcga_xena_signature_recovery_v0_5.csv", root),
            "within_cancer_correlations": rel(out_dir / "tcga_xena_68immune_within_cancer_correlations_v0_5.csv", root),
            "meta_summary": rel(out_dir / "tcga_xena_68immune_meta_summary_v0_5.csv", root),
            "immune_subtype_summary": rel(out_dir / "tcga_xena_immune_subtype_summary_v0_5.csv", root),
            "immune_subtype_tests": rel(out_dir / "tcga_xena_immune_subtype_tests_v0_5.csv", root),
            "signature_summary": rel(out_dir / "tcga_xena_signature_summary_v0_5.csv", root),
            "report": rel(report_path, root),
        },
        "boundary": "Independent TCGA RNA sensitivity only; no external validation, tumor-cell-intrinsic immune-evasion, or source-origin completion claim.",
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
