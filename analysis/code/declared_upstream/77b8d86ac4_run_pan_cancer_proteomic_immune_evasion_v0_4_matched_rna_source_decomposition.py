#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import csv
import gzip
import json
import platform
import re
import sys
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
from scipy.stats import spearmanr

from run_pan_cancer_proteomic_immune_evasion_discovery_v0_1 import (
    SIGNATURES,
    ensg_stable,
    rel,
    resolve_path,
    safe_str,
    zscore_rows,
)
from run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap import (
    DEPMAP_GENE_RE,
    SIGNATURE_LABELS,
)
from run_pan_cancer_proteomic_immune_evasion_v0_3_phosphosite_depmap_mechanism import (
    bh_fdr,
)


PROJECT_ROOT = Path(
    expand_legacy_paths("${SCI27_DATA_ROOT}/14_pan_cancer_proteogenomic_immune_evasion")
)
VERSION = "v0_4"
MIN_CORRELATION_N = 20
MIN_RESIDUAL_N = 20
MIN_RNA_MAPPED_GENES = 3
RANDOM_SEED = 20260624
CONTEXT_VARS = ["estimate_immune_score", "estimate_tumor_purity", "cibersort_t_cell", "cibersort_myeloid"]


GENE_SOURCE_RULES = {
    "lymphoid_t_cell_or_cytotoxic": {
        "CD8A",
        "CD8B",
        "GZMA",
        "GZMB",
        "GZMH",
        "GZMK",
        "PRF1",
        "NKG7",
        "GNLY",
        "CTSW",
        "CST7",
        "CCL5",
        "CXCR3",
        "TRAC",
        "TRBC1",
        "TRBC2",
        "PDCD1",
        "CTLA4",
        "LAG3",
        "HAVCR2",
        "TIGIT",
        "ICOS",
        "TNFRSF9",
    },
    "myeloid_or_professional_apc": {
        "LST1",
        "LYZ",
        "TYROBP",
        "AIF1",
        "FCGR3A",
        "FCGR1A",
        "CSF1R",
        "ITGAM",
        "CD68",
        "CD163",
        "MSR1",
        "C1QA",
        "C1QB",
        "C1QC",
        "IL1B",
        "TNF",
        "CXCL8",
        "S100A8",
        "S100A9",
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-DQA1",
        "HLA-DQB1",
        "HLA-DMA",
        "HLA-DMB",
        "CD74",
        "CIITA",
        "LGMN",
        "CTSS",
        "CD80",
        "CD86",
    },
    "stromal_emt_context": {
        "TGFB1",
        "TGFBR1",
        "TGFBR2",
        "SMAD2",
        "SMAD3",
        "VIM",
        "FN1",
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "ACTA2",
        "TAGLN",
        "SERPINE1",
        "THBS1",
        "ITGA5",
        "ITGB1",
        "LOXL2",
        "ZEB1",
        "SNAI1",
        "SNAI2",
    },
    "antigen_processing_or_interferon_mixed": {
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "TAP1",
        "TAP2",
        "TAPBP",
        "NLRC5",
        "ERAP1",
        "ERAP2",
        "CALR",
        "CANX",
        "PDIA3",
        "IFNG",
        "STAT1",
        "IRF1",
        "IRF9",
        "JAK1",
        "JAK2",
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "GBP1",
        "GBP2",
        "GBP5",
        "IDO1",
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "MX1",
        "OAS1",
        "OAS2",
        "OAS3",
        "CD274",
        "PDCD1LG2",
        "TOX",
        "ENTPD1",
    },
    "proteasome_common_essential_risk": {
        "PSMB8",
        "PSMB9",
        "PSMB10",
        "PSME1",
        "PSME2",
        "PSME3",
        "PSMA1",
        "PSMA2",
        "PSMA3",
        "PSMA4",
        "PSMA5",
        "PSMA6",
        "PSMA7",
        "PSMB1",
        "PSMB2",
        "PSMB3",
        "PSMB4",
        "PSMB5",
        "PSMB6",
        "PSMB7",
        "UBC",
        "UBB",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="v0.4 matched CPTAC BCM RNA-protein-phosphoprotein residual and source-aware context analysis."
    )
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--version", default=VERSION)
    parser.add_argument("--qc-dir", default="04_processed/first_pass_qc")
    parser.add_argument("--harmonization-dir", default="04_processed/harmonization_v0_1")
    parser.add_argument("--discovery-v0-1-dir", default="04_processed/discovery_v0_1")
    parser.add_argument("--discovery-v0-2-dir", default="04_processed/discovery_v0_2")
    parser.add_argument("--discovery-v0-3-dir", default="04_processed/discovery_v0_3")
    parser.add_argument("--rna-manifest", default="04_processed/discovery_v0_4/cptac_bcm_transcriptomics_download_manifest_v0_4.csv")
    parser.add_argument("--out-dir", default="04_processed/discovery_v0_4")
    parser.add_argument("--figure-dir", default="06_figures/discovery_v0_4")
    parser.add_argument("--report", default="05_reports/discovery_v0_4_report.md")
    parser.add_argument("--summary", default="05_reports/discovery_v0_4_summary.json")
    parser.add_argument("--min-correlation-n", type=int, default=MIN_CORRELATION_N)
    parser.add_argument("--min-residual-n", type=int, default=MIN_RESIDUAL_N)
    parser.add_argument("--min-rna-mapped-genes", type=int, default=MIN_RNA_MAPPED_GENES)
    parser.add_argument("--random-seed", type=int, default=RANDOM_SEED)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def read_first_row(path: Path) -> list[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        return [safe_str(x).strip() for x in next(reader)]


def normalize_case_id(sample_id: str) -> str:
    sid = safe_str(sample_id).strip()
    if not sid:
        return ""
    if sid.upper().startswith(("QC", "JHU-QC")) or "QC" in sid.upper():
        return ""
    sid = re.sub(r"(_T|_A|_N)$", "", sid, flags=re.I)
    sid = re.sub(r"-(T|N)$", "", sid, flags=re.I)
    sid = re.sub(r"_D\d+$", "", sid, flags=re.I)
    return sid


def infer_rna_tumor_normal(sample_id: str) -> str:
    sid = safe_str(sample_id).strip().upper()
    if not sid:
        return "unknown"
    if re.search(r"_T$", sid):
        return "tumor"
    if re.search(r"(_A|_N)$", sid):
        return "adjacent_or_normal"
    if re.search(r"-(T)$", sid):
        return "tumor"
    if re.search(r"-(N)$", sid):
        return "adjacent_or_normal"
    return "unknown_tumor_like"


def build_signature_ensembl_map(signature_table: pd.DataFrame) -> dict[str, dict[str, set[str]]]:
    mapping: dict[str, dict[str, set[str]]] = {}
    for signature, group in signature_table.groupby("signature", dropna=False):
        symbol_to_ensg: dict[str, set[str]] = defaultdict(set)
        for _, row in group.iterrows():
            gene = safe_str(row.get("gene_symbol"))
            for token in safe_str(row.get("ensembl_stable_ids_from_cptac", "")).split(";"):
                stable = ensg_stable(token)
                if gene and stable:
                    symbol_to_ensg[gene].add(stable)
        mapping[safe_str(signature)] = symbol_to_ensg
    return mapping


def read_bcm_rna_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    sample_ids = read_first_row(path)
    names = ["gene_id"] + sample_ids
    matrix = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        names=names,
        na_values=["NA", "NaN", "nan", ""],
        keep_default_na=True,
        low_memory=False,
    )
    if matrix.empty:
        return sample_ids, [], np.empty((0, len(sample_ids)))
    gene_ids = matrix["gene_id"].astype(str).tolist()
    numeric = matrix[sample_ids].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    return sample_ids, gene_ids, numeric


def score_rna_file(
    item: pd.Series,
    root: Path,
    signature_table: pd.DataFrame,
    signature_map: dict[str, dict[str, set[str]]],
    min_mapped_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cancer = safe_str(item["cancer"])
    source = safe_str(item.get("source", "bcm"))
    path = resolve_path(safe_str(item["project_local_path"]), root)
    dataset_id = f"{source}-{cancer}-transcriptomics::{path.name}"
    sample_ids, gene_ids, numeric = read_bcm_rna_matrix(path)
    stable_to_rows: dict[str, list[int]] = defaultdict(list)
    for i, raw in enumerate(gene_ids):
        stable = ensg_stable(raw)
        if stable:
            stable_to_rows[stable].append(i)

    sample_info = pd.DataFrame(
        {
            "dataset_id": dataset_id,
            "source": source,
            "cancer": cancer,
            "datatype": "transcriptomics",
            "sample_id": sample_ids,
            "normalized_case_id": [normalize_case_id(x) for x in sample_ids],
            "rna_tumor_normal_hint": [infer_rna_tumor_normal(x) for x in sample_ids],
        }
    )
    sample_info["is_tumor_like_for_matching"] = sample_info["rna_tumor_normal_hint"].isin(["tumor", "unknown_tumor_like"]).astype(int)
    sample_info["is_duplicate_case_within_rna_file"] = sample_info.duplicated(["normalized_case_id"], keep=False).astype(int)

    score_rows: list[dict[str, Any]] = []
    recovery_rows: list[dict[str, Any]] = []
    z = zscore_rows(numeric)
    requested_counts = signature_table.groupby("signature")["gene_symbol"].nunique().to_dict()
    for signature in SIGNATURES:
        symbol_to_ensg = signature_map.get(signature, {})
        row_indices: set[int] = set()
        recovered_symbols = set()
        recovered_ensg = set()
        for gene, ens_set in symbol_to_ensg.items():
            gene_hit = False
            for stable in ens_set:
                for idx in stable_to_rows.get(stable, []):
                    row_indices.add(idx)
                    recovered_ensg.add(stable)
                    gene_hit = True
            if gene_hit:
                recovered_symbols.add(gene)
        requested = int(requested_counts.get(signature, 0))
        recovered_genes = len(recovered_symbols)
        feature_count = len(row_indices)
        recovery_fraction = recovered_genes / requested if requested else np.nan
        scoreable = int(recovered_genes >= min_mapped_genes and feature_count > 0)
        recovery_rows.append(
            {
                "dataset_id": dataset_id,
                "source": source,
                "cancer": cancer,
                "datatype": "transcriptomics",
                "signature": signature,
                "requested_gene_count": requested,
                "recovered_gene_count": recovered_genes,
                "recovered_feature_count": feature_count,
                "recovery_fraction": recovery_fraction,
                "scoreable": scoreable,
                "recovered_gene_symbols": ";".join(sorted(recovered_symbols)),
                "recovered_ensembl_ids": ";".join(sorted(recovered_ensg)),
            }
        )
        if not scoreable:
            continue
        row_order = sorted(row_indices)
        values = np.nanmean(z[row_order, :], axis=0)
        for i, sample_id in enumerate(sample_ids):
            score_rows.append(
                {
                    "dataset_id": dataset_id,
                    "source": source,
                    "cancer": cancer,
                    "datatype": "transcriptomics",
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "sample_id": sample_id,
                    "normalized_case_id": normalize_case_id(sample_id),
                    "rna_tumor_normal_hint": infer_rna_tumor_normal(sample_id),
                    "is_tumor_like_for_matching": int(infer_rna_tumor_normal(sample_id) in {"tumor", "unknown_tumor_like"}),
                    "signature_score": values[i],
                    "n_features_used": feature_count,
                    "n_recovered_genes": recovered_genes,
                    "recovery_fraction": recovery_fraction,
                }
            )
    return pd.DataFrame(score_rows), pd.DataFrame(recovery_rows), sample_info


def prepare_existing_cptac_scores(scores: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    out = out[(out["source"].astype(str) == "bcm") & (out["tumor_normal_hint"].astype(str) != "normal")].copy()
    out = out[out["datatype"].astype(str).isin(["proteomics", "phosphoproteomics"])].copy()
    out["normalized_case_id"] = out["inferred_case_id"].map(normalize_case_id)
    out = out[out["normalized_case_id"].astype(str) != ""].copy()
    return out


def build_matched_table(rna_scores: pd.DataFrame, cptac_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rna = rna_scores[rna_scores["is_tumor_like_for_matching"].astype(int).eq(1)].copy()
    rna = rna[rna["normalized_case_id"].astype(str) != ""].copy()
    rows = []
    flow_rows = []
    for (cancer, signature, datatype), layer in cptac_scores.groupby(["cancer", "signature", "datatype"], dropna=False):
        layer_small = layer[
            [
                "source",
                "cancer",
                "datatype",
                "tumor_normal_hint",
                "dataset_id",
                "signature",
                "sample_id",
                "inferred_case_id",
                "normalized_case_id",
                "signature_score",
                "n_features_used",
                "n_recovered_gene_or_ensembl",
                "recovery_fraction",
            ]
        ].copy()
        layer_small = layer_small.rename(
            columns={
                "dataset_id": "protein_layer_dataset_id",
                "sample_id": "protein_layer_sample_id",
                "signature_score": "protein_layer_signature_score",
                "n_features_used": "protein_layer_n_features_used",
                "n_recovered_gene_or_ensembl": "protein_layer_n_recovered_genes",
                "recovery_fraction": "protein_layer_recovery_fraction",
            }
        )
        rna_small = rna[(rna["cancer"].astype(str) == str(cancer)) & (rna["signature"].astype(str) == str(signature))].copy()
        merged = rna_small.merge(
            layer_small,
            on=["cancer", "signature", "normalized_case_id"],
            how="inner",
            validate="many_to_many",
            suffixes=("_rna", "_protein_layer"),
        )
        if not merged.empty:
            merged["comparison_layer"] = datatype
            rows.append(merged)
        flow_rows.append(
            {
                "cancer": cancer,
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "comparison_layer": datatype,
                "rna_tumor_like_score_rows": int(len(rna_small)),
                "rna_tumor_like_unique_cases": int(rna_small["normalized_case_id"].nunique()),
                "protein_layer_score_rows": int(len(layer_small)),
                "protein_layer_unique_cases": int(layer_small["normalized_case_id"].nunique()),
                "matched_pair_rows": int(len(merged)),
                "matched_unique_cases": int(merged["normalized_case_id"].nunique()) if not merged.empty else 0,
                "rna_duplicate_case_rows": int(rna_small.duplicated(["normalized_case_id"], keep=False).sum()) if not rna_small.empty else 0,
                "protein_layer_duplicate_case_rows": int(layer_small.duplicated(["normalized_case_id"], keep=False).sum()) if not layer_small.empty else 0,
            }
        )
    matched = pd.concat(rows, ignore_index=True, sort=False) if rows else pd.DataFrame()
    if not matched.empty:
        matched = matched.rename(
            columns={
                "dataset_id": "rna_dataset_id",
                "sample_id": "rna_sample_id",
                "signature_score": "rna_signature_score",
                "source_rna": "rna_source",
                "datatype_rna": "rna_datatype",
            }
        )
        matched["signature_label"] = matched["signature"].map(lambda x: SIGNATURE_LABELS.get(x, x))
    return matched, pd.DataFrame(flow_rows)


def spearman_pair(x: pd.Series, y: pd.Series, min_n: int) -> tuple[float, float, int]:
    data = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(data) < min_n or data["x"].nunique() < 3 or data["y"].nunique() < 3:
        return np.nan, np.nan, int(len(data))
    result = spearmanr(data["x"], data["y"])
    return float(result.statistic), float(result.pvalue), int(len(data))


def matched_correlation_summary(matched: pd.DataFrame, min_n: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if matched.empty:
        return pd.DataFrame()
    for keys, group in matched.groupby(["comparison_layer", "cancer", "signature"], dropna=False):
        layer, cancer, signature = keys
        rho, p_value, n = spearman_pair(group["rna_signature_score"], group["protein_layer_signature_score"], min_n)
        rows.append(
            {
                "comparison": f"rna_vs_{layer}",
                "comparison_scope": "within_cancer",
                "comparison_layer": layer,
                "cancer": cancer,
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "n": n,
                "spearman_rho": rho,
                "p_value": p_value,
                "fdr_family": "within_cancer_matched_rna_protein_phospho_correlations_v0_4",
            }
        )
    for keys, group in matched.groupby(["comparison_layer", "signature"], dropna=False):
        layer, signature = keys
        rho, p_value, n = spearman_pair(group["rna_signature_score"], group["protein_layer_signature_score"], min_n)
        rows.append(
            {
                "comparison": f"rna_vs_{layer}",
                "comparison_scope": "pooled_bcm_cancers",
                "comparison_layer": layer,
                "cancer": "pooled_bcm_cancers",
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "n": n,
                "spearman_rho": rho,
                "p_value": p_value,
                "fdr_family": "pooled_matched_rna_protein_phospho_correlations_v0_4",
            }
        )
    corr = pd.DataFrame(rows)
    if not corr.empty:
        corr["fdr_bh"] = np.nan
        for family, idx in corr.groupby("fdr_family").groups.items():
            corr.loc[list(idx), "fdr_bh"] = bh_fdr(pd.to_numeric(corr.loc[list(idx), "p_value"], errors="coerce").to_numpy())
    return corr


def protein_phospho_matched_summary(matched: pd.DataFrame, min_n: int) -> pd.DataFrame:
    if matched.empty:
        return pd.DataFrame()
    prot = matched[matched["comparison_layer"].astype(str) == "proteomics"][
        ["cancer", "signature", "normalized_case_id", "protein_layer_signature_score"]
    ].rename(columns={"protein_layer_signature_score": "protein_score"})
    phos = matched[matched["comparison_layer"].astype(str) == "phosphoproteomics"][
        ["cancer", "signature", "normalized_case_id", "protein_layer_signature_score"]
    ].rename(columns={"protein_layer_signature_score": "phosphoprotein_score"})
    pp = prot.merge(phos, on=["cancer", "signature", "normalized_case_id"], how="inner", validate="many_to_many")
    rows: list[dict[str, Any]] = []
    for keys, group in pp.groupby(["cancer", "signature"], dropna=False):
        cancer, signature = keys
        rho, p_value, n = spearman_pair(group["protein_score"], group["phosphoprotein_score"], min_n)
        rows.append(
            {
                "comparison": "proteomics_vs_phosphoproteomics",
                "comparison_scope": "within_cancer",
                "comparison_layer": "protein_phosphoprotein",
                "cancer": cancer,
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "n": n,
                "spearman_rho": rho,
                "p_value": p_value,
                "fdr_family": "protein_phosphoprotein_matched_case_correlations_v0_4",
            }
        )
    for signature, group in pp.groupby("signature", dropna=False):
        rho, p_value, n = spearman_pair(group["protein_score"], group["phosphoprotein_score"], min_n)
        rows.append(
            {
                "comparison": "proteomics_vs_phosphoproteomics",
                "comparison_scope": "pooled_bcm_cancers",
                "comparison_layer": "protein_phosphoprotein",
                "cancer": "pooled_bcm_cancers",
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "n": n,
                "spearman_rho": rho,
                "p_value": p_value,
                "fdr_family": "protein_phosphoprotein_matched_case_pooled_correlations_v0_4",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh"] = np.nan
        for family, idx in out.groupby("fdr_family").groups.items():
            out.loc[list(idx), "fdr_bh"] = bh_fdr(pd.to_numeric(out.loc[list(idx), "p_value"], errors="coerce").to_numpy())
    return out


def fixed_effect_correlation_meta(corr: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    subset = corr[(corr["comparison_scope"] == "within_cancer") & pd.to_numeric(corr["spearman_rho"], errors="coerce").notna()].copy()
    for keys, group in subset.groupby(["comparison", "comparison_layer", "signature"], dropna=False):
        comparison, layer, signature = keys
        valid = group[pd.to_numeric(group["n"], errors="coerce") > 3].copy()
        if valid.empty:
            continue
        rho = pd.to_numeric(valid["spearman_rho"], errors="coerce").clip(-0.999, 0.999).to_numpy(dtype=float)
        n = pd.to_numeric(valid["n"], errors="coerce").to_numpy(dtype=float)
        z = np.arctanh(rho)
        w = np.maximum(n - 3, 1)
        pooled_z = np.sum(w * z) / np.sum(w)
        pooled_rho = np.tanh(pooled_z)
        q = np.sum(w * (z - pooled_z) ** 2)
        df = max(len(valid) - 1, 1)
        i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
        rows.append(
            {
                "comparison": comparison,
                "comparison_layer": layer,
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "cancer_strata": int(len(valid)),
                "total_matched_n": int(np.nansum(n)),
                "fixed_effect_fisher_z_rho": float(pooled_rho),
                "median_within_cancer_rho": float(np.nanmedian(rho)),
                "positive_cancer_fraction": float(np.mean(rho > 0)),
                "heterogeneity_q_descriptive": float(q),
                "heterogeneity_i2_descriptive": float(i2),
                "interpretation_boundary": "Descriptive fixed-effect summary of within-cancer Spearman correlations; not an external validation or causal meta-analysis.",
            }
        )
    return pd.DataFrame(rows)


def residualize_against_matched_rna(matched: pd.DataFrame, min_n: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    residual_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []
    if matched.empty:
        return pd.DataFrame(), pd.DataFrame()
    for keys, group in matched.groupby(["comparison_layer", "cancer", "signature"], dropna=False):
        layer, cancer, signature = keys
        model_df = group[
            [
                "normalized_case_id",
                "rna_sample_id",
                "protein_layer_sample_id",
                "rna_signature_score",
                "protein_layer_signature_score",
            ]
        ].copy()
        model_df["rna_signature_score"] = pd.to_numeric(model_df["rna_signature_score"], errors="coerce")
        model_df["protein_layer_signature_score"] = pd.to_numeric(model_df["protein_layer_signature_score"], errors="coerce")
        model_df = model_df.dropna()
        if len(model_df) < min_n or model_df["rna_signature_score"].nunique() < 3:
            model_rows.append(
                {
                    "comparison_layer": layer,
                    "cancer": cancer,
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "model_n": int(len(model_df)),
                    "status": "insufficient_matched_samples",
                }
            )
            continue
        y = model_df["protein_layer_signature_score"].to_numpy(dtype=float)
        x = model_df["rna_signature_score"].to_numpy(dtype=float)
        design = np.column_stack([np.ones(len(x)), x])
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        fitted = design @ beta
        residual = y - fitted
        ss_total = np.sum((y - np.mean(y)) ** 2)
        ss_resid = np.sum(residual**2)
        r2 = 1 - ss_resid / ss_total if ss_total > 0 else np.nan
        model_rows.append(
            {
                "comparison_layer": layer,
                "cancer": cancer,
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "model_n": int(len(model_df)),
                "intercept": float(beta[0]),
                "rna_beta": float(beta[1]),
                "rna_adjustment_r2": float(r2),
                "status": "ok",
                "interpretation_boundary": "Linear adjustment of protein/phosphoprotein score for matched RNA score; residuals are descriptive and not mechanistic proof.",
            }
        )
        source_rows = group.loc[model_df.index].copy()
        source_rows["rna_adjusted_residual_score"] = residual
        source_rows["rna_adjustment_intercept"] = float(beta[0])
        source_rows["rna_adjustment_beta"] = float(beta[1])
        source_rows["rna_adjustment_r2"] = float(r2)
        residual_rows.extend(source_rows.to_dict("records"))
    return pd.DataFrame(residual_rows), pd.DataFrame(model_rows)


def residual_context_summary(residuals: pd.DataFrame, context: pd.DataFrame, min_n: int) -> pd.DataFrame:
    if residuals.empty or context.empty:
        return pd.DataFrame()
    ctx = context.copy()
    ctx["normalized_case_id"] = ctx["inferred_case_id"].map(normalize_case_id)
    merged = residuals.merge(ctx, on="normalized_case_id", how="left", validate="many_to_one")
    rows: list[dict[str, Any]] = []
    for keys, group in merged.groupby(["comparison_layer", "cancer", "signature"], dropna=False):
        layer, cancer, signature = keys
        for var in CONTEXT_VARS:
            if var not in group.columns:
                continue
            rho, p_value, n = spearman_pair(group["rna_adjusted_residual_score"], group[var], min_n)
            rows.append(
                {
                    "comparison_layer": layer,
                    "cancer": cancer,
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "context_variable": var,
                    "n": n,
                    "spearman_rho": rho,
                    "p_value": p_value,
                    "fdr_family": "rna_adjusted_residual_context_correlations_v0_4",
                    "interpretation_boundary": "Residual-context association after matched RNA adjustment; not proof of tumor-cell-intrinsic state.",
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh"] = bh_fdr(pd.to_numeric(out["p_value"], errors="coerce").to_numpy())
    return out


def build_sample_inventory(rna_sample_info: pd.DataFrame, cptac_scores: pd.DataFrame, matched_flow: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cancer, group in rna_sample_info.groupby("cancer", dropna=False):
        row = {
            "cancer": cancer,
            "rna_total_samples": int(len(group)),
            "rna_tumor_like_samples": int(group["is_tumor_like_for_matching"].sum()),
            "rna_adjacent_or_normal_samples": int((group["rna_tumor_normal_hint"].astype(str) == "adjacent_or_normal").sum()),
            "rna_unique_cases": int(group["normalized_case_id"].nunique()),
            "rna_duplicate_case_rows": int(group.duplicated(["normalized_case_id"], keep=False).sum()),
        }
        for datatype in ["proteomics", "phosphoproteomics"]:
            layer = cptac_scores[(cptac_scores["cancer"].astype(str) == str(cancer)) & (cptac_scores["datatype"].astype(str) == datatype)]
            flow = matched_flow[(matched_flow["cancer"].astype(str) == str(cancer)) & (matched_flow["comparison_layer"].astype(str) == datatype)]
            row[f"{datatype}_score_rows"] = int(len(layer))
            row[f"{datatype}_unique_cases"] = int(layer["normalized_case_id"].nunique()) if not layer.empty else 0
            row[f"{datatype}_max_matched_unique_cases"] = int(flow["matched_unique_cases"].max()) if not flow.empty else 0
            row[f"{datatype}_median_matched_unique_cases"] = float(flow["matched_unique_cases"].median()) if not flow.empty else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def build_gene_source_proxy(signature_table: pd.DataFrame, depmap_gene: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dep = depmap_gene.copy()
    dep["gene_symbol"] = dep["gene_symbol"].astype(str)
    rows: list[dict[str, Any]] = []
    for _, row in signature_table.iterrows():
        gene = safe_str(row.get("gene_symbol"))
        signature = safe_str(row.get("signature"))
        labels = [label for label, genes in GENE_SOURCE_RULES.items() if gene in genes]
        if not labels:
            labels = ["mixed_or_unassigned_public_marker"]
        dep_row = dep[(dep["signature"].astype(str) == signature) & (dep["gene_symbol"].astype(str) == gene)]
        common_flag = int(dep_row["common_dependency_probability_ge_0_5"].max()) if not dep_row.empty else 0
        global_flag = int(dep_row["global_dependency_top5_percent"].max()) if not dep_row.empty else 0
        rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "gene_symbol": gene,
                "source_proxy_class": ";".join(labels),
                "depmap_common_dependency_probability_ge_0_5": common_flag,
                "depmap_global_dependency_top5_percent": global_flag,
                "evidence_boundary": "Curated marker-class plus DepMap context proxy; not scRNA cell-type deconvolution or tumor-cell-intrinsic validation.",
            }
        )
    gene_source = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for signature, group in gene_source.groupby("signature", dropna=False):
        class_counts = group["source_proxy_class"].str.split(";").explode().value_counts().to_dict()
        summary_rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "signature_gene_count": int(group["gene_symbol"].nunique()),
                "dominant_source_proxy": max(class_counts, key=class_counts.get) if class_counts else "",
                "source_proxy_counts": ";".join(f"{k}:{v}" for k, v in sorted(class_counts.items())),
                "depmap_common_or_global_gene_count": int(
                    (
                        group["depmap_common_dependency_probability_ge_0_5"].astype(int).eq(1)
                        | group["depmap_global_dependency_top5_percent"].astype(int).eq(1)
                    ).sum()
                ),
                "interpretation_boundary": "Source-aware supportive context only; a pan-cancer scRNA expression matrix was not used in this v0.4 scoring layer.",
            }
        )
    return gene_source, pd.DataFrame(summary_rows)


def depmap_lineage_dependency_summary(
    qc_dir: Path,
    harmonization_dir: Path,
    signature_table: pd.DataFrame,
) -> pd.DataFrame:
    inv = pd.read_csv(qc_dir / "depmap_inventory_v0_1.csv")
    model_map = pd.read_csv(harmonization_dir / "depmap_model_map_v0_1.csv")
    effect_path = Path(inv[inv["file_role"].astype(str) == "crispr_gene_effect"]["local_path"].iloc[0])
    dep_path = Path(inv[inv["file_role"].astype(str) == "crispr_dependency"]["local_path"].iloc[0])
    effect_header = pd.read_csv(effect_path, nrows=0).columns.tolist()
    dep_header = pd.read_csv(dep_path, nrows=0).columns.tolist()
    first_col = effect_header[0]
    dep_first_col = dep_header[0]
    genes = set(signature_table["gene_symbol"].astype(str))
    col_to_symbol = {}
    for col in effect_header[1:]:
        match = DEPMAP_GENE_RE.match(col)
        if match and match.group("symbol") in genes:
            col_to_symbol[col] = match.group("symbol")
    selected_cols = [col for col in col_to_symbol if col in effect_header]
    selected_dep_cols = [col for col in selected_cols if col in dep_header]
    effect = pd.read_csv(effect_path, usecols=[first_col] + selected_cols).rename(columns={first_col: "ModelID"})
    dependency = pd.read_csv(dep_path, usecols=[dep_first_col] + selected_dep_cols).rename(columns={dep_first_col: "ModelID"})
    model_cols = ["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"]
    meta = model_map[[c for c in model_cols if c in model_map.columns]].copy()
    effect = effect.merge(meta, on="ModelID", how="left", validate="many_to_one")
    dependency = dependency.merge(meta, on="ModelID", how="left", validate="many_to_one")
    gene_to_signatures = signature_table.groupby("gene_symbol")["signature"].agg(lambda x: sorted(set(map(str, x)))).to_dict()
    rows: list[dict[str, Any]] = []
    for lineage, group in effect.groupby("OncotreeLineage", dropna=False):
        dep_group = dependency.loc[group.index] if len(dependency) == len(effect) else dependency[dependency["ModelID"].isin(group["ModelID"])]
        for signature in SIGNATURES:
            sig_genes = set(signature_table[signature_table["signature"].astype(str) == signature]["gene_symbol"].astype(str))
            cols = [col for col, symbol in col_to_symbol.items() if symbol in sig_genes and col in group.columns]
            if not cols:
                continue
            values = group[cols].apply(pd.to_numeric, errors="coerce")
            dep_values = dep_group[[col for col in cols if col in dep_group.columns]].apply(pd.to_numeric, errors="coerce")
            rows.append(
                {
                    "oncotree_lineage": safe_str(lineage) or "Unknown",
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "model_n": int(len(group)),
                    "mapped_gene_count": int(len(cols)),
                    "median_gene_effect": float(values.stack().median(skipna=True)) if not values.empty else np.nan,
                    "median_dependency_probability": float(dep_values.stack().median(skipna=True)) if not dep_values.empty else np.nan,
                    "genes": ";".join(sorted({col_to_symbol[c] for c in cols})),
                    "interpretation_boundary": "Lineage-stratified DepMap supportive context; not target validation.",
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["signature", "median_gene_effect"], ascending=[True, True]).reset_index(drop=True)


def scrna_support_inventory(qc_dir: Path) -> pd.DataFrame:
    path = qc_dir / "scrna_reference_inventory_v0_1.csv"
    if not path.exists():
        return pd.DataFrame()
    inv = pd.read_csv(path)
    inv["v0_4_use"] = np.where(
        inv["file_kind"].astype(str).eq("h5ad"),
        "inventory_only_not_loaded_for_pan_cancer_source_deconvolution",
        "provenance_or_candidate_index_context",
    )
    inv["interpretation_boundary"] = "scRNA/CELLxGENE reference inventory only in v0.4; no pan-cancer scRNA expression matrix was used for validation."
    return inv


def build_claim_tiers(
    corr: pd.DataFrame,
    meta: pd.DataFrame,
    residual_context: pd.DataFrame,
    depmap_exclusion: pd.DataFrame,
    source_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for signature in SIGNATURES:
        sig_corr = corr[(corr["signature"].astype(str) == signature) & (corr["comparison_scope"] == "within_cancer")]
        sig_meta = meta[meta["signature"].astype(str) == signature] if not meta.empty else pd.DataFrame()
        sig_resid = residual_context[residual_context["signature"].astype(str) == signature] if not residual_context.empty else pd.DataFrame()
        sig_dep = depmap_exclusion[depmap_exclusion["signature"].astype(str) == signature] if not depmap_exclusion.empty else pd.DataFrame()
        sig_source = source_summary[source_summary["signature"].astype(str) == signature] if not source_summary.empty else pd.DataFrame()
        max_n = int(pd.to_numeric(sig_corr["n"], errors="coerce").max()) if not sig_corr.empty else 0
        median_rna_protein = float(
            sig_corr[sig_corr["comparison_layer"].astype(str) == "proteomics"]["spearman_rho"].median()
        ) if not sig_corr.empty else np.nan
        median_rna_phospho = float(
            sig_corr[sig_corr["comparison_layer"].astype(str) == "phosphoproteomics"]["spearman_rho"].median()
        ) if not sig_corr.empty else np.nan
        residual_immune = sig_resid[
            sig_resid["context_variable"].astype(str).isin(["estimate_immune_score", "estimate_tumor_purity"])
        ]
        max_resid_context_abs = float(pd.to_numeric(residual_immune["spearman_rho"], errors="coerce").abs().max()) if not residual_immune.empty else np.nan
        dep_excluded = int(sig_dep["common_or_global_excluded_gene_count"].iloc[0]) if not sig_dep.empty and "common_or_global_excluded_gene_count" in sig_dep.columns else 0
        dep_retained = int(sig_dep["remaining_after_common_or_global_exclusion"].iloc[0]) if not sig_dep.empty and "remaining_after_common_or_global_exclusion" in sig_dep.columns else 0
        dep_total = dep_excluded + dep_retained
        common_fraction = dep_excluded / dep_total if dep_total else 0.0
        if max_n < MIN_CORRELATION_N:
            tier = "insufficient_matched_samples"
            allowed = "Use only as inventory/context; matched RNA-protein evidence is underpowered for manuscript claims."
        elif common_fraction >= 0.5:
            tier = "common_essential_limited_dependency_context"
            allowed = "Matched multi-omic signal can be described, but DepMap dependency interpretation is limited by common/global essential genes."
        elif np.nanmedian([median_rna_protein, median_rna_phospho]) >= 0.35:
            tier = "matched_rna_concordant_context"
            allowed = "Describe as matched RNA-concordant proteogenomic immune-ecology context."
        elif np.isfinite(max_resid_context_abs) and max_resid_context_abs >= 0.25:
            tier = "matched_rna_residual_context"
            allowed = "Describe as cautious RNA-adjusted residual context, with immune/purity overlap explicitly stated."
        else:
            tier = "matched_rna_residual_context"
            allowed = "Describe as modest matched multi-omic residual context; keep in supplement or exploratory results."
        rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "v0_4_claim_tier": tier,
                "max_matched_n": max_n,
                "median_within_cancer_rna_protein_rho": median_rna_protein,
                "median_within_cancer_rna_phosphoprotein_rho": median_rna_phospho,
                "max_abs_rna_adjusted_residual_context_rho": max_resid_context_abs,
                "depmap_common_or_global_excluded_gene_count": dep_excluded,
                "depmap_common_or_global_fraction": common_fraction,
                "dominant_source_proxy": safe_str(sig_source["dominant_source_proxy"].iloc[0]) if not sig_source.empty else "",
                "allowed_claim": allowed,
                "prohibited_claims": "No validation, biomarker, tumor-cell-intrinsic immune-evasion, therapeutic vulnerability, kinase mechanism, or causal claim.",
                "recommended_wording": "matched RNA-concordant context; RNA-adjusted residual context; common-essential dependency context; supportive source-aware triangulation.",
            }
        )
    return pd.DataFrame(rows)


def write_figures(
    fig_dir: Path,
    corr: pd.DataFrame,
    meta: pd.DataFrame,
    residual_context: pd.DataFrame,
    claim: pd.DataFrame,
) -> list[Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    within = corr[(corr["comparison_scope"] == "within_cancer") & corr["comparison_layer"].astype(str).isin(["proteomics", "phosphoproteomics"])]
    for layer, filename, title in [
        ("proteomics", "fig11_bcm_matched_rna_protein_concordance_v0_4.png", "Matched RNA-protein concordance"),
        ("phosphoproteomics", "fig12_bcm_matched_rna_phosphoprotein_concordance_v0_4.png", "Matched RNA-phosphoprotein concordance"),
    ]:
        sub = within[within["comparison_layer"].astype(str) == layer].copy()
        if sub.empty:
            continue
        heat = sub.pivot_table(index="signature", columns="cancer", values="spearman_rho", aggfunc="median")
        heat = heat.reindex([sig for sig in SIGNATURES if sig in heat.index])
        if heat.empty:
            continue
        path = fig_dir / filename
        fig, ax = plt.subplots(figsize=(9.2, 5.2))
        data = heat.to_numpy(dtype=float)
        im = ax.imshow(data, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels([SIGNATURE_LABELS.get(sig, sig) for sig in heat.index], fontsize=8)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if np.isfinite(data[i, j]):
                    ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center", fontsize=7)
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax, shrink=0.82)
        cbar.set_label("Spearman rho")
        fig.tight_layout()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)

    if not claim.empty:
        plot = claim.copy()
        plot["signature_label"] = plot["signature"].map(lambda x: SIGNATURE_LABELS.get(x, x))
        x = pd.to_numeric(plot["median_within_cancer_rna_protein_rho"], errors="coerce")
        y = pd.to_numeric(plot["max_abs_rna_adjusted_residual_context_rho"], errors="coerce")
        size = 60 + 240 * pd.to_numeric(plot["depmap_common_or_global_fraction"], errors="coerce").fillna(0)
        label_offsets = {
            "antigen_presentation_mhc_i": (6, 5),
            "antigen_presentation_mhc_ii": (-38, 9),
            "checkpoint_exhaustion_context": (6, 5),
            "cytolytic_t_cell_context": (8, -10),
            "ifn_gamma_response": (8, -2),
            "myeloid_inflammatory_context": (6, 5),
            "proteasome_antigen_processing": (7, -4),
            "tgfb_emt_exclusion_context": (6, 6),
        }
        path = fig_dir / "fig13_claim_map_rna_concordance_residual_context_v0_4.png"
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        ax.scatter(x, y, s=size, alpha=0.72, color="#4C78A8", edgecolor="black", linewidth=0.6)
        for _, row in plot.iterrows():
            if np.isfinite(row["median_within_cancer_rna_protein_rho"]) and np.isfinite(row["max_abs_rna_adjusted_residual_context_rho"]):
                dx, dy = label_offsets.get(safe_str(row["signature"]), (5, 4))
                ax.annotate(
                    SIGNATURE_LABELS.get(row["signature"], row["signature"]),
                    (row["median_within_cancer_rna_protein_rho"], row["max_abs_rna_adjusted_residual_context_rho"]),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    fontsize=8,
                )
        ax.axvline(0.35, color="grey", linestyle="--", linewidth=1)
        ax.axhline(0.25, color="grey", linestyle="--", linewidth=1)
        finite_x = x[np.isfinite(x)]
        finite_y = y[np.isfinite(y)]
        if len(finite_x):
            ax.set_xlim(max(0.25, float(finite_x.min()) - 0.04), min(0.95, float(finite_x.max()) + 0.07))
        if len(finite_y):
            ax.set_ylim(max(0.22, float(finite_y.min()) - 0.03), min(0.58, float(finite_y.max()) + 0.04))
        ax.set_xlabel("Median within-cancer RNA-protein rho")
        ax.set_ylabel("Max |RNA-adjusted residual context rho|")
        ax.set_title("Claim map for matched multi-omic context")
        fig.tight_layout()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)

    if not meta.empty:
        sub = meta[meta["comparison_layer"].astype(str).isin(["proteomics", "phosphoproteomics"])].copy()
        if not sub.empty:
            sub["signature_label"] = sub["signature"].map(lambda x: SIGNATURE_LABELS.get(x, x))
            pivot = sub.pivot_table(index="signature_label", columns="comparison_layer", values="fixed_effect_fisher_z_rho", aggfunc="median")
            path = fig_dir / "fig14_fixed_effect_matched_correlation_summary_v0_4.png"
            fig, ax = plt.subplots(figsize=(7.4, 5.2))
            pivot = pivot.reindex([SIGNATURE_LABELS.get(sig, sig) for sig in SIGNATURES if SIGNATURE_LABELS.get(sig, sig) in pivot.index])
            pivot.plot(kind="barh", ax=ax, color=["#72B7B2", "#F58518"], width=0.74)
            ax.axvline(0, color="black", linewidth=0.8)
            ax.set_xlabel("Fixed-effect descriptive rho")
            ax.set_ylabel("")
            ax.set_title("Cancer-level matched correlation summary")
            ax.legend(title="")
            fig.tight_layout()
            fig.savefig(path, dpi=300)
            plt.close(fig)
            figures.append(path)
    return figures


def write_report(
    report_path: Path,
    summary: dict[str, Any],
    inventory: pd.DataFrame,
    corr: pd.DataFrame,
    meta: pd.DataFrame,
    residual_models: pd.DataFrame,
    residual_context: pd.DataFrame,
    claim: pd.DataFrame,
) -> None:
    lines = [
        "# Discovery v0.4 matched RNA residual and source-aware context",
        "",
        "This run integrates downloaded CPTAC BCM transcriptomics with existing BCM protein and phosphoprotein signature scores. It uses explicit case-ID normalization and match-flow tables before any correlation or residual analysis. Results are supportive context only.",
        "",
        "## Run Summary",
    ]
    for key in [
        "rna_manifest_rows",
        "rna_score_rows",
        "matched_pair_rows",
        "matched_correlation_rows",
        "protein_phospho_correlation_rows",
        "residual_model_rows",
        "residual_context_rows",
        "claim_tier_rows",
    ]:
        lines.append(f"- {key}: {summary.get(key)}")
    lines.append("")
    lines.append("## Sample Inventory")
    lines.append(inventory.to_markdown(index=False) if not inventory.empty else "No inventory rows.")
    lines.append("")
    lines.append("## Matched Correlation Summary")
    if not corr.empty:
        show = corr[(corr["comparison_scope"] == "pooled_bcm_cancers")][
            ["comparison", "signature_label", "n", "spearman_rho", "fdr_bh"]
        ]
        lines.append(show.to_markdown(index=False))
    else:
        lines.append("No matched correlations.")
    lines.append("")
    lines.append("## Descriptive Cancer-Level Meta Summary")
    if not meta.empty:
        lines.append(
            meta[
                [
                    "comparison",
                    "signature_label",
                    "cancer_strata",
                    "total_matched_n",
                    "fixed_effect_fisher_z_rho",
                    "positive_cancer_fraction",
                    "heterogeneity_i2_descriptive",
                ]
            ].to_markdown(index=False)
        )
    else:
        lines.append("No meta summary rows.")
    lines.append("")
    lines.append("## RNA Adjustment")
    ok_models = residual_models[residual_models["status"].astype(str) == "ok"] if not residual_models.empty else pd.DataFrame()
    lines.append(f"- RNA-adjusted residual models with status ok: {len(ok_models)}")
    if not residual_context.empty:
        top = residual_context.sort_values("spearman_rho", key=lambda s: pd.to_numeric(s, errors="coerce").abs(), ascending=False).head(20)
        lines.append(top[["comparison_layer", "cancer", "signature_label", "context_variable", "n", "spearman_rho", "fdr_bh"]].to_markdown(index=False))
    lines.append("")
    lines.append("## Claim Tiers")
    if not claim.empty:
        lines.append(claim[["signature_label", "v0_4_claim_tier", "allowed_claim"]].to_markdown(index=False))
    lines.append("")
    lines.append("## Boundary")
    lines.append("Do not write this as validation, biomarker discovery, tumor-cell-intrinsic immune evasion, therapeutic vulnerability, kinase mechanism, or causal evidence. The safest wording is matched RNA-concordant proteogenomic immune ecology with RNA-adjusted residual and common-essential dependency context.")
    lines.append("")
    lines.append("## Outputs")
    for name, value in summary.get("outputs", {}).items():
        lines.append(f"- {name}: `{value}`")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    start = time.time()
    args = parse_args()
    np.random.seed(args.random_seed)
    root = Path(args.project_root)
    qc_dir = resolve_path(args.qc_dir, root)
    harmonization_dir = resolve_path(args.harmonization_dir, root)
    discovery_v0_1_dir = resolve_path(args.discovery_v0_1_dir, root)
    discovery_v0_2_dir = resolve_path(args.discovery_v0_2_dir, root)
    discovery_v0_3_dir = resolve_path(args.discovery_v0_3_dir, root)
    rna_manifest_path = resolve_path(args.rna_manifest, root)
    out_dir = resolve_path(args.out_dir, root)
    fig_dir = resolve_path(args.figure_dir, root)
    report_path = resolve_path(args.report, root)
    summary_path = resolve_path(args.summary, root)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    signature_table = pd.read_csv(discovery_v0_1_dir / "signature_gene_sets_v0_1.csv")
    signature_map = build_signature_ensembl_map(signature_table)
    rna_manifest = pd.read_csv(rna_manifest_path)
    cptac_scores = prepare_existing_cptac_scores(pd.read_csv(discovery_v0_1_dir / "cptac_signature_scores_long_v0_1.csv"))
    immune_context = pd.read_csv(discovery_v0_1_dir / "cptac_immune_context_scores_v0_1.csv")
    depmap_gene = pd.read_csv(discovery_v0_3_dir / "depmap_signature_gene_level_common_essential_audit_v0_3.csv")
    depmap_exclusion = pd.read_csv(discovery_v0_3_dir / "depmap_common_essential_exclusion_summary_v0_3.csv")

    rna_score_parts = []
    rna_recovery_parts = []
    rna_sample_parts = []
    for _, row in rna_manifest.iterrows():
        print(f"Scoring RNA signatures for {row['cancer']}: {row['data_file']}")
        scores, recovery, samples = score_rna_file(row, root, signature_table, signature_map, args.min_rna_mapped_genes)
        rna_score_parts.append(scores)
        rna_recovery_parts.append(recovery)
        rna_sample_parts.append(samples)
    rna_scores = pd.concat(rna_score_parts, ignore_index=True, sort=False) if rna_score_parts else pd.DataFrame()
    rna_recovery = pd.concat(rna_recovery_parts, ignore_index=True, sort=False) if rna_recovery_parts else pd.DataFrame()
    rna_samples = pd.concat(rna_sample_parts, ignore_index=True, sort=False) if rna_sample_parts else pd.DataFrame()

    matched, match_flow = build_matched_table(rna_scores, cptac_scores)
    sample_inventory = build_sample_inventory(rna_samples, cptac_scores, match_flow)
    corr = matched_correlation_summary(matched, args.min_correlation_n)
    pp_corr = protein_phospho_matched_summary(matched, args.min_correlation_n)
    all_corr = pd.concat([corr, pp_corr], ignore_index=True, sort=False) if not pp_corr.empty else corr
    meta = fixed_effect_correlation_meta(all_corr)
    residuals, residual_models = residualize_against_matched_rna(matched, args.min_residual_n)
    residual_ctx = residual_context_summary(residuals, immune_context, args.min_correlation_n)
    gene_source, source_summary = build_gene_source_proxy(signature_table, depmap_gene)
    lineage_depmap = depmap_lineage_dependency_summary(qc_dir, harmonization_dir, signature_table)
    scrna_inventory = scrna_support_inventory(qc_dir)
    claim = build_claim_tiers(all_corr, meta, residual_ctx, depmap_exclusion, source_summary)
    figures = write_figures(fig_dir, all_corr, meta, residual_ctx, claim)

    outputs = {
        "rna_signature_scores": out_dir / "bcm_rna_signature_scores_long_v0_4.csv",
        "rna_signature_recovery": out_dir / "bcm_rna_signature_recovery_v0_4.csv",
        "rna_sample_inventory": out_dir / "bcm_rna_sample_inventory_v0_4.csv",
        "sample_flow": out_dir / "bcm_matched_sample_flow_v0_4.csv",
        "matched_scores": out_dir / "bcm_matched_rna_protein_phosphoprotein_scores_v0_4.csv",
        "matched_correlations": out_dir / "bcm_matched_correlation_summary_v0_4.csv",
        "matched_correlation_meta": out_dir / "bcm_matched_correlation_meta_summary_v0_4.csv",
        "rna_adjusted_residuals": out_dir / "bcm_rna_adjusted_residual_scores_v0_4.csv",
        "rna_adjusted_residual_models": out_dir / "bcm_rna_adjusted_residual_models_v0_4.csv",
        "rna_adjusted_residual_context": out_dir / "bcm_rna_adjusted_residual_context_summary_v0_4.csv",
        "gene_source_proxy": out_dir / "signature_gene_source_proxy_v0_4.csv",
        "source_proxy_summary": out_dir / "signature_source_proxy_summary_v0_4.csv",
        "depmap_lineage_dependency": out_dir / "depmap_lineage_dependency_summary_v0_4.csv",
        "scrna_support_inventory": out_dir / "scrna_support_inventory_v0_4.csv",
        "claim_tiers": out_dir / "discovery_claim_tier_table_v0_4.csv",
    }
    rna_scores.to_csv(outputs["rna_signature_scores"], index=False)
    rna_recovery.to_csv(outputs["rna_signature_recovery"], index=False)
    sample_inventory.to_csv(outputs["rna_sample_inventory"], index=False)
    match_flow.to_csv(outputs["sample_flow"], index=False)
    matched.to_csv(outputs["matched_scores"], index=False)
    all_corr.to_csv(outputs["matched_correlations"], index=False)
    meta.to_csv(outputs["matched_correlation_meta"], index=False)
    residuals.to_csv(outputs["rna_adjusted_residuals"], index=False)
    residual_models.to_csv(outputs["rna_adjusted_residual_models"], index=False)
    residual_ctx.to_csv(outputs["rna_adjusted_residual_context"], index=False)
    gene_source.to_csv(outputs["gene_source_proxy"], index=False)
    source_summary.to_csv(outputs["source_proxy_summary"], index=False)
    lineage_depmap.to_csv(outputs["depmap_lineage_dependency"], index=False)
    scrna_inventory.to_csv(outputs["scrna_support_inventory"], index=False)
    claim.to_csv(outputs["claim_tiers"], index=False)

    png_qa = []
    for fig in figures:
        png_qa.append(
            {
                "figure": rel(fig, root),
                "exists": fig.exists(),
                "size_bytes": fig.stat().st_size if fig.exists() else 0,
                "dpi_target": 300,
            }
        )

    summary = {
        "version": args.version,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "elapsed_seconds": round(time.time() - start, 2),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "random_seed": args.random_seed,
        "min_correlation_n": args.min_correlation_n,
        "min_residual_n": args.min_residual_n,
        "rna_manifest_rows": int(len(rna_manifest)),
        "rna_score_rows": int(len(rna_scores)),
        "rna_recovery_rows": int(len(rna_recovery)),
        "matched_pair_rows": int(len(matched)),
        "matched_unique_cases": int(matched["normalized_case_id"].nunique()) if not matched.empty else 0,
        "matched_correlation_rows": int(len(all_corr)),
        "protein_phospho_correlation_rows": int(len(pp_corr)),
        "meta_summary_rows": int(len(meta)),
        "residual_score_rows": int(len(residuals)),
        "residual_model_rows": int(len(residual_models)),
        "residual_context_rows": int(len(residual_ctx)),
        "gene_source_proxy_rows": int(len(gene_source)),
        "depmap_lineage_rows": int(len(lineage_depmap)),
        "claim_tier_rows": int(len(claim)),
        "claim_tier_counts": claim["v0_4_claim_tier"].value_counts().to_dict() if not claim.empty else {},
        "figures": [rel(fig, root) for fig in figures],
        "figure_qa": png_qa,
        "outputs": {key: rel(value, root) for key, value in outputs.items()},
        "boundary": "Supportive matched multi-omic context only; no validation, biomarker, therapeutic vulnerability, tumor-cell-intrinsic immune-evasion, kinase mechanism, or causal claim.",
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    write_report(report_path, summary, sample_inventory, all_corr, meta, residual_models, residual_ctx, claim)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
