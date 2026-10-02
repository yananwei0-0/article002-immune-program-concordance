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
import math
import platform
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, spearmanr

from run_pan_cancer_proteomic_immune_evasion_discovery_v0_1 import (
    CONTEXT_VARS,
    SIGNATURES,
    build_immune_context,
    ensg_stable,
    infer_cptac_case_id,
    parse_feature_gene,
    read_header,
    rel,
    resolve_path,
    safe_str,
    zscore_rows,
)


PROJECT_ROOT = Path(
    _workspace_os.environ.get("ARTICLE002_DATA_ROOT", str(Path.cwd()))
).expanduser().resolve()
VERSION = "v0_2"
RANDOM_CONTEXT_VARS = ["estimate_immune_score", "estimate_tumor_purity"]
TCGA_ANCHOR_VARS = [
    "IFNG_score_21050467",
    "IFN_21978456",
    "Module3_IFN_score",
    "MHC1_21978456",
    "MHC2_21978456",
    "TcClassII_score",
    "Tcell_21978456",
    "CD8A",
    "CD68",
    "CSF1_response",
    "TGFB_score_21050467",
    "TGFB_PCA_17349583",
    "PD1_PDL1_score",
]
SIGNATURE_ANCHORS = {
    "antigen_presentation_mhc_i": ["MHC1_21978456", "IFNG_score_21050467"],
    "antigen_presentation_mhc_ii": ["MHC2_21978456", "TcClassII_score"],
    "ifn_gamma_response": ["IFNG_score_21050467", "IFN_21978456", "Module3_IFN_score"],
    "cytolytic_t_cell_context": ["CD8A", "Tcell_21978456"],
    "checkpoint_exhaustion_context": ["PD1_PDL1_score", "Tcell_21978456"],
    "myeloid_inflammatory_context": ["CD68", "CSF1_response"],
    "tgfb_emt_exclusion_context": ["TGFB_score_21050467", "TGFB_PCA_17349583"],
    "proteasome_antigen_processing": ["MHC1_21978456", "IFN_21978456"],
}
SIGNATURE_LABELS = {
    "antigen_presentation_mhc_i": "MHC-I",
    "antigen_presentation_mhc_ii": "MHC-II",
    "ifn_gamma_response": "IFN-gamma",
    "cytolytic_t_cell_context": "Cytolytic",
    "checkpoint_exhaustion_context": "Checkpoint",
    "myeloid_inflammatory_context": "Myeloid",
    "tgfb_emt_exclusion_context": "TGF/EMT",
    "proteasome_antigen_processing": "Proteasome",
}
ENTREZ_RE = re.compile(r"^(\d+)$")
DEPMAP_GENE_RE = re.compile(r"^(?P<symbol>.+?) \((?P<entrez>\d+)\)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="v0.2 controls, TCGA RNA comparator, and DepMap supportive context."
    )
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--version", default=VERSION)
    parser.add_argument("--harmonization-dir", default="04_processed/harmonization_v0_1")
    parser.add_argument("--qc-dir", default="04_processed/first_pass_qc")
    parser.add_argument("--discovery-v0-1-dir", default="04_processed/discovery_v0_1")
    parser.add_argument("--out-dir", default="04_processed/discovery_v0_2")
    parser.add_argument("--figure-dir", default="06_figures/discovery_v0_2")
    parser.add_argument("--report", default="05_reports/discovery_v0_2_report.md")
    parser.add_argument("--summary", default="05_reports/discovery_v0_2_summary.json")
    parser.add_argument("--random-seed", type=int, default=20260624)
    parser.add_argument("--n-random-controls", type=int, default=500)
    parser.add_argument("--min-correlation-n", type=int, default=25)
    parser.add_argument("--min-valid-feature-fraction", type=float, default=0.5)
    return parser.parse_args()


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def read_tabular_header(path: Path) -> list[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        return next(reader)


def signature_maps(signature_table: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    symbol_map: dict[str, set[str]] = defaultdict(set)
    ensembl_map: dict[str, set[str]] = defaultdict(set)
    for _, row in signature_table.iterrows():
        signature = safe_str(row["signature"])
        gene = safe_str(row["gene_symbol"])
        if gene:
            symbol_map[gene].add(signature)
        for stable in safe_str(row.get("ensembl_stable_ids_from_cptac", "")).split(";"):
            stable = ensg_stable(stable)
            if stable:
                ensembl_map[stable].add(signature)
    return symbol_map, ensembl_map


def signature_to_ensembl(signature_table: pd.DataFrame) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for _, row in signature_table.iterrows():
        signature = safe_str(row["signature"])
        for stable in safe_str(row.get("ensembl_stable_ids_from_cptac", "")).split(";"):
            stable = ensg_stable(stable)
            if stable:
                out[signature].add(stable)
    return out


def infer_ann_count(item: pd.Series) -> int:
    source = safe_str(item.get("source"))
    ann_count = int(item.get("annotation_column_count", 0) or 0)
    if source == "bcm" and ann_count < 1:
        ann_count = 1
    if source == "umich" and ann_count < 5:
        ann_count = 5
    return ann_count


def load_cptac_z_matrix(
    root: Path,
    item: pd.Series,
    min_valid_feature_fraction: float,
) -> tuple[list[str], np.ndarray, list[set[str]], np.ndarray]:
    path = resolve_path(safe_str(item["local_path"]), root)
    source = safe_str(item["source"])
    ann_count = infer_ann_count(item)
    header = read_header(path, "\t")
    if not header:
        raise ValueError(f"empty header: {path}")
    sample_ids = [safe_str(x) for x in header[ann_count:]]
    names = list(range(len(header)))
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
    if matrix.empty or not sample_ids:
        raise ValueError(f"empty matrix: {path}")
    raw = matrix.iloc[:, 0].astype(str).tolist()
    gene_keys: list[set[str]] = []
    for idx, feature_id in enumerate(raw):
        row_values = matrix.iloc[idx, :ann_count].tolist()
        symbol, stable, _phosphosite = parse_feature_gene(feature_id, row_values, header, source)
        keys = set()
        if symbol:
            keys.add(symbol)
        if stable:
            keys.add(stable)
        gene_keys.append(keys)
    numeric = matrix.iloc[:, ann_count:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    valid_fraction = np.isfinite(numeric).mean(axis=1)
    z = zscore_rows(numeric)
    valid_rows = (
        (valid_fraction >= min_valid_feature_fraction)
        & np.isfinite(z).any(axis=1)
        & (np.nanstd(z, axis=1) > 0)
    )
    return sample_ids, z, gene_keys, valid_rows


def signature_masks(
    gene_keys: list[set[str]],
    signature_symbol_map: dict[str, set[str]],
    signature_ensembl_map: dict[str, set[str]],
) -> dict[str, np.ndarray]:
    masks = {signature: np.zeros(len(gene_keys), dtype=bool) for signature in SIGNATURES}
    for idx, keys in enumerate(gene_keys):
        hit_signatures = set()
        for key in keys:
            hit_signatures.update(signature_symbol_map.get(key, set()))
            hit_signatures.update(signature_ensembl_map.get(key, set()))
        for signature in hit_signatures:
            if signature in masks:
                masks[signature][idx] = True
    return masks


def vectorized_spearman(score_matrix: np.ndarray, context: np.ndarray) -> np.ndarray:
    finite_context = np.isfinite(context)
    if finite_context.sum() < 3:
        return np.full(score_matrix.shape[0], np.nan)
    scores = score_matrix[:, finite_context]
    ctx = context[finite_context]
    out = np.full(scores.shape[0], np.nan)
    ctx_rank = rankdata(ctx)
    ctx_centered = ctx_rank - ctx_rank.mean()
    ctx_ss = np.sqrt(np.sum(ctx_centered**2))
    if ctx_ss == 0:
        return out
    for i in range(scores.shape[0]):
        row = scores[i, :]
        finite = np.isfinite(row)
        if finite.sum() < 3:
            continue
        if finite.all():
            row_rank = rankdata(row)
            row_centered = row_rank - row_rank.mean()
            row_ss = np.sqrt(np.sum(row_centered**2))
            if row_ss > 0:
                out[i] = float(np.sum(row_centered * ctx_centered) / (row_ss * ctx_ss))
        else:
            local_ctx = ctx[finite]
            local_row = row[finite]
            if len(np.unique(local_ctx)) < 3 or len(np.unique(local_row)) < 3:
                continue
            out[i] = float(spearmanr(local_row, local_ctx).correlation)
    return out


def empirical_p_and_percentile(observed: float, controls: np.ndarray) -> tuple[float, float, float, float]:
    controls = np.asarray(controls, dtype=float)
    controls = controls[np.isfinite(controls)]
    if not np.isfinite(observed) or len(controls) == 0:
        return np.nan, np.nan, np.nan, np.nan
    abs_obs = abs(observed)
    abs_controls = np.abs(controls)
    p = (1 + np.sum(abs_controls >= abs_obs)) / (len(abs_controls) + 1)
    percentile = 100.0 * np.mean(abs_controls <= abs_obs)
    z = (abs_obs - np.nanmean(abs_controls)) / np.nanstd(abs_controls, ddof=1) if len(abs_controls) > 1 else np.nan
    delta = abs_obs - np.nanmedian(abs_controls)
    return float(p), float(percentile), float(z), float(delta)


def run_random_control_calibration(
    root: Path,
    qc_dir: Path,
    discovery_dir: Path,
    harmonization_dir: Path,
    out_dir: Path,
    n_controls: int,
    random_seed: int,
    min_n: int,
    min_valid_feature_fraction: float,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    cptac_inv = pd.read_csv(qc_dir / "cptac_processed_inventory_v0_1.csv")
    signature_table = pd.read_csv(discovery_dir / "signature_gene_sets_v0_1.csv")
    signature_symbol_map, stable_to_signature = signature_maps(signature_table)
    signature_stable_map = signature_to_ensembl(signature_table)
    recovery = pd.read_csv(discovery_dir / "cptac_signature_recovery_v0_1.csv")
    observed_corr = pd.read_csv(discovery_dir / "cptac_signature_context_correlations_v0_1.csv")
    context = pd.read_csv(discovery_dir / "cptac_immune_context_scores_v0_1.csv")
    context = context.drop_duplicates("inferred_case_id")

    discovery_mats = cptac_inv[
        cptac_inv["source"].isin(["bcm", "umich"])
        & cptac_inv["datatype"].isin(["proteomics", "phosphoproteomics"])
        & cptac_inv["tumor_normal_hint"].isin(["tumor", "unknown"])
    ].copy()
    rows: list[dict[str, Any]] = []
    for _, item in discovery_mats.iterrows():
        dataset_id = safe_str(item["dataset_id"])
        print(f"random-control matrix {dataset_id}", flush=True)
        try:
            sample_ids, z, gene_keys, valid_rows = load_cptac_z_matrix(root, item, min_valid_feature_fraction)
        except Exception as exc:
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "source": safe_str(item["source"]),
                    "cancer": safe_str(item["cancer"]),
                    "datatype": safe_str(item["datatype"]),
                    "tumor_normal_hint": safe_str(item["tumor_normal_hint"]),
                    "signature": "",
                    "context_variable": "",
                    "status": f"matrix_skipped:{type(exc).__name__}",
                    "n_controls_requested": n_controls,
                    "n_controls_evaluated": 0,
                    "random_seed": random_seed,
                }
            )
            continue
        masks = signature_masks(gene_keys, signature_symbol_map, stable_to_signature)
        sample_context = pd.DataFrame(
            {
                "sample_id": sample_ids,
                "inferred_case_id": [infer_cptac_case_id(x) for x in sample_ids],
            }
        ).merge(context, on="inferred_case_id", how="left", validate="many_to_one")
        valid_pool = np.where(valid_rows)[0]
        for signature in SIGNATURES:
            true_mask = masks[signature] & valid_rows
            k = int(true_mask.sum())
            recovery_row = recovery[
                (recovery["dataset_id"].astype(str) == dataset_id)
                & (recovery["signature"].astype(str) == signature)
            ]
            requested_gene_count = (
                int(recovery_row["requested_gene_count"].iloc[0]) if len(recovery_row) else len(SIGNATURES[signature]["genes"])
            )
            recovered_gene_count = (
                int(recovery_row["recovered_gene_or_ensembl_count"].iloc[0]) if len(recovery_row) else np.nan
            )
            if k < 3:
                continue
            pool = np.setdiff1d(valid_pool, np.where(true_mask)[0], assume_unique=False)
            if len(pool) < max(k * 2, 50):
                pool = valid_pool
            sampled_scores = []
            for _i in range(n_controls):
                chosen = rng.choice(pool, size=k, replace=False if len(pool) >= k else True)
                sampled_scores.append(np.nanmean(z[chosen, :], axis=0))
            control_matrix = np.vstack(sampled_scores)
            obs_sub = observed_corr[
                (observed_corr["source"].astype(str) == safe_str(item["source"]))
                & (observed_corr["cancer"].astype(str) == safe_str(item["cancer"]))
                & (observed_corr["datatype"].astype(str) == safe_str(item["datatype"]))
                & (observed_corr["tumor_normal_hint"].astype(str) == safe_str(item["tumor_normal_hint"]))
                & (observed_corr["signature"].astype(str) == signature)
            ]
            for var in RANDOM_CONTEXT_VARS:
                if var not in sample_context.columns:
                    continue
                context_values = pd.to_numeric(sample_context[var], errors="coerce").to_numpy(dtype=float)
                valid_context = np.isfinite(context_values)
                if valid_context.sum() < min_n or len(np.unique(context_values[valid_context])) < 3:
                    status = "insufficient_context"
                    true_rho = np.nan
                    controls = np.array([], dtype=float)
                else:
                    controls = vectorized_spearman(control_matrix, context_values)
                    true_row = obs_sub[obs_sub["context_variable"].astype(str) == var]
                    true_rho = float(true_row["spearman_rho"].iloc[0]) if len(true_row) else np.nan
                    status = "ok" if np.isfinite(true_rho) and np.isfinite(controls).sum() else "missing_observed_or_controls"
                p_emp, pct, z_emp, delta = empirical_p_and_percentile(true_rho, controls)
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "source": safe_str(item["source"]),
                        "cancer": safe_str(item["cancer"]),
                        "datatype": safe_str(item["datatype"]),
                        "tumor_normal_hint": safe_str(item["tumor_normal_hint"]),
                        "signature": signature,
                        "context_variable": var,
                        "status": status,
                        "requested_gene_count": requested_gene_count,
                        "recovered_gene_or_ensembl_count": recovered_gene_count,
                        "observed_feature_count": k,
                        "valid_feature_pool_size": int(len(valid_pool)),
                        "control_pool_size": int(len(pool)),
                        "observed_spearman_rho": true_rho,
                        "control_median_abs_rho": float(np.nanmedian(np.abs(controls))) if len(controls) else np.nan,
                        "control_mean_abs_rho": float(np.nanmean(np.abs(controls))) if len(controls) else np.nan,
                        "control_sd_abs_rho": float(np.nanstd(np.abs(controls), ddof=1)) if len(controls) > 1 else np.nan,
                        "observed_minus_control_median_abs_rho": delta,
                        "empirical_two_sided_p_abs_rho": p_emp,
                        "observed_abs_rho_percentile_vs_controls": pct,
                        "observed_abs_rho_z_vs_controls": z_emp,
                        "n_context_samples": int(valid_context.sum()) if var in sample_context.columns else 0,
                        "n_controls_requested": n_controls,
                        "n_controls_evaluated": int(np.isfinite(controls).sum()) if len(controls) else 0,
                        "random_seed": random_seed,
                        "control_definition": "feature_count_matched_within_same_cptac_matrix_excluding_curated_features_when_possible",
                    }
                )
    result = pd.DataFrame(rows)
    if not result.empty and "empirical_two_sided_p_abs_rho" in result.columns:
        ok = result["status"].eq("ok")
        result["empirical_fdr_bh_within_random_calibration"] = np.nan
        result.loc[ok, "empirical_fdr_bh_within_random_calibration"] = bh_fdr(
            result.loc[ok, "empirical_two_sided_p_abs_rho"].to_numpy()
        )
    result.to_csv(out_dir / "cptac_random_control_calibration_v0_2.csv", index=False)
    return result


def depmap_symbol_entrez_map(depmap_feature_map: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    symbol_to_entrez: dict[str, set[str]] = defaultdict(set)
    entrez_to_symbol: dict[str, set[str]] = defaultdict(set)
    for _, row in depmap_feature_map.iterrows():
        symbol = safe_str(row.get("gene_symbol"))
        entrez = safe_str(row.get("entrez_id"))
        if "." in entrez:
            try:
                entrez = str(int(float(entrez)))
            except Exception:
                pass
        if symbol and entrez and entrez.lower() != "nan":
            symbol_to_entrez[symbol].add(entrez)
            entrez_to_symbol[entrez].add(symbol)
    return symbol_to_entrez, entrez_to_symbol


def read_selected_tcga_rows(path: Path, wanted_row_ids: set[str], sample_ids: list[str]) -> tuple[pd.DataFrame, list[str]]:
    header = read_tabular_header(path)
    sample_to_index = {sample: idx for idx, sample in enumerate(header[1:])}
    indices = [sample_to_index[sample] for sample in sample_ids if sample in sample_to_index]
    kept_samples = [sample for sample in sample_ids if sample in sample_to_index]
    rows = []
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if not row:
                continue
            row_id = safe_str(row[0])
            if row_id not in wanted_row_ids:
                continue
            values = []
            for idx in indices:
                try:
                    values.append(float(row[idx + 1]))
                except Exception:
                    values.append(np.nan)
            rows.append((row_id, values))
    data = pd.DataFrame([r[1] for r in rows], index=[r[0] for r in rows], columns=kept_samples)
    return data, kept_samples


def read_tcga_immune_anchors(path: Path, sample_ids: list[str]) -> pd.DataFrame:
    data, kept_samples = read_selected_tcga_rows(path, set(TCGA_ANCHOR_VARS), sample_ids)
    if data.empty:
        return pd.DataFrame({"sample_id": kept_samples})
    out = data.T.reset_index().rename(columns={"index": "sample_id"})
    out.columns = [safe_str(c) for c in out.columns]
    return out


def run_tcga_rna_comparator(
    root: Path,
    harmonization_dir: Path,
    out_dir: Path,
    min_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tcga_map = pd.read_csv(harmonization_dir / "tcga_sample_map_v0_1.csv")
    depmap_feature_map = pd.read_csv(harmonization_dir / "depmap_gene_feature_map_v0_1.csv")
    symbol_to_entrez, entrez_to_symbol = depmap_symbol_entrez_map(depmap_feature_map)
    rna_map = tcga_map[
        (tcga_map["layer"].astype(str) == "rna_expression")
        & (tcga_map["sample_type_code"].astype(int) == 1)
    ].copy()
    rna_map = rna_map.sort_values(["patient_id", "sample_id"]).drop_duplicates("patient_id", keep="first")
    sample_ids = rna_map["sample_id"].astype(str).tolist()
    sample_meta = rna_map[
        ["sample_id", "patient_id", "cancer_hint", "primary_disease_hint", "sample_type_code", "sample_type_label"]
    ].copy()
    target_row_ids: set[str] = set()
    signature_row_ids: dict[str, set[str]] = defaultdict(set)
    row_id_to_symbol: dict[str, set[str]] = defaultdict(set)
    for signature, spec in SIGNATURES.items():
        for gene in spec["genes"]:
            target_row_ids.add(gene)
            signature_row_ids[signature].add(gene)
            row_id_to_symbol[gene].add(gene)
            for entrez in symbol_to_entrez.get(gene, set()):
                target_row_ids.add(entrez)
                signature_row_ids[signature].add(entrez)
                row_id_to_symbol[entrez].add(gene)
    rna_path = root / "03_downloads/tcga_xena_pancanatlas/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    expr, kept_samples = read_selected_tcga_rows(rna_path, target_row_ids, sample_ids)
    sample_meta = sample_meta[sample_meta["sample_id"].isin(kept_samples)].set_index("sample_id").loc[kept_samples].reset_index()
    expr = expr.loc[:, kept_samples]
    expr_z = pd.DataFrame(zscore_rows(expr.to_numpy(dtype=float)), index=expr.index, columns=expr.columns)
    score_rows: list[dict[str, Any]] = []
    recovery_rows: list[dict[str, Any]] = []
    for signature, spec in SIGNATURES.items():
        row_ids = sorted(signature_row_ids.get(signature, set()).intersection(set(expr_z.index.astype(str))))
        mapped_symbols = sorted({s for row_id in row_ids for s in row_id_to_symbol.get(row_id, set()) if s in spec["genes"]})
        recovery_rows.append(
            {
                "signature": signature,
                "requested_gene_count": len(spec["genes"]),
                "mapped_expression_row_count": len(row_ids),
                "mapped_gene_symbol_count": len(mapped_symbols),
                "mapped_gene_symbols": ";".join(mapped_symbols),
                "mapping_source": "TCGA RNA row IDs matched by gene symbol first and DepMap-derived Entrez IDs as fallback",
            }
        )
        if not row_ids:
            continue
        scores = expr_z.loc[row_ids].mean(axis=0, skipna=True)
        for sample_id, value in scores.items():
            score_rows.append({"sample_id": sample_id, "signature": signature, "rna_signature_score": value})
    scores_long = pd.DataFrame(score_rows).merge(sample_meta, on="sample_id", how="left", validate="many_to_one")
    scores_long.to_csv(out_dir / "tcga_rna_signature_scores_v0_2.csv", index=False)
    pd.DataFrame(recovery_rows).to_csv(out_dir / "tcga_rna_signature_gene_recovery_v0_2.csv", index=False)

    immune_path = root / "03_downloads/tcga_xena_pancanatlas/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz"
    anchors = read_tcga_immune_anchors(immune_path, kept_samples)
    anchors = anchors.merge(sample_meta, on="sample_id", how="left", validate="one_to_one")
    corr_rows: list[dict[str, Any]] = []
    merged = scores_long.merge(
        anchors.drop(columns=["patient_id", "cancer_hint", "primary_disease_hint", "sample_type_code", "sample_type_label"], errors="ignore"),
        on="sample_id",
        how="left",
        validate="many_to_one",
    )
    for (signature, cancer), group in merged.groupby(["signature", "cancer_hint"], dropna=False):
        if not safe_str(cancer):
            continue
        for var in TCGA_ANCHOR_VARS:
            if var not in group.columns:
                continue
            sub = group[["rna_signature_score", var]].dropna()
            if len(sub) < min_n or sub["rna_signature_score"].nunique() < 3 or sub[var].nunique() < 3:
                continue
            rho, p = spearmanr(sub["rna_signature_score"], sub[var])
            corr_rows.append(
                {
                    "signature": signature,
                    "cancer": cancer,
                    "tcga_anchor_variable": var,
                    "n": len(sub),
                    "spearman_rho": rho,
                    "p_value": p,
                    "p_family": "tcga_rna_signature_to_pancanatlas_immune_anchor_correlations_v0_2",
                    "interpretation_boundary": "RNA comparator only; not CPTAC validation and not sample-overlap analysis.",
                }
            )
    corr = pd.DataFrame(corr_rows)
    if not corr.empty:
        corr["fdr_bh"] = bh_fdr(corr["p_value"].to_numpy())
    corr.to_csv(out_dir / "tcga_rna_signature_context_correlations_v0_2.csv", index=False)
    return scores_long, corr, pd.DataFrame(recovery_rows)


def run_depmap_context(
    harmonization_dir: Path,
    qc_dir: Path,
    out_dir: Path,
) -> pd.DataFrame:
    inv = pd.read_csv(qc_dir / "depmap_inventory_v0_1.csv")
    model_map = pd.read_csv(harmonization_dir / "depmap_model_map_v0_1.csv")
    depmap_feature_map = pd.read_csv(harmonization_dir / "depmap_gene_feature_map_v0_1.csv")
    symbol_to_entrez, _entrez_to_symbol = depmap_symbol_entrez_map(depmap_feature_map)
    effect_path = Path(inv[inv["file_role"].astype(str) == "crispr_gene_effect"]["local_path"].iloc[0])
    dependency_path = Path(inv[inv["file_role"].astype(str) == "crispr_dependency"]["local_path"].iloc[0])
    effect_header = pd.read_csv(effect_path, nrows=0).columns.tolist()
    dep_header = pd.read_csv(dependency_path, nrows=0).columns.tolist()
    col_to_symbol = {}
    for col in effect_header[1:]:
        match = DEPMAP_GENE_RE.match(col)
        if match:
            col_to_symbol[col] = match.group("symbol")
    signature_symbols = {g for spec in SIGNATURES.values() for g in spec["genes"]}
    wanted_effect_cols = ["Unnamed: 0"] + [c for c, s in col_to_symbol.items() if s in signature_symbols]
    background_cols = effect_header[1:]
    effect = pd.read_csv(effect_path, usecols=wanted_effect_cols)
    dependency_cols = ["Unnamed: 0"] + [c for c in dep_header[1:] if c in wanted_effect_cols]
    dependency = pd.read_csv(dependency_path, usecols=dependency_cols)
    background_effect = pd.read_csv(effect_path, usecols=["Unnamed: 0"] + background_cols)
    model_lineage = model_map[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "AgeCategory"]].copy()
    effect = effect.rename(columns={"Unnamed: 0": "ModelID"}).merge(model_lineage, on="ModelID", how="left", validate="many_to_one")
    dependency = dependency.rename(columns={"Unnamed: 0": "ModelID"}).merge(model_lineage, on="ModelID", how="left", validate="many_to_one")
    background_effect = background_effect.rename(columns={"Unnamed: 0": "ModelID"})
    bg_matrix = background_effect[background_cols].apply(pd.to_numeric, errors="coerce")
    bg_gene_median = bg_matrix.median(axis=0, skipna=True)

    rows: list[dict[str, Any]] = []
    for signature, spec in SIGNATURES.items():
        genes = set(spec["genes"])
        cols = [c for c, s in col_to_symbol.items() if s in genes and c in effect.columns]
        if not cols:
            rows.append(
                {
                    "signature": signature,
                    "mapped_depmap_gene_count": 0,
                    "status": "no_depmap_gene_columns",
                    "allowed_claim": "No DepMap context available.",
                }
            )
            continue
        all_eff = effect[cols].apply(pd.to_numeric, errors="coerce")
        all_dep = dependency[[c for c in cols if c in dependency.columns]].apply(pd.to_numeric, errors="coerce")
        cns_mask = effect["OncotreeLineage"].astype(str).str.contains("CNS|Brain", case=False, na=False)
        gene_medians = all_eff.median(axis=0, skipna=True)
        gene_dep_prob = all_dep.mean(axis=0, skipna=True) if not all_dep.empty else pd.Series(dtype=float)
        background_less_equal = float(np.mean(bg_gene_median <= gene_medians.median())) if len(bg_gene_median) else np.nan
        try:
            mwu_p = mannwhitneyu(gene_medians.dropna(), bg_gene_median.dropna(), alternative="two-sided").pvalue
        except Exception:
            mwu_p = np.nan
        top_gene_effect = []
        for col, value in gene_medians.sort_values().head(8).items():
            top_gene_effect.append(f"{col_to_symbol.get(col, col)}:{value:.3f}")
        common_dep_genes = []
        for col, value in gene_dep_prob.sort_values(ascending=False).items():
            if np.isfinite(value) and value >= 0.5:
                common_dep_genes.append(f"{col_to_symbol.get(col, col)}:{value:.2f}")
        rows.append(
            {
                "signature": signature,
                "status": "ok",
                "requested_gene_count": len(genes),
                "mapped_depmap_gene_count": len(cols),
                "mapped_depmap_genes": ";".join(sorted({col_to_symbol.get(c, c) for c in cols})),
                "all_model_n": int(all_eff.shape[0]),
                "cns_brain_model_n": int(cns_mask.sum()),
                "median_gene_effect_all_models": float(np.nanmedian(all_eff.to_numpy(dtype=float))),
                "median_gene_effect_cns_brain_models": float(np.nanmedian(all_eff.loc[cns_mask].to_numpy(dtype=float))) if cns_mask.any() else np.nan,
                "median_gene_dependency_probability": float(np.nanmedian(all_dep.to_numpy(dtype=float))) if not all_dep.empty else np.nan,
                "gene_level_median_effect_vs_background_percentile_lower_is_more_dependent": background_less_equal * 100 if np.isfinite(background_less_equal) else np.nan,
                "gene_level_mannwhitney_p_vs_background": mwu_p,
                "top_more_dependent_signature_genes_by_median_effect": ";".join(top_gene_effect),
                "common_dependency_probability_ge_0_5_genes": ";".join(common_dep_genes[:10]),
                "common_dependency_gene_count_prob_ge_0_5": len(common_dep_genes),
                "allowed_claim": "Supportive cell-line dependency context only; not a therapeutic target or validation claim.",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty and "gene_level_mannwhitney_p_vs_background" in out.columns:
        out["gene_level_mannwhitney_fdr_bh_vs_background"] = bh_fdr(
            pd.to_numeric(out["gene_level_mannwhitney_p_vs_background"], errors="coerce").to_numpy()
        )
    out.to_csv(out_dir / "depmap_signature_dependency_context_v0_2.csv", index=False)
    return out


def build_platform_comparator(
    discovery_dir: Path,
    random_calibration: pd.DataFrame,
    tcga_corr: pd.DataFrame,
    depmap_context: pd.DataFrame,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cptac_corr = pd.read_csv(discovery_dir / "cptac_signature_context_correlations_v0_1.csv")
    residualized = pd.read_csv(discovery_dir / "cptac_signature_residualized_scores_v0_1.csv")
    rows: list[dict[str, Any]] = []
    for signature in SIGNATURES:
        cptac_immune = cptac_corr[
            (cptac_corr["signature"].astype(str) == signature)
            & (cptac_corr["context_variable"].astype(str) == "estimate_immune_score")
        ]
        cptac_purity = cptac_corr[
            (cptac_corr["signature"].astype(str) == signature)
            & (cptac_corr["context_variable"].astype(str) == "estimate_tumor_purity")
        ]
        rand = random_calibration[
            (random_calibration["signature"].astype(str) == signature)
            & (random_calibration["status"].astype(str) == "ok")
        ] if not random_calibration.empty else pd.DataFrame()
        anchors = SIGNATURE_ANCHORS.get(signature, [])
        tcga = tcga_corr[
            (tcga_corr["signature"].astype(str) == signature)
            & (tcga_corr["tcga_anchor_variable"].astype(str).isin(anchors))
        ] if not tcga_corr.empty else pd.DataFrame()
        dep = depmap_context[depmap_context["signature"].astype(str) == signature] if not depmap_context.empty else pd.DataFrame()
        rows.append(
            {
                "signature": signature,
                "cptac_immune_score_corr_strata": int(len(cptac_immune)),
                "cptac_median_immune_score_rho": float(cptac_immune["spearman_rho"].median()) if len(cptac_immune) else np.nan,
                "cptac_tumor_purity_corr_strata": int(len(cptac_purity)),
                "cptac_median_tumor_purity_rho": float(cptac_purity["spearman_rho"].median()) if len(cptac_purity) else np.nan,
                "cptac_residualized_strata": int(
                    residualized[residualized["signature"].astype(str) == signature][
                        ["source", "cancer", "datatype", "tumor_normal_hint"]
                    ].drop_duplicates().shape[0]
                ),
                "random_control_ok_tests": int(len(rand)),
                "random_control_empirical_p_lt_0_05_tests": int((rand["empirical_two_sided_p_abs_rho"] < 0.05).sum()) if len(rand) else 0,
                "random_control_median_percentile_abs_rho": float(rand["observed_abs_rho_percentile_vs_controls"].median()) if len(rand) else np.nan,
                "tcga_anchor_variables": ";".join(anchors),
                "tcga_anchor_corr_tests": int(len(tcga)),
                "tcga_median_anchor_rho": float(tcga["spearman_rho"].median()) if len(tcga) else np.nan,
                "tcga_positive_anchor_fraction": float((tcga["spearman_rho"] > 0).mean()) if len(tcga) else np.nan,
                "depmap_mapped_gene_count": int(dep["mapped_depmap_gene_count"].iloc[0]) if len(dep) and "mapped_depmap_gene_count" in dep else 0,
                "depmap_common_dependency_gene_count_prob_ge_0_5": int(dep["common_dependency_gene_count_prob_ge_0_5"].iloc[0]) if len(dep) and "common_dependency_gene_count_prob_ge_0_5" in dep else 0,
                "depmap_median_gene_effect_all_models": float(dep["median_gene_effect_all_models"].iloc[0]) if len(dep) and "median_gene_effect_all_models" in dep else np.nan,
            }
        )
    platform = pd.DataFrame(rows)
    platform.to_csv(out_dir / "platform_comparator_summary_v0_2.csv", index=False)

    claim_rows = []
    for _, row in platform.iterrows():
        signature = row["signature"]
        context_dominated = (
            abs(float(row["cptac_median_immune_score_rho"])) >= 0.35
            or abs(float(row["cptac_median_tumor_purity_rho"])) >= 0.35
        )
        random_supported = int(row["random_control_empirical_p_lt_0_05_tests"]) >= 4
        tcga_concordant = (
            int(row["tcga_anchor_corr_tests"]) >= 10
            and np.isfinite(row["tcga_median_anchor_rho"])
            and abs(float(row["tcga_median_anchor_rho"])) >= 0.25
        )
        residualized_ready = int(row["cptac_residualized_strata"]) >= 10
        if residualized_ready and random_supported and not tcga_concordant:
            tier = "possible_protein_or_phosphosite_prioritized_signal"
            allowed = "Prioritize for v0.3 phosphosite/kinase analysis; requires RNA comparator caveats and no validation wording."
        elif residualized_ready and tcga_concordant and context_dominated:
            tier = "rna_concordant_immune_context_signal"
            allowed = "Consistent immune-ecology program across CPTAC and TCGA RNA; likely context-dominated."
        elif residualized_ready and context_dominated:
            tier = "immune_purity_context_dominated_cptac_signal"
            allowed = "Use as immune/purity-context proteomic signal, not tumor-cell-intrinsic mechanism."
        elif residualized_ready:
            tier = "analysis_ready_low_context_signal"
            allowed = "Exploratory candidate for deeper v0.3 analysis."
        else:
            tier = "limited_support"
            allowed = "Keep in supplement or use only as negative/context result."
        claim_rows.append(
            {
                "signature": signature,
                "claim_tier_v0_2": tier,
                "allowed_claim": allowed,
                "forbidden_claim": "Do not write validated biomarker, causal mechanism, therapeutic target, or external validation.",
                "evidence_inputs": "CPTAC v0.1 residualization; feature-count random controls; TCGA RNA comparator; DepMap supportive context",
            }
        )
    claim = pd.DataFrame(claim_rows)
    claim.to_csv(out_dir / "discovery_claim_tier_table_v0_2.csv", index=False)
    return platform, claim


def write_figures(
    random_calibration: pd.DataFrame,
    tcga_corr: pd.DataFrame,
    platform: pd.DataFrame,
    depmap_context: pd.DataFrame,
    fig_dir: Path,
) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figures: list[str] = []
    ok = random_calibration[random_calibration["status"].astype(str) == "ok"].copy() if not random_calibration.empty else pd.DataFrame()
    if not ok.empty:
        plot = ok.groupby(["signature", "context_variable"])["observed_abs_rho_percentile_vs_controls"].median().reset_index()
        pivot = plot.pivot(index="signature", columns="context_variable", values="observed_abs_rho_percentile_vs_controls")
        pivot = pivot.rename(index=SIGNATURE_LABELS)
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        im = ax.imshow(pivot.fillna(0).to_numpy(), aspect="auto", cmap="magma", vmin=0, vmax=100)
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=25, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=8)
        ax.set_title("CPTAC observed context correlation percentile vs matched random controls")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Median percentile of abs(rho)")
        fig.tight_layout()
        out = fig_dir / "fig03_random_control_context_calibration_v0_2.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        figures.append(str(out))
    if not tcga_corr.empty:
        anchors = []
        for signature, group in tcga_corr.groupby("signature"):
            sub = group[group["tcga_anchor_variable"].isin(SIGNATURE_ANCHORS.get(signature, []))]
            anchors.append(sub)
        anchor_df = pd.concat(anchors, ignore_index=True) if anchors else pd.DataFrame()
        if not anchor_df.empty:
            plot = anchor_df.groupby(["signature", "tcga_anchor_variable"])["spearman_rho"].median().reset_index()
            pivot = plot.pivot(index="signature", columns="tcga_anchor_variable", values="spearman_rho")
            pivot = pivot.rename(index=SIGNATURE_LABELS)
            fig, ax = plt.subplots(figsize=(9.5, 4.8))
            im = ax.imshow(pivot.fillna(0).to_numpy(), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
            ax.set_xticks(range(pivot.shape[1]))
            ax.set_xticklabels(pivot.columns, rotation=35, ha="right", fontsize=8)
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels(pivot.index, fontsize=8)
            ax.set_title("TCGA RNA signature correlation with PanCanAtlas immune anchors")
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("Median Spearman rho across cancers")
            fig.tight_layout()
            out = fig_dir / "fig04_tcga_rna_anchor_correlation_heatmap_v0_2.png"
            fig.savefig(out, dpi=300)
            plt.close(fig)
            figures.append(str(out))
    if not platform.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        x = platform["cptac_median_immune_score_rho"].to_numpy(dtype=float)
        y = platform["tcga_median_anchor_rho"].to_numpy(dtype=float)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.scatter(x, y, s=55, color="#2c7fb8")
        for _, row in platform.iterrows():
            ax.annotate(
                SIGNATURE_LABELS.get(row["signature"], row["signature"]),
                (row["cptac_median_immune_score_rho"], row["tcga_median_anchor_rho"]),
                xytext=(5, 4),
                textcoords="offset points",
                fontsize=8,
            )
        ax.set_xlabel("CPTAC median rho with ESTIMATE immune score")
        ax.set_ylabel("TCGA RNA median rho with matched immune anchors")
        ax.set_title("Program-level CPTAC versus TCGA RNA comparator")
        fig.tight_layout()
        out = fig_dir / "fig05_cptac_tcga_program_comparator_v0_2.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        figures.append(str(out))
    if not depmap_context.empty and "median_gene_effect_all_models" in depmap_context.columns:
        plot = depmap_context.sort_values("median_gene_effect_all_models")
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.barh(plot["signature"].map(lambda x: SIGNATURE_LABELS.get(x, x)), plot["median_gene_effect_all_models"], color="#7fcdbb")
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Median DepMap CRISPR gene effect across signature genes")
        ax.set_title("DepMap supportive dependency context, lower is more dependent")
        fig.tight_layout()
        out = fig_dir / "fig06_depmap_signature_dependency_context_v0_2.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        figures.append(str(out))
    return figures


def write_report(
    report_path: Path,
    summary: dict[str, Any],
    platform: pd.DataFrame,
    claim: pd.DataFrame,
    random_calibration: pd.DataFrame,
    tcga_corr: pd.DataFrame,
    depmap_context: pd.DataFrame,
) -> None:
    lines = [
        "# Pan-cancer proteogenomic immune-evasion discovery v0.2",
        "",
        f"Generated: {summary['generated_at']}",
        "",
        "## Scope",
        "",
        "- Added feature-count-matched CPTAC random controls for immune/purity correlation calibration.",
        "- Added TCGA PanCanAtlas primary-tumor RNA comparator using the same curated signatures.",
        "- Added DepMap CRISPR gene-effect/dependency summaries as supportive context only.",
        "- CPTAC and TCGA were compared at signature/cancer/program level; no sample-level CPTAC-TCGA merge was performed.",
        "- No survival model, clinical prediction model, biomarker validation, causal mechanism, or therapeutic target claim was performed.",
        "",
        "## Output Counts",
        "",
    ]
    for key in [
        "random_control_rows",
        "random_controls_per_test_requested",
        "tcga_rna_score_rows",
        "tcga_rna_correlation_rows",
        "depmap_context_rows",
        "platform_summary_rows",
        "claim_tier_rows",
    ]:
        lines.append(f"- {key}: {summary.get(key)}")
    lines.extend(["", "## Key Program-Level Results", ""])
    if platform.empty:
        lines.append("- No platform comparator summary was generated.")
    else:
        for _, row in platform.sort_values("random_control_median_percentile_abs_rho", ascending=False).iterrows():
            lines.append(
                f"- `{row['signature']}`: CPTAC immune rho={row['cptac_median_immune_score_rho']:.3f}, "
                f"purity rho={row['cptac_median_tumor_purity_rho']:.3f}, "
                f"random-control percentile={row['random_control_median_percentile_abs_rho']:.1f}, "
                f"TCGA anchor rho={row['tcga_median_anchor_rho']:.3f}"
            )
    lines.extend(["", "## Claim Tiers", ""])
    if not claim.empty:
        for tier, n in claim["claim_tier_v0_2"].value_counts().items():
            lines.append(f"- `{tier}`: {n}")
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- v0.2 strengthens calibration and platform-context comparison, but it still does not validate a biomarker.",
            "- Strong immune/purity associations should be written as immune-ecology or TME-context signals.",
            "- DepMap findings are supportive cell-line dependency context only.",
            "- Any higher claim requires v0.3 phosphosite/kinase analysis and, ideally, external proteogenomic support.",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    start = time.time()
    args = parse_args()
    root = Path(args.project_root)
    harmonization_dir = resolve_path(args.harmonization_dir, root)
    qc_dir = resolve_path(args.qc_dir, root)
    discovery_dir = resolve_path(args.discovery_v0_1_dir, root)
    out_dir = resolve_path(args.out_dir, root)
    fig_dir = resolve_path(args.figure_dir, root)
    report_path = resolve_path(args.report, root)
    summary_path = resolve_path(args.summary, root)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    random_calibration = run_random_control_calibration(
        root=root,
        qc_dir=qc_dir,
        discovery_dir=discovery_dir,
        harmonization_dir=harmonization_dir,
        out_dir=out_dir,
        n_controls=args.n_random_controls,
        random_seed=args.random_seed,
        min_n=args.min_correlation_n,
        min_valid_feature_fraction=args.min_valid_feature_fraction,
    )
    tcga_scores, tcga_corr, tcga_recovery = run_tcga_rna_comparator(
        root=root,
        harmonization_dir=harmonization_dir,
        out_dir=out_dir,
        min_n=args.min_correlation_n,
    )
    depmap_context = run_depmap_context(harmonization_dir=harmonization_dir, qc_dir=qc_dir, out_dir=out_dir)
    platform_summary, claim = build_platform_comparator(
        discovery_dir=discovery_dir,
        random_calibration=random_calibration,
        tcga_corr=tcga_corr,
        depmap_context=depmap_context,
        out_dir=out_dir,
    )
    figures = write_figures(random_calibration, tcga_corr, platform_summary, depmap_context, fig_dir)

    summary = {
        "version": args.version,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "project_root": str(root),
        "random_control_rows": int(len(random_calibration)),
        "random_controls_per_test_requested": int(args.n_random_controls),
        "tcga_rna_score_rows": int(len(tcga_scores)),
        "tcga_rna_correlation_rows": int(len(tcga_corr)),
        "tcga_rna_recovery_rows": int(len(tcga_recovery)),
        "depmap_context_rows": int(len(depmap_context)),
        "platform_summary_rows": int(len(platform_summary)),
        "claim_tier_rows": int(len(claim)),
        "figures": [rel(Path(f), root) for f in figures],
        "outputs": {
            "random_control_calibration": rel(out_dir / "cptac_random_control_calibration_v0_2.csv", root),
            "tcga_rna_scores": rel(out_dir / "tcga_rna_signature_scores_v0_2.csv", root),
            "tcga_rna_correlations": rel(out_dir / "tcga_rna_signature_context_correlations_v0_2.csv", root),
            "platform_summary": rel(out_dir / "platform_comparator_summary_v0_2.csv", root),
            "depmap_context": rel(out_dir / "depmap_signature_dependency_context_v0_2.csv", root),
            "claim_tier": rel(out_dir / "discovery_claim_tier_table_v0_2.csv", root),
        },
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "elapsed_seconds": round(time.time() - start, 2),
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    write_report(report_path, summary, platform_summary, claim, random_calibration, tcga_corr, depmap_context)
    print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
