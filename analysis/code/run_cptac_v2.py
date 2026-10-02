#!/usr/bin/env python3
"""Locked BIO-05 v2 CPTAC cross-sectional concordance analysis.

All scores in this script are recomputed from the local BCM matrices. Historical
numeric outputs are never read.  The patient/case is the independent unit and
the primary protein/phosphoprotein estimand gives equal weight to genes after
within-gene feature collapse.
"""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", expand_legacy_paths("${SCI27_TMP_ROOT}/bio05_mpl_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, spearmanr


PROTOCOL_SHA256 = "f419f5ac589bf6d6ec61abb6663a47a5114b86a7d82980327b298adccb46aefc"
CANCERS = ["brca", "ccrcc", "coad", "gbm", "hnscc", "lscc", "luad", "ov", "pdac", "ucec"]
LAYERS = ["proteomics", "phosphoproteomics"]
PROGRAMS = [
    "antigen_presentation_mhc_i",
    "antigen_presentation_mhc_ii",
    "ifn_gamma_response",
    "cytolytic_t_cell_context",
    "checkpoint_exhaustion_context",
    "myeloid_inflammatory_context",
    "tgfb_emt_exclusion_context",
    "proteasome_antigen_processing",
]
PROGRAM_LABELS = {
    "antigen_presentation_mhc_i": "MHC-I",
    "antigen_presentation_mhc_ii": "MHC-II",
    "checkpoint_exhaustion_context": "Checkpoint",
    "cytolytic_t_cell_context": "Cytolytic",
    "ifn_gamma_response": "IFN-gamma",
    "myeloid_inflammatory_context": "Myeloid",
    "proteasome_antigen_processing": "Proteasome",
    "tgfb_emt_exclusion_context": "TGF/EMT",
}
CONTEXT_VARIABLES = [
    "estimate_tumor_purity",
    "estimate_immune_score",
    "cibersort_t_cell",
    "cibersort_myeloid",
]
MIN_GENES = 3
MIN_N = 20
Z95 = 1.959963984540054
PRIMARY_BOOTSTRAP_RESAMPLES = 20_000
PRIMARY_PERMUTATION_RESAMPLES = 999_999
PRIMARY_PERMUTATION_N_MAX = 30


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument(
        "--signature-table",
        default="",
        help="Optional signature CSV; defaults to DATA_ROOT/04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv.",
    )
    return parser.parse_args()


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if pd.isna(value):
        return None
    return value


def stable_ensg(value: Any) -> str:
    match = re.search(r"(ENSG\d+)(?:\.\d+)?", str(value))
    return match.group(1) if match else ""


def normalize_case_id(sample_id: Any) -> str:
    sid = str(sample_id).strip()
    if not sid or "QC" in sid.upper():
        return ""
    sid = re.sub(r"(_T|_A|_N)$", "", sid, flags=re.I)
    sid = re.sub(r"-(T|N)$", "", sid, flags=re.I)
    sid = re.sub(r"_D\d+$", "", sid, flags=re.I)
    return sid


def rna_tumor_hint(sample_id: Any) -> str:
    sid = str(sample_id).strip().upper()
    if re.search(r"(_T|-T)$", sid):
        return "tumor"
    if re.search(r"(_A|_N|-N)$", sid):
        return "adjacent_or_normal"
    return "unknown_tumor_like"


def tumor_priority(sample_id: Any) -> tuple[int, str]:
    sid = str(sample_id).strip()
    explicit = int(bool(re.search(r"(_T|-T)$", sid, flags=re.I)))
    return (-explicit, sid)


def select_unique_samples(sample_ids: Iterable[str], layer: str, cancer: str) -> tuple[list[int], list[dict[str, Any]]]:
    registry = pd.DataFrame({"sample_id": list(map(str, sample_ids))})
    registry["normalized_case_id"] = registry["sample_id"].map(normalize_case_id)
    registry["eligible"] = registry["normalized_case_id"].ne("")
    if layer == "transcriptomics":
        registry["tumor_hint"] = registry["sample_id"].map(rna_tumor_hint)
        registry["eligible"] &= registry["tumor_hint"].isin(["tumor", "unknown_tumor_like"])
    else:
        registry["tumor_hint"] = "tumor"
    keep: list[int] = []
    audit: list[dict[str, Any]] = []
    for case_id, group in registry[registry["eligible"]].groupby("normalized_case_id", sort=True):
        ordered = sorted(group.index.tolist(), key=lambda idx: tumor_priority(registry.loc[idx, "sample_id"]))
        keep.append(ordered[0])
        for rank, idx in enumerate(ordered, start=1):
            audit.append(
                {
                    "cancer": cancer,
                    "layer": layer,
                    "normalized_case_id": case_id,
                    "sample_id": registry.loc[idx, "sample_id"],
                    "candidate_count": len(ordered),
                    "priority_rank": rank,
                    "selected": int(rank == 1),
                    "decision_rule": "explicit_tumor_suffix_then_lexicographic_sample_id",
                }
            )
    return sorted(keep), audit


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def read_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with open_text(path) as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
    if header[0].strip().lower() in {"idx", "index", "gene", "gene_id"}:
        samples = header[1:]
    else:
        samples = header
    names = ["feature_id"] + samples
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        names=names,
        na_values=["NA", "NaN", "nan", ""],
        keep_default_na=True,
        low_memory=False,
    )
    values = frame[samples].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    return samples, frame["feature_id"].astype(str).tolist(), values


def zscore_rows(values: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.nanmean(values, axis=1, keepdims=True)
        std = np.nanstd(values, axis=1, ddof=1, keepdims=True)
        std[(std == 0) | ~np.isfinite(std)] = np.nan
        return (values - means) / std


def load_program_maps(signature_path: Path) -> tuple[pd.DataFrame, dict[str, dict[str, set[str]]], dict[str, int]]:
    signature_table = pd.read_csv(signature_path)
    observed = list(signature_table["signature"].drop_duplicates().astype(str))
    if set(observed) != set(PROGRAMS) or len(observed) != len(PROGRAMS):
        raise RuntimeError(f"Locked program order/identity mismatch: {observed}")
    maps: dict[str, dict[str, set[str]]] = {}
    for program, group in signature_table.groupby("signature", sort=False):
        gene_map: dict[str, set[str]] = defaultdict(set)
        for row in group.itertuples(index=False):
            gene = str(row.gene_symbol).strip().upper()
            raw_ids = getattr(row, "ensembl_stable_ids_from_cptac")
            if pd.isna(raw_ids):
                continue
            for token in str(raw_ids).split(";"):
                ensg = stable_ensg(token)
                if ensg:
                    gene_map[gene].add(ensg)
        maps[str(program)] = dict(gene_map)
    requested = signature_table.groupby("signature", sort=False)["gene_symbol"].nunique().astype(int).to_dict()
    return signature_table, maps, requested


def score_matrix(
    path: Path,
    cancer: str,
    layer: str,
    program_maps: dict[str, dict[str, set[str]]],
    requested_gene_counts: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    samples, feature_ids, values = read_matrix(path)
    keep_indices, duplicate_audit = select_unique_samples(samples, layer, cancer)
    samples_selected = [samples[i] for i in keep_indices]
    values = values[:, keep_indices]
    z = zscore_rows(values)
    stable_to_rows: dict[str, list[int]] = defaultdict(list)
    for idx, feature_id in enumerate(feature_ids):
        ensg = stable_ensg(feature_id)
        if ensg:
            stable_to_rows[ensg].append(idx)

    score_rows: list[dict[str, Any]] = []
    recovery_rows: list[dict[str, Any]] = []
    for program in PROGRAMS:
        gene_map = program_maps[program]
        gene_arrays: list[np.ndarray] = []
        gene_names: list[str] = []
        feature_rows: list[int] = []
        for gene in sorted(gene_map):
            rows = sorted({idx for ensg in gene_map[gene] for idx in stable_to_rows.get(ensg, [])})
            if not rows:
                continue
            with np.errstate(invalid="ignore"):
                gene_arrays.append(np.nanmean(z[rows, :], axis=0))
            gene_names.append(gene)
            feature_rows.extend(rows)
        requested_mapped = len(gene_map)
        requested_total = int(requested_gene_counts[program])
        recovered = len(gene_arrays)
        unique_feature_rows = sorted(set(feature_rows))
        recovery_rows.append(
            {
                "cancer": cancer,
                "layer": layer,
                "program": program,
                "program_label": PROGRAM_LABELS[program],
                "requested_gene_count": requested_total,
                "requested_gene_count_with_ensembl_mapping": requested_mapped,
                "recovered_gene_count": recovered,
                "recovered_feature_row_count": len(unique_feature_rows),
                "recovery_fraction_of_all_requested_genes": recovered / requested_total if requested_total else np.nan,
                "recovery_fraction_of_mapped_genes": recovered / requested_mapped if requested_mapped else np.nan,
                "recovered_genes": ";".join(gene_names),
                "source_file": str(path),
            }
        )
        if recovered == 0:
            continue
        gene_matrix = np.vstack(gene_arrays)
        gene_nonmissing = np.isfinite(gene_matrix).sum(axis=0)
        with np.errstate(invalid="ignore"):
            gene_scores = np.nanmean(gene_matrix, axis=0)
            feature_scores = np.nanmean(z[unique_feature_rows, :], axis=0)
        gene_scores[gene_nonmissing < MIN_GENES] = np.nan
        feature_scores[gene_nonmissing < MIN_GENES] = np.nan
        threshold1_scores = np.nanmean(gene_matrix, axis=0)
        threshold1_scores[gene_nonmissing < 1] = np.nan
        for col, sample_id in enumerate(samples_selected):
            score_rows.append(
                {
                    "cancer": cancer,
                    "layer": layer,
                    "program": program,
                    "program_label": PROGRAM_LABELS[program],
                    "sample_id": sample_id,
                    "normalized_case_id": normalize_case_id(sample_id),
                    "gene_equal_score": gene_scores[col],
                    "feature_row_equal_score": feature_scores[col],
                    "threshold1_gene_equal_score": threshold1_scores[col],
                    "nonmissing_gene_count": int(gene_nonmissing[col]),
                    "recovered_gene_count": recovered,
                    "recovered_feature_row_count": len(unique_feature_rows),
                }
            )
    sample_audit = pd.DataFrame(duplicate_audit)
    selected = pd.DataFrame(
        {
            "cancer": cancer,
            "layer": layer,
            "sample_id": samples_selected,
            "normalized_case_id": [normalize_case_id(x) for x in samples_selected],
        }
    )
    return pd.DataFrame(score_rows), pd.DataFrame(recovery_rows), sample_audit, selected


def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out


def spearman_ci(rho: float, n: int) -> tuple[float, float]:
    if not np.isfinite(rho) or n <= 3:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / math.sqrt(n - 3)
    return float(np.tanh(z - Z95 * se)), float(np.tanh(z + Z95 * se))


def stable_seed(*parts: str) -> int:
    token = "::".join(map(str, parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(token).digest()[:8], "little")


def rowwise_spearman(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Vectorized Spearman rho for bootstrap rows, including induced ties."""
    rx = rankdata(x, axis=1, method="average")
    ry = rankdata(y, axis=1, method="average")
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    numerator = np.sum(rx * ry, axis=1)
    denominator = np.sqrt(np.sum(rx**2, axis=1) * np.sum(ry**2, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return numerator / denominator


def spearman_bca_ci(
    x: np.ndarray,
    y: np.ndarray,
    rho: float,
    seed: int,
    n_resamples: int = PRIMARY_BOOTSTRAP_RESAMPLES,
) -> tuple[float, float, int, str]:
    """Paired-case BCa bootstrap interval for Spearman rho."""
    n = len(x)
    rng = np.random.default_rng(seed)
    boot = np.empty(n_resamples, dtype=float)
    batch = 2_000
    for start in range(0, n_resamples, batch):
        stop = min(start + batch, n_resamples)
        indices = rng.integers(0, n, size=(stop - start, n))
        boot[start:stop] = rowwise_spearman(x[indices], y[indices])
    boot = boot[np.isfinite(boot)]
    if len(boot) < int(0.99 * n_resamples):
        raise RuntimeError(f"Too few finite Spearman bootstrap replicates: {len(boot)}/{n_resamples}")

    less = np.sum(boot < rho)
    equal = np.sum(boot == rho)
    proportion = (less + 0.5 * equal) / len(boot)
    proportion = float(np.clip(proportion, 0.5 / len(boot), 1 - 0.5 / len(boot)))
    z0 = float(norm.ppf(proportion))

    jack = np.array(
        [spearmanr(np.delete(x, idx), np.delete(y, idx)).statistic for idx in range(n)],
        dtype=float,
    )
    jack_mean = float(np.mean(jack))
    centered = jack_mean - jack
    denominator = 6.0 * float(np.sum(centered**2) ** 1.5)
    acceleration = float(np.sum(centered**3) / denominator) if denominator > 0 else 0.0

    adjusted = []
    for alpha in (0.025, 0.975):
        za = float(norm.ppf(alpha))
        denom = 1.0 - acceleration * (z0 + za)
        adjusted.append(float(norm.cdf(z0 + (z0 + za) / denom)))
    adjusted = np.clip(adjusted, 0, 1)
    low, high = np.quantile(boot, adjusted)
    return float(low), float(high), int(len(boot)), "paired_case_bca_bootstrap"


def spearman_permutation_p(
    x: np.ndarray,
    y: np.ndarray,
    rho: float,
    seed: int,
    n_resamples: int = PRIMARY_PERMUTATION_RESAMPLES,
) -> tuple[float, float]:
    """Fixed-seed Monte Carlo two-sided permutation P value."""
    rng = np.random.default_rng(seed)
    rx = rankdata(x, method="average")
    ry = rankdata(y, method="average")
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denominator = float(np.sqrt(np.sum(rx**2) * np.sum(ry**2)))
    exceed = 0
    batch = 25_000
    for start in range(0, n_resamples, batch):
        size = min(batch, n_resamples - start)
        random_order = np.argsort(rng.random((size, len(ry))), axis=1)
        permuted = ry[random_order]
        permuted_rho = (permuted @ rx) / denominator
        exceed += int(np.sum(np.abs(permuted_rho) >= abs(rho) - 1e-15))
    p_value = (exceed + 1.0) / (n_resamples + 1.0)
    mcse = math.sqrt(p_value * (1 - p_value) / (n_resamples + 1.0))
    return float(p_value), float(mcse)


def leave_one_out_influence(x: np.ndarray, y: np.ndarray, rho: float) -> tuple[float, str]:
    if len(x) < 5 or not np.isfinite(rho):
        return np.nan, ""
    deltas: list[float] = []
    for idx in range(len(x)):
        test_rho = spearmanr(np.delete(x, idx), np.delete(y, idx)).statistic
        deltas.append(abs(float(test_rho) - rho) if np.isfinite(test_rho) else np.nan)
    max_idx = int(np.nanargmax(deltas))
    return float(deltas[max_idx]), str(max_idx)


def evaluate_pair(
    x_raw: pd.Series,
    y_raw: pd.Series,
    case_ids: pd.Series,
    *,
    calibrated_primary: bool = False,
    seed_key: str = "",
) -> dict[str, Any]:
    x = pd.to_numeric(x_raw, errors="coerce")
    y = pd.to_numeric(y_raw, errors="coerce")
    mask = x.notna() & y.notna()
    xv = x[mask].to_numpy(dtype=float)
    yv = y[mask].to_numpy(dtype=float)
    cases = case_ids[mask].astype(str).to_numpy()
    n = len(xv)
    x_distinct = int(pd.Series(xv).nunique())
    y_distinct = int(pd.Series(yv).nunique())
    base = {
        "eligible_pair_rows": int(len(x_raw)),
        "complete_pair_n": n,
        "missing_pair_n": int(len(x_raw) - n),
        "x_distinct_n": x_distinct,
        "y_distinct_n": y_distinct,
        "x_tie_fraction": 1 - x_distinct / n if n else np.nan,
        "y_tie_fraction": 1 - y_distinct / n if n else np.nan,
        "spearman_rho": np.nan,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p_value": np.nan,
        "p_value_asymptotic": np.nan,
        "p_value_method": "",
        "permutation_resamples": 0,
        "permutation_p_mcse": np.nan,
        "ci_method": "",
        "ci_resamples": 0,
        "max_abs_leave_one_out_rho_delta": np.nan,
        "most_influential_case_id": "",
    }
    if n < MIN_N:
        base.update(status="not_evaluable", reason="complete_pair_n_below_20")
        return base
    if x_distinct < 3 or y_distinct < 3:
        base.update(status="not_evaluable", reason="fewer_than_3_distinct_values")
        return base
    result = spearmanr(xv, yv)
    rho, asymptotic_p = float(result.statistic), float(result.pvalue)
    if calibrated_primary:
        ci_low, ci_high, ci_resamples, ci_method = spearman_bca_ci(
            xv, yv, rho, stable_seed("BIO-05", "primary-ci", seed_key)
        )
        if n <= PRIMARY_PERMUTATION_N_MAX or base["x_tie_fraction"] > 0 or base["y_tie_fraction"] > 0:
            p_value, p_mcse = spearman_permutation_p(
                xv, yv, rho, stable_seed("BIO-05", "primary-p", seed_key)
            )
            p_method = "two_sided_monte_carlo_permutation"
            permutation_resamples = PRIMARY_PERMUTATION_RESAMPLES
        else:
            p_value, p_mcse = asymptotic_p, np.nan
            p_method = "scipy_spearman_asymptotic_no_ties_n_gt_30"
            permutation_resamples = 0
    else:
        ci_low, ci_high = spearman_ci(rho, n)
        ci_resamples, ci_method = 0, "fisher_z_approximation_sensitivity_only"
        p_value, p_mcse = asymptotic_p, np.nan
        p_method, permutation_resamples = "scipy_spearman_asymptotic", 0
    influence, influence_idx = leave_one_out_influence(xv, yv, rho)
    base.update(
        status="evaluable",
        reason="",
        spearman_rho=rho,
        ci_lower=ci_low,
        ci_upper=ci_high,
        p_value=p_value,
        p_value_asymptotic=asymptotic_p,
        p_value_method=p_method,
        permutation_resamples=permutation_resamples,
        permutation_p_mcse=p_mcse,
        ci_method=ci_method,
        ci_resamples=ci_resamples,
        max_abs_leave_one_out_rho_delta=influence,
        most_influential_case_id=cases[int(influence_idx)] if influence_idx else "",
    )
    return base


def make_matched(scores: pd.DataFrame, score_col: str) -> pd.DataFrame:
    rna = scores[scores["layer"].eq("transcriptomics")][
        ["cancer", "program", "program_label", "sample_id", "normalized_case_id", score_col, "nonmissing_gene_count"]
    ].rename(columns={"sample_id": "rna_sample_id", score_col: "rna_score", "nonmissing_gene_count": "rna_nonmissing_genes"})
    frames: list[pd.DataFrame] = []
    for layer in LAYERS:
        omic = scores[scores["layer"].eq(layer)][
            ["cancer", "program", "sample_id", "normalized_case_id", score_col, "nonmissing_gene_count"]
        ].rename(columns={"sample_id": "omic_sample_id", score_col: "omic_score", "nonmissing_gene_count": "omic_nonmissing_genes"})
        keys = ["cancer", "program", "normalized_case_id"]
        if rna.duplicated(keys).any() or omic.duplicated(keys).any():
            raise RuntimeError(f"Uniqueness assertion failed before {layer} matching")
        merged = rna.merge(omic, on=keys, how="inner", validate="one_to_one")
        merged["comparison_layer"] = layer
        frames.append(merged)
    out = pd.concat(frames, ignore_index=True)
    key = ["cancer", "program", "comparison_layer", "normalized_case_id"]
    if out.duplicated(key).any():
        raise RuntimeError("Matched table violates one-case-per-registered-stratum assertion")
    return out


def build_registry(matched: pd.DataFrame, family: str, *, calibrated_primary: bool = False) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cancer in CANCERS:
        for program in PROGRAMS:
            for layer in LAYERS:
                group = matched[
                    matched["cancer"].eq(cancer)
                    & matched["program"].eq(program)
                    & matched["comparison_layer"].eq(layer)
                ]
                result = evaluate_pair(
                    group["rna_score"],
                    group["omic_score"],
                    group["normalized_case_id"],
                    calibrated_primary=calibrated_primary,
                    seed_key=f"{cancer}::{program}::{layer}",
                )
                result.update(
                    comparison=f"rna_vs_{layer}",
                    comparison_scope="within_cancer",
                    cancer=cancer,
                    program=program,
                    signature=program,
                    program_label=PROGRAM_LABELS[program],
                    signature_label=PROGRAM_LABELS[program],
                    comparison_layer=layer,
                    fdr_family=family,
                )
                rows.append(result)
    registry = pd.DataFrame(rows)
    registry["q_value_bh"] = bh_fdr(registry["p_value"])
    registry["fdr_bh"] = registry["q_value_bh"]
    return registry


def build_heterogeneity(registry: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    loo_rows: list[dict[str, Any]] = []
    valid = registry[registry["status"].eq("evaluable")].copy()
    for (layer, program), group in valid.groupby(["comparison_layer", "program"], sort=False):
        rho = group["spearman_rho"].to_numpy(dtype=float)
        n = group["complete_pair_n"].to_numpy(dtype=float)
        z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
        w = np.maximum(n - 3, 1)
        fixed_z = np.sum(w * z) / np.sum(w)
        q_stat = float(np.sum(w * (z - fixed_z) ** 2))
        df = max(len(group) - 1, 0)
        i2 = max(0.0, (q_stat - df) / q_stat) if q_stat > 0 else 0.0
        rows.append(
            {
                "comparison_layer": layer,
                "program": program,
                "program_label": PROGRAM_LABELS[program],
                "evaluable_cancer_n": len(group),
                "total_complete_pairs_across_strata": int(n.sum()),
                "median_within_cancer_rho": float(np.median(rho)),
                "q25_within_cancer_rho": float(np.quantile(rho, 0.25)),
                "q75_within_cancer_rho": float(np.quantile(rho, 0.75)),
                "minimum_within_cancer_rho": float(np.min(rho)),
                "maximum_within_cancer_rho": float(np.max(rho)),
                "positive_cancer_n": int((rho > 0).sum()),
                "positive_cancer_fraction": float((rho > 0).mean()),
                "descriptive_cochran_q_fisher_z": q_stat,
                "descriptive_i2": i2,
                "interpretation_boundary": "Cancer-specific effects are primary; Q and I2 are descriptive and no single pan-cancer effect is an estimand.",
            }
        )
        for excluded in group["cancer"]:
            retained = group[~group["cancer"].eq(excluded)]["spearman_rho"].to_numpy(dtype=float)
            loo_rows.append(
                {
                    "comparison_layer": layer,
                    "program": program,
                    "program_label": PROGRAM_LABELS[program],
                    "excluded_cancer": excluded,
                    "retained_cancer_n": len(retained),
                    "leave_one_cancer_out_median_rho": float(np.median(retained)),
                    "leave_one_cancer_out_positive_fraction": float((retained > 0).mean()),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(loo_rows)


def build_context(data_root: Path, inventory: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    file_rows: list[dict[str, Any]] = []
    subset = inventory[inventory["source"].astype(str).eq("washu")]
    for item in subset.itertuples(index=False):
        path = data_root / str(item.local_path)
        if not path.exists():
            raise FileNotFoundError(path)
        file_rows.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path), "role": str(item.datatype)})
        datatype = str(item.datatype)
        cancer = str(item.cancer)
        if datatype == "tumor_purity":
            raw = pd.read_csv(path, sep="\t")
            frame = pd.DataFrame(
                {
                    "sample_id": raw["Sample_ID"].astype(str),
                    "estimate_tumor_purity": pd.to_numeric(raw["TumorPurity"], errors="coerce"),
                    "estimate_immune_score": pd.to_numeric(raw["ImmuneScore"], errors="coerce"),
                }
            )
        elif datatype == "cibersort":
            raw = pd.read_csv(path, sep="\t")
            numeric = raw.copy()
            for col in numeric.columns:
                if col != "Input Sample":
                    numeric[col] = pd.to_numeric(numeric[col], errors="coerce")
            t_cols = [c for c in numeric if c.startswith("T cells")]
            myeloid_cols = [c for c in numeric if c.startswith(("Macrophages", "Monocytes", "Dendritic cells"))]
            frame = pd.DataFrame(
                {
                    "sample_id": raw["Input Sample"].astype(str),
                    "cibersort_t_cell": numeric[t_cols].sum(axis=1, min_count=1),
                    "cibersort_myeloid": numeric[myeloid_cols].sum(axis=1, min_count=1),
                }
            )
        else:
            continue
        frame["normalized_case_id"] = frame["sample_id"].map(normalize_case_id)
        frame["context_cancer"] = cancer
        frame["context_source"] = datatype
        frames.append(frame)
    stacked = pd.concat(frames, ignore_index=True, sort=False)
    numeric_cols = [c for c in CONTEXT_VARIABLES if c in stacked]
    context = stacked.groupby("normalized_case_id", as_index=False)[numeric_cols].mean(numeric_only=True)
    if context.duplicated("normalized_case_id").any():
        raise RuntimeError("Context case map is not unique")
    return context, pd.DataFrame(file_rows)


def residualize(matched: pd.DataFrame, context: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    residual_rows: list[dict[str, Any]] = []
    diagnostic_rows: list[dict[str, Any]] = []
    for keys, group in matched.groupby(["cancer", "program", "comparison_layer"], sort=False):
        cancer, program, layer = keys
        work = group[["normalized_case_id", "rna_score", "omic_score"]].dropna().copy()
        n = len(work)
        status = "evaluable" if n >= MIN_N and work["rna_score"].nunique() >= 3 and work["omic_score"].nunique() >= 3 else "not_evaluable"
        diag: dict[str, Any] = {
            "cancer": cancer,
            "program": program,
            "program_label": PROGRAM_LABELS[program],
            "comparison_layer": layer,
            "model_n": n,
            "status": status,
        }
        if status != "evaluable":
            diag["reason"] = "model_n_or_distinct_value_gate_failed"
            diagnostic_rows.append(diag)
            continue
        x = work["rna_score"].to_numpy(dtype=float)
        y = work["omic_score"].to_numpy(dtype=float)
        design = np.column_stack([np.ones(n), x])
        beta = np.linalg.lstsq(design, y, rcond=None)[0]
        fitted = design @ beta
        residual = y - fitted
        ss_total = float(np.sum((y - np.mean(y)) ** 2))
        r2 = 1 - float(np.sum(residual**2)) / ss_total if ss_total > 0 else np.nan
        rx, ry = rankdata(x, method="average"), rankdata(y, method="average")
        rank_design = np.column_stack([np.ones(n), rx])
        rank_beta = np.linalg.lstsq(rank_design, ry, rcond=None)[0]
        rank_residual = ry - rank_design @ rank_beta
        quadratic = np.column_stack([np.ones(n), x, x**2])
        quad_beta = np.linalg.lstsq(quadratic, y, rcond=None)[0]
        quad_residual = y - quadratic @ quad_beta
        quad_r2 = 1 - float(np.sum(quad_residual**2)) / ss_total if ss_total > 0 else np.nan
        diag.update(
            reason="",
            intercept=float(beta[0]),
            rna_beta=float(beta[1]),
            linear_r2=r2,
            quadratic_r2=quad_r2,
            quadratic_incremental_r2=quad_r2 - r2 if np.isfinite(quad_r2) and np.isfinite(r2) else np.nan,
            residual_mean=float(np.mean(residual)),
            residual_sd=float(np.std(residual, ddof=1)),
            residual_skew=float(pd.Series(residual).skew()),
            residual_excess_kurtosis=float(pd.Series(residual).kurt()),
            residual_vs_rna_spearman=float(spearmanr(residual, x).statistic),
            rank_adjustment_r2=float(1 - np.sum(rank_residual**2) / np.sum((ry - np.mean(ry)) ** 2)),
        )
        diagnostic_rows.append(diag)
        for idx, row in enumerate(work.itertuples(index=False)):
            residual_rows.append(
                {
                    "cancer": cancer,
                    "program": program,
                    "program_label": PROGRAM_LABELS[program],
                    "comparison_layer": layer,
                    "normalized_case_id": row.normalized_case_id,
                    "rna_score": row.rna_score,
                    "omic_score": row.omic_score,
                    "linear_residual": float(residual[idx]),
                    "rank_residual": float(rank_residual[idx]),
                }
            )
    residuals = pd.DataFrame(residual_rows)
    diagnostics = pd.DataFrame(diagnostic_rows)
    joined = residuals.merge(context, on="normalized_case_id", how="left", validate="many_to_one")
    context_rows: list[dict[str, Any]] = []
    for cancer in CANCERS:
        for program in PROGRAMS:
            for layer in LAYERS:
                group = joined[
                    joined["cancer"].eq(cancer)
                    & joined["program"].eq(program)
                    & joined["comparison_layer"].eq(layer)
                ]
                for variable in CONTEXT_VARIABLES:
                    linear = evaluate_pair(group["linear_residual"], group.get(variable, pd.Series(index=group.index, dtype=float)), group["normalized_case_id"])
                    rank = evaluate_pair(group["rank_residual"], group.get(variable, pd.Series(index=group.index, dtype=float)), group["normalized_case_id"])
                    context_rows.append(
                        {
                            "cancer": cancer,
                            "program": program,
                            "program_label": PROGRAM_LABELS[program],
                            "comparison_layer": layer,
                            "context_variable": variable,
                            "linear_status": linear["status"],
                            "linear_reason": linear["reason"],
                            "linear_n": linear["complete_pair_n"],
                            "linear_rho": linear["spearman_rho"],
                            "linear_p_value": linear["p_value"],
                            "rank_status": rank["status"],
                            "rank_reason": rank["reason"],
                            "rank_n": rank["complete_pair_n"],
                            "rank_rho": rank["spearman_rho"],
                            "rank_p_value": rank["p_value"],
                        }
                    )
    context_registry = pd.DataFrame(context_rows)
    context_registry["linear_q_value_bh"] = bh_fdr(context_registry["linear_p_value"])
    context_registry["rank_q_value_bh"] = bh_fdr(context_registry["rank_p_value"])
    return residuals, diagnostics, context_registry


def standardized_pooled(matched: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (layer, program), group in matched.groupby(["comparison_layer", "program"], sort=False):
        frames: list[pd.DataFrame] = []
        for cancer, cg in group.groupby("cancer"):
            comp = cg[["rna_score", "omic_score"]].dropna().copy()
            if len(comp) < MIN_N:
                continue
            for col in ["rna_score", "omic_score"]:
                sd = comp[col].std(ddof=1)
                comp[col] = (comp[col] - comp[col].mean()) / sd if np.isfinite(sd) and sd > 0 else np.nan
            comp["cancer"] = cancer
            frames.append(comp.dropna())
        pool = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        if len(pool) >= MIN_N:
            result = spearmanr(pool["rna_score"], pool["omic_score"])
            rho, p = float(result.statistic), float(result.pvalue)
        else:
            rho, p = np.nan, np.nan
        rows.append(
            {
                "comparison_layer": layer,
                "program": program,
                "program_label": PROGRAM_LABELS[program],
                "contributing_cancer_n": int(pool["cancer"].nunique()) if not pool.empty else 0,
                "pooled_complete_n": len(pool),
                "within_cancer_standardized_pooled_spearman_rho": rho,
                "descriptive_p_value": p,
                "interpretation_boundary": "Descriptive standardized pooled sensitivity; cancer-specific estimates remain primary.",
            }
        )
    return pd.DataFrame(rows)


def make_figures(registry: pd.DataFrame, sensitivity: pd.DataFrame, flow: pd.DataFrame, figure_dir: Path) -> list[str]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    made: list[str] = []
    plot = registry.copy()
    plot["stratum"] = plot["cancer"].str.upper() + " " + plot["comparison_layer"].map({"proteomics": "P", "phosphoproteomics": "PP"})
    pivot = plot.pivot(index="program_label", columns="stratum", values="spearman_rho").reindex([PROGRAM_LABELS[p] for p in PROGRAMS])
    fig, ax = plt.subplots(figsize=(12, 5.2))
    image = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(range(len(pivot.index)), pivot.index, fontsize=9)
    ax.set_title("CPTAC within-cancer RNA–protein/phosphoprotein concordance")
    fig.colorbar(image, ax=ax, label="Spearman rho", fraction=0.03, pad=0.02)
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        out = figure_dir / f"fig01_primary_concordance_heatmap.{ext}"
        fig.savefig(out, dpi=300)
        made.append(str(out))
    plt.close(fig)

    comp = registry[["cancer", "program", "comparison_layer", "spearman_rho"]].merge(
        sensitivity[["cancer", "program", "comparison_layer", "spearman_rho"]],
        on=["cancer", "program", "comparison_layer"],
        suffixes=("_gene_equal", "_feature_equal"),
        validate="one_to_one",
    ).dropna()
    fig, ax = plt.subplots(figsize=(6, 6))
    for layer, group in comp.groupby("comparison_layer"):
        ax.scatter(group["spearman_rho_gene_equal"], group["spearman_rho_feature_equal"], s=20, alpha=0.7, label=layer)
    ax.axline((-1, -1), (1, 1), color="black", linewidth=0.8, linestyle="--")
    ax.set(xlim=(-1, 1), ylim=(-1, 1), xlabel="Gene-equal primary rho", ylabel="Feature-row sensitivity rho", title="Primary versus feature-row scoring")
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        out = figure_dir / f"fig02_gene_vs_feature_sensitivity.{ext}"
        fig.savefig(out, dpi=300)
        made.append(str(out))
    plt.close(fig)

    counts = flow.groupby(["cancer", "comparison_layer"], as_index=False)["matched_unique_cases"].max()
    pivot_n = counts.pivot(index="cancer", columns="comparison_layer", values="matched_unique_cases").reindex(CANCERS)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    pivot_n.plot(kind="bar", ax=ax)
    ax.set(xlabel="Cancer", ylabel="Maximum matched unique cases", title="Matched CPTAC case denominators")
    ax.legend(title="Layer", fontsize=8)
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        out = figure_dir / f"fig03_matched_case_denominators.{ext}"
        fig.savefig(out, dpi=300)
        made.append(str(out))
    plt.close(fig)
    return made


def main() -> int:
    args = parse_args()
    started = time.time()
    analysis_root = Path(args.analysis_root).resolve()
    data_root = Path(args.data_root).resolve()
    protocol_path = Path(args.protocol).resolve()
    if sha256_file(protocol_path) != PROTOCOL_SHA256:
        raise RuntimeError("Locked protocol SHA-256 mismatch")
    results_dir = analysis_root / "results" / "cptac"
    tables_dir = analysis_root / "tables"
    figures_dir = analysis_root / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    signature_path = (
        Path(args.signature_table).expanduser().resolve()
        if args.signature_table
        else data_root / "04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv"
    )
    inventory_path = data_root / "04_processed/first_pass_qc/cptac_processed_inventory_v0_1.csv"
    rna_manifest_path = data_root / "04_processed/discovery_v0_4/cptac_bcm_transcriptomics_download_manifest_v0_4.csv"
    sample_map_path = data_root / "04_processed/harmonization_v0_1/cptac_sample_map_v0_1.csv"
    for path in [signature_path, inventory_path, rna_manifest_path, sample_map_path]:
        if not path.exists():
            raise FileNotFoundError(path)
    _, program_maps, requested_gene_counts = load_program_maps(signature_path)
    inventory = pd.read_csv(inventory_path)
    rna_manifest = pd.read_csv(rna_manifest_path)

    matrix_specs: list[tuple[str, str, Path]] = []
    for cancer in CANCERS:
        row = rna_manifest[rna_manifest["cancer"].astype(str).eq(cancer)]
        if len(row) != 1:
            raise RuntimeError(f"RNA manifest must contain exactly one row for {cancer}")
        matrix_specs.append((cancer, "transcriptomics", data_root / str(row.iloc[0]["project_local_path"])))
        for layer in LAYERS:
            hit = inventory[
                inventory["source"].astype(str).eq("bcm")
                & inventory["cancer"].astype(str).eq(cancer)
                & inventory["datatype"].astype(str).eq(layer)
                & inventory["tumor_normal_hint"].astype(str).eq("tumor")
            ]
            if len(hit) != 1:
                raise RuntimeError(f"Expected one BCM tumor {layer} matrix for {cancer}; observed {len(hit)}")
            matrix_specs.append((cancer, layer, data_root / str(hit.iloc[0]["local_path"])))

    score_frames: list[pd.DataFrame] = []
    recovery_frames: list[pd.DataFrame] = []
    duplicate_frames: list[pd.DataFrame] = []
    sample_frames: list[pd.DataFrame] = []
    input_rows: list[dict[str, Any]] = []
    for cancer, layer, path in matrix_specs:
        if not path.exists():
            raise FileNotFoundError(path)
        score, recovery, duplicate, selected = score_matrix(path, cancer, layer, program_maps, requested_gene_counts)
        score_frames.append(score)
        recovery_frames.append(recovery)
        duplicate_frames.append(duplicate)
        sample_frames.append(selected)
        input_rows.append(
            {
                "cancer": cancer,
                "layer": layer,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    scores = pd.concat(score_frames, ignore_index=True)
    recovery = pd.concat(recovery_frames, ignore_index=True)
    duplicate_audit = pd.concat(duplicate_frames, ignore_index=True)
    selected_samples = pd.concat(sample_frames, ignore_index=True)

    sample_map = pd.read_csv(sample_map_path)
    map_check_rows: list[dict[str, Any]] = []
    for (cancer, layer), group in selected_samples[selected_samples["layer"].isin(LAYERS)].groupby(["cancer", "layer"]):
        ref = sample_map[
            sample_map["source"].astype(str).eq("bcm")
            & sample_map["cancer"].astype(str).eq(cancer)
            & sample_map["datatype"].astype(str).eq(layer)
            & sample_map["tumor_normal_hint"].astype(str).eq("tumor")
        ]
        raw_pairs = set(zip(group["sample_id"].astype(str), group["normalized_case_id"].astype(str)))
        ref_pairs = set(zip(ref["sample_id"].astype(str), ref["inferred_case_id"].map(normalize_case_id).astype(str)))
        map_check_rows.append(
            {
                "cancer": cancer,
                "layer": layer,
                "raw_selected_sample_n": len(raw_pairs),
                "reference_mapping_sample_n": len(ref_pairs),
                "raw_pairs_absent_from_reference_n": len(raw_pairs - ref_pairs),
                "reference_pairs_absent_from_raw_n": len(ref_pairs - raw_pairs),
                "status": "pass" if raw_pairs == ref_pairs else "fail",
            }
        )
    mapping_audit = pd.DataFrame(map_check_rows)
    if mapping_audit["status"].ne("pass").any():
        raise RuntimeError("Raw tumor sample headers do not reproduce the frozen sample-to-case map")

    scores.to_csv(results_dir / "cptac_program_scores_v2.csv", index=False)
    recovery.to_csv(results_dir / "cptac_gene_recovery_v2.csv", index=False)
    duplicate_audit.to_csv(results_dir / "cptac_duplicate_priority_audit_v2.csv", index=False)
    selected_samples.to_csv(results_dir / "cptac_selected_sample_map_v2.csv", index=False)
    mapping_audit.to_csv(results_dir / "cptac_sample_mapping_crosscheck_v2.csv", index=False)
    pd.DataFrame(input_rows).to_csv(results_dir / "cptac_input_files_v2.csv", index=False)

    primary_matched = make_matched(scores, "gene_equal_score")
    feature_matched = make_matched(scores, "feature_row_equal_score")
    threshold1_matched = make_matched(scores, "threshold1_gene_equal_score")
    primary_matched.to_csv(results_dir / "cptac_matched_gene_equal_scores_v2.csv", index=False)

    primary = build_registry(
        primary_matched,
        "cptac_primary_10_cancers_x_8_programs_x_2_layers",
        calibrated_primary=True,
    )
    feature = build_registry(feature_matched, "cptac_feature_row_sensitivity_10x8x2")
    threshold1 = build_registry(threshold1_matched, "cptac_min_1_gene_sensitivity_10x8x2")
    primary.to_csv(results_dir / "cptac_primary_concordance_registry_v2.csv", index=False)
    feature.to_csv(results_dir / "cptac_feature_row_sensitivity_registry_v2.csv", index=False)
    threshold1.to_csv(results_dir / "cptac_threshold1_sensitivity_registry_v2.csv", index=False)

    heterogeneity, leave_one_cancer = build_heterogeneity(primary)
    standardized = standardized_pooled(primary_matched)
    heterogeneity.to_csv(results_dir / "cptac_heterogeneity_summary_v2.csv", index=False)
    leave_one_cancer.to_csv(results_dir / "cptac_leave_one_cancer_out_v2.csv", index=False)
    standardized.to_csv(results_dir / "cptac_within_cancer_standardized_pooled_sensitivity_v2.csv", index=False)

    context, context_inputs = build_context(data_root, inventory)
    context.to_csv(results_dir / "cptac_context_features_recomputed_v2.csv", index=False)
    context_inputs.to_csv(results_dir / "cptac_context_input_files_v2.csv", index=False)
    residual_scores, residual_diagnostics, residual_context = residualize(primary_matched, context)
    residual_scores.to_csv(results_dir / "cptac_residual_scores_v2.csv", index=False)
    residual_diagnostics.to_csv(results_dir / "cptac_residual_model_diagnostics_v2.csv", index=False)
    residual_context.to_csv(results_dir / "cptac_residual_context_registry_v2.csv", index=False)

    flow = primary[
        ["cancer", "program", "program_label", "comparison_layer", "eligible_pair_rows", "complete_pair_n", "missing_pair_n", "status", "reason"]
    ].rename(columns={"eligible_pair_rows": "matched_unique_cases"})
    flow.to_csv(results_dir / "cptac_complete_denominator_flow_v2.csv", index=False)
    figures = make_figures(primary, feature, flow, figures_dir)

    table1 = flow.groupby(["cancer", "comparison_layer"], as_index=False).agg(
        registered_program_rows=("program", "size"),
        maximum_matched_unique_cases=("matched_unique_cases", "max"),
        minimum_complete_pair_n=("complete_pair_n", "min"),
        evaluable_program_rows=("status", lambda x: int((x == "evaluable").sum())),
    )
    table2 = heterogeneity.copy()
    table1.to_csv(tables_dir / "TABLE1_CPTAC_CASE_FLOW_V2.csv", index=False)
    table2.to_csv(tables_dir / "TABLE2_CPTAC_PROGRAM_SUMMARY_V2.csv", index=False)
    primary.to_csv(tables_dir / "SUPPLEMENT_TABLE_S1_PRIMARY_REGISTRY_V2.csv", index=False)
    feature.to_csv(tables_dir / "SUPPLEMENT_TABLE_S2_SCORING_SENSITIVITY_V2.csv", index=False)
    residual_context.to_csv(tables_dir / "SUPPLEMENT_TABLE_S3_RESIDUAL_CONTEXT_V2.csv", index=False)

    checks = {
        "protocol_sha256_match": True,
        "registered_primary_rows": int(len(primary)),
        "registered_primary_rows_expected": 160,
        "primary_finite_p_n": int(primary["p_value"].notna().sum()),
        "primary_evaluable_n": int(primary["status"].eq("evaluable").sum()),
        "primary_not_evaluable_n": int(primary["status"].ne("evaluable").sum()),
        "matched_duplicate_key_rows": int(primary_matched.duplicated(["cancer", "program", "comparison_layer", "normalized_case_id"], keep=False).sum()),
        "sample_mapping_crosscheck_failures": int(mapping_audit["status"].ne("pass").sum()),
        "rho_out_of_range_n": int(((primary["spearman_rho"].dropna() < -1) | (primary["spearman_rho"].dropna() > 1)).sum()),
        "p_out_of_range_n": int(((primary["p_value"].dropna() < 0) | (primary["p_value"].dropna() > 1)).sum()),
        "q_out_of_range_n": int(((primary["q_value_bh"].dropna() < 0) | (primary["q_value_bh"].dropna() > 1)).sum()),
        "ci_order_failures": int((primary.dropna(subset=["ci_lower", "spearman_rho", "ci_upper"])["ci_lower"] > primary.dropna(subset=["ci_lower", "spearman_rho", "ci_upper"])["spearman_rho"]).sum() + (primary.dropna(subset=["ci_lower", "spearman_rho", "ci_upper"])["spearman_rho"] > primary.dropna(subset=["ci_lower", "spearman_rho", "ci_upper"])["ci_upper"]).sum()),
        "primary_bca_ci_rows": int(primary["ci_method"].eq("paired_case_bca_bootstrap").sum()),
        "primary_permutation_verified_rows": int(primary["p_value_method"].eq("two_sided_monte_carlo_permutation").sum()),
        "primary_evaluable_tied_rows": int(((primary["status"].eq("evaluable")) & ((primary["x_tie_fraction"] > 0) | (primary["y_tie_fraction"] > 0))).sum()),
        "residual_context_registry_rows": int(len(residual_context)),
        "residual_context_registry_rows_expected": 640,
    }
    checks["critical_failure_n"] = int(
        len(primary) != 160
        or checks["matched_duplicate_key_rows"] != 0
        or checks["sample_mapping_crosscheck_failures"] != 0
        or checks["rho_out_of_range_n"] != 0
        or checks["p_out_of_range_n"] != 0
        or checks["q_out_of_range_n"] != 0
        or checks["ci_order_failures"] != 0
        or len(residual_context) != 640
    )
    write_json(results_dir / "cptac_internal_validation_v2.json", json_safe(checks))
    summary = {
        "generated_at": now_utc(),
        "elapsed_seconds": round(time.time() - started, 3),
        "python": platform.python_version(),
        "protocol_sha256": PROTOCOL_SHA256,
        "program_count": len(PROGRAMS),
        "cancer_count": len(CANCERS),
        "comparison_layer_count": len(LAYERS),
        "unique_cases_any_matched_layer": int(primary_matched["normalized_case_id"].nunique()),
        "primary_registry_rows": len(primary),
        "primary_evaluable_rows": int(primary["status"].eq("evaluable").sum()),
        "primary_fdr_lt_0_05_rows": int((primary["q_value_bh"] < 0.05).sum()),
        "feature_sensitivity_evaluable_rows": int(feature["status"].eq("evaluable").sum()),
        "residual_linear_fdr_lt_0_05_rows": int((residual_context["linear_q_value_bh"] < 0.05).sum()),
        "residual_rank_fdr_lt_0_05_rows": int((residual_context["rank_q_value_bh"] < 0.05).sum()),
        "figures": [str(Path(x).relative_to(analysis_root)) for x in figures],
        "validation": checks,
    }
    write_json(results_dir / "cptac_results_summary_v2.json", json_safe(summary))
    print(json.dumps(json_safe(summary), indent=2, ensure_ascii=False))
    return 0 if checks["critical_failure_n"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
