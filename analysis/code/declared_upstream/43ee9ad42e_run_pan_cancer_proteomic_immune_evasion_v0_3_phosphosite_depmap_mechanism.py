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
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

from run_pan_cancer_proteomic_immune_evasion_discovery_v0_1 import (
    SIGNATURES,
    ensg_stable,
    parse_feature_gene,
    read_header,
    rel,
    resolve_path,
    safe_str,
    zscore_rows,
)
from run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap import (
    DEPMAP_GENE_RE,
    SIGNATURE_ANCHORS,
    SIGNATURE_LABELS,
)


PROJECT_ROOT = Path(
    _workspace_path('${DATA_WORKSPACE}/14_pan_cancer_proteogenomic_immune_evasion')
)
VERSION = "v0_3"
MIN_CORRELATION_N = 20
MOTIF_LABELS = {
    "acidic_ck2_like": "acidic/CK2-like",
    "basophilic_like": "basophilic-like",
    "phosphotyrosine_like": "pY-like",
    "proline_directed_like": "proline-directed",
    "ser_thr_other": "S/T other",
    "ser_thr_sequence_unavailable": "S/T sequence unavailable",
    "ser_thr_site_residue_not_in_peptide": "S/T site not in peptide",
    "unknown_site": "unknown site",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="v0.3 phosphosite, DepMap common-essential, and RNA comparator mechanism audit."
    )
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--version", default=VERSION)
    parser.add_argument("--harmonization-dir", default="04_processed/harmonization_v0_1")
    parser.add_argument("--qc-dir", default="04_processed/first_pass_qc")
    parser.add_argument("--discovery-v0-1-dir", default="04_processed/discovery_v0_1")
    parser.add_argument("--discovery-v0-2-dir", default="04_processed/discovery_v0_2")
    parser.add_argument("--out-dir", default="04_processed/discovery_v0_3")
    parser.add_argument("--figure-dir", default="06_figures/discovery_v0_3")
    parser.add_argument("--report", default="05_reports/discovery_v0_3_report.md")
    parser.add_argument("--summary", default="05_reports/discovery_v0_3_summary.json")
    parser.add_argument("--random-seed", type=int, default=20260624)
    parser.add_argument("--min-correlation-n", type=int, default=MIN_CORRELATION_N)
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


def infer_ann_count(item: pd.Series) -> int:
    raw = item.get("annotation_column_count", 0)
    try:
        ann_count = int(float(raw))
    except Exception:
        ann_count = 0
    source = safe_str(item.get("source"))
    if source == "bcm" and ann_count < 1:
        ann_count = 1
    if source == "umich" and ann_count < 5:
        ann_count = 5
    return ann_count


def signature_lookup(signature_table: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[str]]]:
    symbol_to_sigs: dict[str, set[str]] = defaultdict(set)
    stable_to_sigs: dict[str, set[str]] = defaultdict(set)
    stable_to_symbols: dict[str, set[str]] = defaultdict(set)
    for _, row in signature_table.iterrows():
        signature = safe_str(row.get("signature"))
        gene = safe_str(row.get("gene_symbol"))
        if signature and gene:
            symbol_to_sigs[gene].add(signature)
        for stable in safe_str(row.get("ensembl_stable_ids_from_cptac", "")).split(";"):
            stable = ensg_stable(stable)
            if stable and signature:
                stable_to_sigs[stable].add(signature)
                if gene:
                    stable_to_symbols[stable].add(gene)
    return symbol_to_sigs, stable_to_sigs, stable_to_symbols


def extract_peptide(feature_id: str) -> str:
    raw = safe_str(feature_id)
    pieces = [x.strip() for x in raw.split("|") if x.strip()]
    for token in pieces:
        clean = token.upper()
        clean_no_gap = clean.replace("_", "")
        if 7 <= len(clean_no_gap) <= 60 and re.fullmatch(r"[A-Z_]+", clean):
            if not clean.startswith(("ENSG", "ENSP", "ENST", "OTTHUM")):
                return clean
    match = re.search(r"(?:pepseq|peptide)[:=]([A-Za-z_]{7,60})", raw, flags=re.I)
    if match:
        return match.group(1).upper()
    return ""


def motif_class(phosphosite: str, peptide: str) -> str:
    site = safe_str(phosphosite).upper()
    residue = site[:1]
    if residue == "Y":
        return "phosphotyrosine_like"
    if residue not in {"S", "T"}:
        return "unknown_site"
    seq = safe_str(peptide).upper().replace("_", "")
    if not seq:
        return "ser_thr_sequence_unavailable"
    positions = [i for i, aa in enumerate(seq) if aa == residue]
    if not positions:
        return "ser_thr_site_residue_not_in_peptide"
    center = (len(seq) - 1) / 2
    site_idx = min(positions, key=lambda x: abs(x - center))
    upstream = seq[max(0, site_idx - 3) : site_idx]
    downstream = seq[site_idx + 1 : site_idx + 4]
    next_aa = seq[site_idx + 1 : site_idx + 2]
    if next_aa == "P":
        return "proline_directed_like"
    if any(aa in {"R", "K"} for aa in upstream):
        return "basophilic_like"
    if any(aa in {"D", "E"} for aa in downstream):
        return "acidic_ck2_like"
    return "ser_thr_other"


def spearman_pair(x: np.ndarray, y: np.ndarray, min_n: int) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    n = int(finite.sum())
    if n < min_n:
        return np.nan, np.nan, n
    if len(np.unique(x[finite])) < 3 or len(np.unique(y[finite])) < 3:
        return np.nan, np.nan, n
    rho, p_value = spearmanr(x[finite], y[finite])
    return float(rho), float(p_value), n


def load_phosphosite_inventory(
    harmonization_dir: Path,
    discovery_v0_1_dir: Path,
    out_dir: Path,
) -> pd.DataFrame:
    signature_table = pd.read_csv(discovery_v0_1_dir / "signature_gene_sets_v0_1.csv")
    symbol_to_sigs, stable_to_sigs, stable_to_symbols = signature_lookup(signature_table)
    feature_path = harmonization_dir / "cptac_feature_universe_v0_1.csv"
    rows: list[dict[str, Any]] = []
    usecols = [
        "dataset_id",
        "source",
        "cancer",
        "datatype",
        "tumor_normal_hint",
        "feature_id_raw",
        "gene_symbol_guess",
        "gene_identifier_guess",
        "phosphosite_hint",
        "duplicate_feature_flag",
    ]
    for chunk in pd.read_csv(feature_path, usecols=usecols, chunksize=300_000, low_memory=False):
        chunk = chunk[
            (chunk["datatype"].astype(str) == "phosphoproteomics")
            & (chunk["tumor_normal_hint"].astype(str) != "normal")
        ].copy()
        if chunk.empty:
            continue
        chunk["gene_symbol_guess"] = chunk["gene_symbol_guess"].fillna("").astype(str)
        chunk["gene_identifier_guess"] = chunk["gene_identifier_guess"].fillna("").astype(str)
        chunk["ensembl_stable"] = chunk["gene_identifier_guess"].map(ensg_stable)
        hit = chunk[
            chunk["gene_symbol_guess"].isin(symbol_to_sigs)
            | chunk["ensembl_stable"].isin(stable_to_sigs)
        ].copy()
        if hit.empty:
            continue
        for _, row in hit.iterrows():
            symbol = safe_str(row.get("gene_symbol_guess"))
            stable = ensg_stable(row.get("ensembl_stable", ""))
            signatures = set()
            matched_genes = set()
            if symbol:
                signatures.update(symbol_to_sigs.get(symbol, set()))
            if stable:
                signatures.update(stable_to_sigs.get(stable, set()))
            feature_id = safe_str(row.get("feature_id_raw"))
            peptide = extract_peptide(feature_id)
            phosphosite = safe_str(row.get("phosphosite_hint"))
            mclass = motif_class(phosphosite, peptide)
            for signature in sorted(signatures):
                genes = set()
                if symbol and symbol in SIGNATURES[signature]["genes"]:
                    genes.add(symbol)
                if stable:
                    genes.update(stable_to_symbols.get(stable, set()).intersection(SIGNATURES[signature]["genes"]))
                matched_genes = genes or matched_genes
                rows.append(
                    {
                        "dataset_id": safe_str(row.get("dataset_id")),
                        "source": safe_str(row.get("source")),
                        "cancer": safe_str(row.get("cancer")),
                        "datatype": "phosphoproteomics",
                        "tumor_normal_hint": safe_str(row.get("tumor_normal_hint")),
                        "signature": signature,
                        "signature_label": SIGNATURE_LABELS.get(signature, signature),
                        "feature_id_raw": feature_id,
                        "gene_symbol_guess": symbol,
                        "ensembl_stable": stable,
                        "matched_signature_genes": ";".join(sorted(genes)),
                        "phosphosite_hint": phosphosite,
                        "peptide_sequence_extracted": peptide,
                        "motif_class_exploratory": mclass,
                        "duplicate_feature_flag": int(row.get("duplicate_feature_flag", 0) or 0),
                        "interpretation_boundary": "Exploratory phosphosite inventory; motif classes are sequence-context approximations, not kinase activity validation.",
                    }
                )
    inventory = pd.DataFrame(rows).drop_duplicates(["dataset_id", "feature_id_raw", "signature"])
    inventory.to_csv(out_dir / "signature_phosphosite_feature_inventory_v0_3.csv", index=False)
    return inventory


def protein_vs_phosphoprotein_summary(
    discovery_v0_1_dir: Path,
    discovery_v0_2_dir: Path,
    out_dir: Path,
) -> pd.DataFrame:
    corr = pd.read_csv(discovery_v0_1_dir / "cptac_signature_context_correlations_v0_1.csv")
    scores = pd.read_csv(discovery_v0_1_dir / "cptac_signature_scores_long_v0_1.csv")
    residualized = pd.read_csv(discovery_v0_1_dir / "cptac_signature_residualized_scores_v0_1.csv")
    random_cal = pd.read_csv(discovery_v0_2_dir / "cptac_random_control_calibration_v0_2.csv")
    rows: list[dict[str, Any]] = []
    for signature in SIGNATURES:
        for datatype in ["proteomics", "phosphoproteomics"]:
            sub_corr = corr[
                (corr["signature"].astype(str) == signature)
                & (corr["datatype"].astype(str) == datatype)
                & (corr["tumor_normal_hint"].astype(str) != "normal")
            ].copy()
            immune = sub_corr[sub_corr["context_variable"].astype(str) == "estimate_immune_score"]
            purity = sub_corr[sub_corr["context_variable"].astype(str) == "estimate_tumor_purity"]
            rand = random_cal[
                (random_cal["signature"].astype(str) == signature)
                & (random_cal["datatype"].astype(str) == datatype)
                & (random_cal["status"].astype(str) == "ok")
            ].copy()
            sc = scores[
                (scores["signature"].astype(str) == signature)
                & (scores["datatype"].astype(str) == datatype)
                & (scores["tumor_normal_hint"].astype(str) != "normal")
            ]
            resid = residualized[
                (residualized["signature"].astype(str) == signature)
                & (residualized["datatype"].astype(str) == datatype)
                & (residualized["tumor_normal_hint"].astype(str) != "normal")
            ]
            rows.append(
                {
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "datatype": datatype,
                    "dataset_strata_with_scores": int(sc["dataset_id"].nunique()),
                    "sample_score_rows": int(len(sc)),
                    "immune_score_corr_strata": int(len(immune)),
                    "median_immune_score_rho": float(immune["spearman_rho"].median()) if not immune.empty else np.nan,
                    "median_abs_immune_score_rho": float(immune["spearman_rho"].abs().median()) if not immune.empty else np.nan,
                    "tumor_purity_corr_strata": int(len(purity)),
                    "median_tumor_purity_rho": float(purity["spearman_rho"].median()) if not purity.empty else np.nan,
                    "median_abs_tumor_purity_rho": float(purity["spearman_rho"].abs().median()) if not purity.empty else np.nan,
                    "residualized_dataset_strata": int(resid["dataset_id"].nunique()),
                    "residualized_score_rows": int(len(resid)),
                    "random_control_ok_tests": int(len(rand)),
                    "random_control_empirical_p_lt_0_05_tests": int(
                        (pd.to_numeric(rand["empirical_two_sided_p_abs_rho"], errors="coerce") < 0.05).sum()
                    ),
                    "random_control_median_percentile_abs_rho": float(
                        pd.to_numeric(rand["observed_abs_rho_percentile_vs_controls"], errors="coerce").median()
                    )
                    if not rand.empty
                    else np.nan,
                    "interpretation_boundary": "Platform stratum summary only; protein and phosphoprotein matrices are not paired molecule-level validation experiments.",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "protein_vs_phosphoprotein_stratified_summary_v0_3.csv", index=False)
    return out


def score_phosphosite_coherence(
    root: Path,
    qc_dir: Path,
    discovery_v0_1_dir: Path,
    inventory: pd.DataFrame,
    out_dir: Path,
    min_n: int,
    min_valid_feature_fraction: float,
) -> pd.DataFrame:
    cptac_inventory = pd.read_csv(qc_dir / "cptac_processed_inventory_v0_1.csv")
    cptac_inventory = cptac_inventory.set_index("dataset_id", drop=False)
    scores = pd.read_csv(discovery_v0_1_dir / "cptac_signature_scores_long_v0_1.csv")
    residualized = pd.read_csv(discovery_v0_1_dir / "cptac_signature_residualized_scores_v0_1.csv")
    score_groups = {
        (dataset_id, signature): group[["sample_id", "signature_score"]].dropna()
        for (dataset_id, signature), group in scores.groupby(["dataset_id", "signature"], dropna=False)
    }
    resid_groups = {
        (dataset_id, signature): group[["sample_id", "residualized_signature_score"]].dropna()
        for (dataset_id, signature), group in residualized.groupby(["dataset_id", "signature"], dropna=False)
    }
    rows: list[dict[str, Any]] = []
    if inventory.empty:
        out = pd.DataFrame()
        out.to_csv(out_dir / "signature_phosphosite_coherence_v0_3.csv", index=False)
        return out
    for dataset_id, inv_group in inventory.groupby("dataset_id", dropna=False):
        if dataset_id not in cptac_inventory.index:
            continue
        item = cptac_inventory.loc[dataset_id]
        if isinstance(item, pd.DataFrame):
            item = item.iloc[0]
        if str(item.get("readable")).lower() not in {"true", "1"}:
            continue
        path = resolve_path(safe_str(item["local_path"]), root)
        if not path.exists():
            continue
        ann_count = infer_ann_count(item)
        header = read_header(path, "\t")
        sample_ids = [safe_str(x) for x in header[ann_count:]]
        if not sample_ids:
            continue
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
        if matrix.empty:
            continue
        matrix["_feature_id_raw"] = matrix.iloc[:, 0].astype(str)
        wanted = set(inv_group["feature_id_raw"].astype(str))
        selected = matrix[matrix["_feature_id_raw"].isin(wanted)].copy()
        if selected.empty:
            continue
        numeric = selected.iloc[:, ann_count : ann_count + len(sample_ids)].apply(pd.to_numeric, errors="coerce")
        numeric_values = numeric.to_numpy(dtype=float)
        valid_fraction = np.isfinite(numeric_values).mean(axis=1)
        z_values = zscore_rows(numeric_values)
        inv_by_feature = {
            feature_id: group.copy()
            for feature_id, group in inv_group.groupby("feature_id_raw", dropna=False)
        }
        for local_idx, (_, matrix_row) in enumerate(selected.iterrows()):
            raw = safe_str(matrix_row["_feature_id_raw"])
            if raw not in inv_by_feature:
                continue
            if valid_fraction[local_idx] < min_valid_feature_fraction:
                base_status = "low_valid_feature_fraction"
            elif not np.isfinite(z_values[local_idx]).any():
                base_status = "noninformative_feature"
            else:
                base_status = "ok"
            feature_vector = z_values[local_idx]
            for _, inv_row in inv_by_feature[raw].iterrows():
                signature = safe_str(inv_row["signature"])
                score_df = score_groups.get((dataset_id, signature), pd.DataFrame())
                resid_df = resid_groups.get((dataset_id, signature), pd.DataFrame())
                score_map = dict(zip(score_df.get("sample_id", []), score_df.get("signature_score", [])))
                resid_map = dict(zip(resid_df.get("sample_id", []), resid_df.get("residualized_signature_score", [])))
                score_vector = np.array([score_map.get(sample, np.nan) for sample in sample_ids], dtype=float)
                resid_vector = np.array([resid_map.get(sample, np.nan) for sample in sample_ids], dtype=float)
                if base_status == "ok":
                    rho, p_value, n_pair = spearman_pair(feature_vector, score_vector, min_n)
                    rho_resid, p_resid, n_resid = spearman_pair(feature_vector, resid_vector, min_n)
                    status = "ok" if np.isfinite(rho) else "insufficient_aligned_samples"
                else:
                    rho, p_value, n_pair = np.nan, np.nan, 0
                    rho_resid, p_resid, n_resid = np.nan, np.nan, 0
                    status = base_status
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "source": safe_str(inv_row["source"]),
                        "cancer": safe_str(inv_row["cancer"]),
                        "tumor_normal_hint": safe_str(inv_row.get("tumor_normal_hint")),
                        "signature": signature,
                        "signature_label": SIGNATURE_LABELS.get(signature, signature),
                        "feature_id_raw": raw,
                        "matched_signature_genes": safe_str(inv_row["matched_signature_genes"]),
                        "gene_symbol_guess": safe_str(inv_row["gene_symbol_guess"]),
                        "ensembl_stable": safe_str(inv_row["ensembl_stable"]),
                        "phosphosite_hint": safe_str(inv_row["phosphosite_hint"]),
                        "motif_class_exploratory": safe_str(inv_row["motif_class_exploratory"]),
                        "status": status,
                        "feature_valid_fraction": float(valid_fraction[local_idx]),
                        "feature_nonmissing_n": int(np.isfinite(feature_vector).sum()),
                        "n_signature_score_aligned_samples": n_pair,
                        "rho_to_signature_score": rho,
                        "p_to_signature_score": p_value,
                        "n_residualized_score_aligned_samples": n_resid,
                        "rho_to_residualized_signature_score": rho_resid,
                        "p_to_residualized_signature_score": p_resid,
                        "fdr_family_signature_score": "all_signature_phosphosite_program_coherence_v0_3",
                        "fdr_family_residualized_score": "all_signature_phosphosite_residualized_program_coherence_v0_3",
                        "interpretation_boundary": "Within-program phosphosite coherence. This is partly circular for genes used in the signature and is not independent kinase validation.",
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh_to_signature_score"] = np.nan
        out["fdr_bh_to_residualized_signature_score"] = np.nan
        ok = out["status"].astype(str).eq("ok")
        out.loc[ok, "fdr_bh_to_signature_score"] = bh_fdr(
            pd.to_numeric(out.loc[ok, "p_to_signature_score"], errors="coerce").to_numpy()
        )
        finite_resid = ok & pd.to_numeric(out["p_to_residualized_signature_score"], errors="coerce").notna()
        out.loc[finite_resid, "fdr_bh_to_residualized_signature_score"] = bh_fdr(
            pd.to_numeric(out.loc[finite_resid, "p_to_residualized_signature_score"], errors="coerce").to_numpy()
        )
        out["coherent_abs_rho_ge_0_35"] = (
            pd.to_numeric(out["rho_to_signature_score"], errors="coerce").abs() >= 0.35
        ).astype(int)
        out["residualized_coherent_abs_rho_ge_0_35"] = (
            pd.to_numeric(out["rho_to_residualized_signature_score"], errors="coerce").abs() >= 0.35
        ).astype(int)
    out.to_csv(out_dir / "signature_phosphosite_coherence_v0_3.csv", index=False)
    return out


def summarize_motifs(coherence: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if coherence.empty:
        out = pd.DataFrame()
        out.to_csv(out_dir / "phosphosite_motif_class_summary_v0_3.csv", index=False)
        return out
    ok = coherence[coherence["status"].astype(str) == "ok"].copy()
    ok["abs_rho"] = pd.to_numeric(ok["rho_to_signature_score"], errors="coerce").abs()
    ok["is_coherent"] = (ok["abs_rho"] >= 0.35).astype(int)
    for signature, sig_df in ok.groupby("signature", dropna=False):
        total_coherent = int(sig_df["is_coherent"].sum())
        total_not = int(len(sig_df) - total_coherent)
        for motif, motif_df in sig_df.groupby("motif_class_exploratory", dropna=False):
            a = int(motif_df["is_coherent"].sum())
            b = int(len(motif_df) - a)
            c = total_coherent - a
            d = total_not - b
            try:
                odds, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            except Exception:
                odds, p_value = np.nan, np.nan
            rows.append(
                {
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "motif_class_exploratory": safe_str(motif),
                    "feature_signature_rows": int(len(motif_df)),
                    "coherent_abs_rho_ge_0_35_rows": a,
                    "median_abs_rho_to_signature_score": float(motif_df["abs_rho"].median()) if not motif_df.empty else np.nan,
                    "median_rho_to_signature_score": float(pd.to_numeric(motif_df["rho_to_signature_score"], errors="coerce").median()),
                    "fisher_exact_p_for_coherence_enrichment": p_value,
                    "fisher_exact_odds_ratio": odds,
                    "p_family": "motif_class_coherence_enrichment_within_all_v0_3",
                    "interpretation_boundary": "Motif classes are exploratory sequence-context summaries; no kinase activity or substrate validation is claimed.",
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh_motif_enrichment"] = bh_fdr(
            pd.to_numeric(out["fisher_exact_p_for_coherence_enrichment"], errors="coerce").to_numpy()
        )
    out.to_csv(out_dir / "phosphosite_motif_class_summary_v0_3.csv", index=False)
    return out


def run_depmap_gene_level_audit(
    harmonization_dir: Path,
    qc_dir: Path,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    inv = pd.read_csv(qc_dir / "depmap_inventory_v0_1.csv")
    model_map = pd.read_csv(harmonization_dir / "depmap_model_map_v0_1.csv")
    effect_path = Path(inv[inv["file_role"].astype(str) == "crispr_gene_effect"]["local_path"].iloc[0])
    dependency_path = Path(inv[inv["file_role"].astype(str) == "crispr_dependency"]["local_path"].iloc[0])
    effect_header = pd.read_csv(effect_path, nrows=0).columns.tolist()
    dep_header = pd.read_csv(dependency_path, nrows=0).columns.tolist()
    first_col = effect_header[0]
    dep_first_col = dep_header[0]
    col_to_symbol: dict[str, str] = {}
    for col in effect_header[1:]:
        match = DEPMAP_GENE_RE.match(col)
        if match:
            col_to_symbol[col] = match.group("symbol")
    signature_symbols = {gene for spec in SIGNATURES.values() for gene in spec["genes"]}
    selected_cols = [col for col, symbol in col_to_symbol.items() if symbol in signature_symbols]
    selected_dep_cols = [col for col in selected_cols if col in dep_header]
    effect = pd.read_csv(effect_path)
    effect = effect.rename(columns={first_col: "ModelID"})
    dependency = pd.read_csv(dependency_path, usecols=[dep_first_col] + selected_dep_cols)
    dependency = dependency.rename(columns={dep_first_col: "ModelID"})
    model_cols = ["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "AgeCategory"]
    model_lineage = model_map[[c for c in model_cols if c in model_map.columns]].copy()
    effect = effect.merge(model_lineage, on="ModelID", how="left", validate="many_to_one")
    dependency = dependency.merge(model_lineage, on="ModelID", how="left", validate="many_to_one")
    all_gene_cols = [col for col in effect_header[1:] if col in effect.columns]
    bg_matrix = effect[all_gene_cols].apply(pd.to_numeric, errors="coerce")
    bg_gene_medians = bg_matrix.median(axis=0, skipna=True)
    cns_mask = (
        effect.get("OncotreeLineage", pd.Series(index=effect.index, dtype=str)).astype(str).str.contains("CNS|Brain", case=False, na=False)
        | effect.get("OncotreePrimaryDisease", pd.Series(index=effect.index, dtype=str)).astype(str).str.contains("CNS|Brain", case=False, na=False)
    )
    dep_cns_mask = (
        dependency.get("OncotreeLineage", pd.Series(index=dependency.index, dtype=str)).astype(str).str.contains("CNS|Brain", case=False, na=False)
        | dependency.get("OncotreePrimaryDisease", pd.Series(index=dependency.index, dtype=str)).astype(str).str.contains("CNS|Brain", case=False, na=False)
    )
    rows: list[dict[str, Any]] = []
    for signature, spec in SIGNATURES.items():
        genes = set(spec["genes"])
        cols = [col for col, symbol in col_to_symbol.items() if symbol in genes and col in effect.columns]
        for col in cols:
            symbol = col_to_symbol.get(col, col)
            effect_values = pd.to_numeric(effect[col], errors="coerce")
            dep_values = pd.to_numeric(dependency[col], errors="coerce") if col in dependency.columns else pd.Series(dtype=float)
            median_effect = float(effect_values.median(skipna=True))
            background_percentile = float(np.mean(bg_gene_medians <= median_effect) * 100) if len(bg_gene_medians) else np.nan
            rows.append(
                {
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "gene_symbol": symbol,
                    "depmap_feature_id": col,
                    "all_model_n": int(effect_values.notna().sum()),
                    "cns_brain_model_n": int(effect_values.loc[cns_mask].notna().sum()),
                    "median_gene_effect_all_models": median_effect,
                    "median_gene_effect_cns_brain_models": float(effect_values.loc[cns_mask].median(skipna=True)) if cns_mask.any() else np.nan,
                    "mean_dependency_probability_all_models": float(dep_values.mean(skipna=True)) if not dep_values.empty else np.nan,
                    "mean_dependency_probability_cns_brain_models": float(dep_values.loc[dep_cns_mask].mean(skipna=True)) if (not dep_values.empty and dep_cns_mask.any()) else np.nan,
                    "background_percentile_lower_is_more_dependent": background_percentile,
                    "common_dependency_probability_ge_0_5": int(float(dep_values.mean(skipna=True)) >= 0.5) if not dep_values.empty else 0,
                    "global_dependency_top5_percent": int(background_percentile <= 5) if np.isfinite(background_percentile) else 0,
                    "interpretation_boundary": "Gene-level DepMap supportive context; not target validation and not therapy sensitivity evidence.",
                }
            )
    gene_level = pd.DataFrame(rows)
    gene_level.to_csv(out_dir / "depmap_signature_gene_level_common_essential_audit_v0_3.csv", index=False)

    summary_rows: list[dict[str, Any]] = []
    for signature, spec in SIGNATURES.items():
        sig = gene_level[gene_level["signature"].astype(str) == signature].copy()
        if sig.empty:
            summary_rows.append(
                {
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "mapped_gene_count": 0,
                    "status": "no_depmap_mapping",
                    "interpretation": "No DepMap context available.",
                }
            )
            continue
        common_or_global = (
            sig["common_dependency_probability_ge_0_5"].astype(int).eq(1)
            | sig["global_dependency_top5_percent"].astype(int).eq(1)
        )
        kept = sig[~common_or_global]
        before = pd.to_numeric(sig["median_gene_effect_all_models"], errors="coerce")
        after = pd.to_numeric(kept["median_gene_effect_all_models"], errors="coerce")
        try:
            mwu_p = mannwhitneyu(before.dropna(), bg_gene_medians.dropna(), alternative="two-sided").pvalue
        except Exception:
            mwu_p = np.nan
        summary_rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "mapped_gene_count": int(len(sig)),
                "common_dependency_probability_ge_0_5_gene_count": int(sig["common_dependency_probability_ge_0_5"].sum()),
                "global_dependency_top5_percent_gene_count": int(sig["global_dependency_top5_percent"].sum()),
                "common_or_global_excluded_gene_count": int(common_or_global.sum()),
                "remaining_after_common_or_global_exclusion": int((~common_or_global).sum()),
                "median_gene_effect_all_models_before_exclusion": float(before.median(skipna=True)),
                "median_gene_effect_all_models_after_common_or_global_exclusion": float(after.median(skipna=True)) if len(after) else np.nan,
                "median_dependency_probability_all_models": float(pd.to_numeric(sig["mean_dependency_probability_all_models"], errors="coerce").median(skipna=True)),
                "more_dependent_than_background_mannwhitney_p": mwu_p,
                "excluded_genes_common_or_global": ";".join(sorted(sig.loc[common_or_global, "gene_symbol"].astype(str).unique())),
                "retained_genes_after_exclusion": ";".join(sorted(kept["gene_symbol"].astype(str).unique())),
                "status": "ok",
                "interpretation": "If dependency signal attenuates after excluding common/global dependency genes, interpret as common-essential context rather than specific vulnerability.",
            }
        )
    exclusion = pd.DataFrame(summary_rows)
    if not exclusion.empty:
        exclusion["more_dependent_than_background_fdr_bh"] = bh_fdr(
            pd.to_numeric(exclusion["more_dependent_than_background_mannwhitney_p"], errors="coerce").to_numpy()
        )
    exclusion.to_csv(out_dir / "depmap_common_essential_exclusion_summary_v0_3.csv", index=False)
    return gene_level, exclusion


def rna_residual_comparator_summary(
    discovery_v0_2_dir: Path,
    protein_summary: pd.DataFrame,
    depmap_exclusion: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    tcga_corr = pd.read_csv(discovery_v0_2_dir / "tcga_rna_signature_context_correlations_v0_2.csv")
    platform_summary = pd.read_csv(discovery_v0_2_dir / "platform_comparator_summary_v0_2.csv")
    rows: list[dict[str, Any]] = []
    for signature in SIGNATURES:
        anchors = set(SIGNATURE_ANCHORS.get(signature, []))
        sub = tcga_corr[
            (tcga_corr["signature"].astype(str) == signature)
            & (tcga_corr["tcga_anchor_variable"].astype(str).isin(anchors))
        ].copy()
        cptac = platform_summary[platform_summary["signature"].astype(str) == signature]
        dep = depmap_exclusion[depmap_exclusion["signature"].astype(str) == signature]
        prot = protein_summary[
            (protein_summary["signature"].astype(str) == signature)
            & (protein_summary["datatype"].astype(str) == "proteomics")
        ]
        phospho = protein_summary[
            (protein_summary["signature"].astype(str) == signature)
            & (protein_summary["datatype"].astype(str) == "phosphoproteomics")
        ]
        median_rna = float(sub["spearman_rho"].median()) if not sub.empty else np.nan
        positive_fraction = float((sub["spearman_rho"] > 0).mean()) if not sub.empty else np.nan
        if np.isfinite(median_rna) and median_rna >= 0.6 and positive_fraction >= 0.8:
            comparator_class = "strong_rna_concordant_context"
        elif np.isfinite(median_rna) and median_rna >= 0.25:
            comparator_class = "moderate_rna_concordant_context"
        else:
            comparator_class = "weak_or_signature_specific_rna_context"
        rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "tcga_anchor_variables": ";".join(sorted(anchors)),
                "tcga_anchor_corr_tests": int(len(sub)),
                "tcga_median_anchor_rho": median_rna,
                "tcga_positive_anchor_fraction": positive_fraction,
                "rna_comparator_class": comparator_class,
                "cptac_v0_2_median_immune_score_rho": float(cptac["cptac_median_immune_score_rho"].iloc[0]) if not cptac.empty else np.nan,
                "proteomics_median_immune_score_rho": float(prot["median_immune_score_rho"].iloc[0]) if not prot.empty else np.nan,
                "phosphoproteomics_median_immune_score_rho": float(phospho["median_immune_score_rho"].iloc[0]) if not phospho.empty else np.nan,
                "random_control_median_percentile_abs_rho": float(cptac["random_control_median_percentile_abs_rho"].iloc[0]) if not cptac.empty else np.nan,
                "depmap_remaining_after_common_or_global_exclusion": int(dep["remaining_after_common_or_global_exclusion"].iloc[0]) if not dep.empty else 0,
                "depmap_median_gene_effect_after_exclusion": float(dep["median_gene_effect_all_models_after_common_or_global_exclusion"].iloc[0]) if not dep.empty else np.nan,
                "interpretation_boundary": "TCGA RNA comparator only; not external validation and not sample-overlap analysis.",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "rna_residual_comparator_summary_v0_3.csv", index=False)
    return out


def build_claim_tier_table(
    protein_summary: pd.DataFrame,
    motif_summary: pd.DataFrame,
    depmap_exclusion: pd.DataFrame,
    rna_summary: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for signature in SIGNATURES:
        rna = rna_summary[rna_summary["signature"].astype(str) == signature]
        dep = depmap_exclusion[depmap_exclusion["signature"].astype(str) == signature]
        motif = motif_summary[motif_summary["signature"].astype(str) == signature]
        phospho = protein_summary[
            (protein_summary["signature"].astype(str) == signature)
            & (protein_summary["datatype"].astype(str) == "phosphoproteomics")
        ]
        rna_class = safe_str(rna["rna_comparator_class"].iloc[0]) if not rna.empty else ""
        common_excluded = int(dep["common_or_global_excluded_gene_count"].iloc[0]) if not dep.empty else 0
        remaining_depmap = int(dep["remaining_after_common_or_global_exclusion"].iloc[0]) if not dep.empty else 0
        motif_rows = int(motif["feature_signature_rows"].sum()) if not motif.empty else 0
        coherent_rows = int(motif["coherent_abs_rho_ge_0_35_rows"].sum()) if not motif.empty else 0
        phospho_percentile = float(phospho["random_control_median_percentile_abs_rho"].iloc[0]) if not phospho.empty else np.nan
        if signature == "proteasome_antigen_processing":
            tier = "common_essential_rich_dependency_context"
            allowed = "Proteasome program shows low immune-context coupling and strong common/global dependency context; use as supportive dependency context only."
        elif rna_class == "strong_rna_concordant_context":
            tier = "rna_concordant_immune_ecology_program"
            allowed = "Program is consistently reflected by TCGA RNA immune anchors and CPTAC protein/phosphoprotein context; describe as RNA-concordant immune ecology."
        else:
            tier = "exploratory_multiomic_context"
            allowed = "Exploratory cross-platform context only; avoid biomarker or therapeutic claims."
        rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "v0_3_claim_tier": tier,
                "rna_comparator_class": rna_class,
                "phosphoproteomics_random_control_median_percentile_abs_rho": phospho_percentile,
                "phosphosite_motif_feature_signature_rows": motif_rows,
                "phosphosite_coherent_abs_rho_ge_0_35_rows": coherent_rows,
                "depmap_common_or_global_excluded_gene_count": common_excluded,
                "depmap_remaining_after_exclusion": remaining_depmap,
                "allowed_claim": allowed,
                "prohibited_claims": "Do not claim validation, biomarker performance, kinase activity, tumor-cell-intrinsic immune evasion, therapeutic vulnerability, or causal mechanism.",
                "recommended_manuscript_wording": "supportive context; exploratory motif class; RNA-concordant immune ecology; common-essential-rich dependency context",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "discovery_claim_tier_table_v0_3.csv", index=False)
    return out


def write_figures(
    fig_dir: Path,
    protein_summary: pd.DataFrame,
    motif_summary: pd.DataFrame,
    depmap_exclusion: pd.DataFrame,
    rna_summary: pd.DataFrame,
) -> list[Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    heat = protein_summary.pivot(index="signature_label", columns="datatype", values="median_immune_score_rho")
    if not heat.empty:
        path = fig_dir / "fig07_protein_vs_phosphoprotein_context_summary_v0_3.png"
        fig, ax = plt.subplots(figsize=(6.6, 4.8))
        data = heat.reindex([SIGNATURE_LABELS.get(sig, sig) for sig in SIGNATURES])
        im = ax.imshow(data.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
        ax.set_xticks(range(len(data.columns)))
        ax.set_xticklabels([c.replace("phosphoproteomics", "phospho").replace("proteomics", "protein") for c in data.columns], rotation=0)
        ax.set_yticks(range(len(data.index)))
        ax.set_yticklabels(data.index)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                value = data.iloc[i, j]
                if np.isfinite(value):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
        ax.set_title("CPTAC immune-score correlation by platform")
        fig.colorbar(im, ax=ax, shrink=0.8, label="Median Spearman rho")
        fig.tight_layout()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)

    if not motif_summary.empty:
        path = fig_dir / "fig08_phosphosite_motif_class_summary_v0_3.png"
        motif_plot = motif_summary.copy()
        motif_plot["motif_label"] = motif_plot["motif_class_exploratory"].map(
            lambda x: MOTIF_LABELS.get(safe_str(x), safe_str(x))
        )
        motif_counts = motif_plot.pivot_table(
            index="signature_label",
            columns="motif_label",
            values="feature_signature_rows",
            aggfunc="sum",
            fill_value=0,
        )
        motif_counts = motif_counts.reindex([SIGNATURE_LABELS.get(sig, sig) for sig in SIGNATURES]).fillna(0)
        motif_frac = motif_counts.div(motif_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        bottom = np.zeros(len(motif_frac))
        for col in motif_frac.columns:
            values = motif_frac[col].to_numpy(dtype=float)
            ax.barh(motif_frac.index, values, left=bottom, label=col)
            bottom += values
        ax.set_xlabel("Fraction of signature-matched phosphosite rows")
        ax.set_title("Exploratory phosphosite motif-class composition")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
        fig.tight_layout()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)

    if not depmap_exclusion.empty:
        path = fig_dir / "fig09_depmap_common_essential_sensitivity_v0_3.png"
        data = depmap_exclusion.copy()
        data["signature_label"] = data["signature"].map(lambda x: SIGNATURE_LABELS.get(x, x))
        data = data.set_index("signature_label").reindex([SIGNATURE_LABELS.get(sig, sig) for sig in SIGNATURES])
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        x = np.arange(len(data))
        ax.bar(x - 0.18, data["common_or_global_excluded_gene_count"].fillna(0), width=0.36, label="Excluded common/global")
        ax.bar(x + 0.18, data["remaining_after_common_or_global_exclusion"].fillna(0), width=0.36, label="Retained")
        ax.set_xticks(x)
        ax.set_xticklabels(data.index, rotation=35, ha="right")
        ax.set_ylabel("DepMap-mapped genes")
        ax.set_title("DepMap common/global dependency sensitivity")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)

    if not rna_summary.empty:
        path = fig_dir / "fig10_rna_cptac_dependency_claim_map_v0_3.png"
        data = rna_summary.copy()
        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        x = pd.to_numeric(data["tcga_median_anchor_rho"], errors="coerce")
        y = pd.to_numeric(data["random_control_median_percentile_abs_rho"], errors="coerce")
        sizes = 40 + 12 * pd.to_numeric(data["depmap_remaining_after_common_or_global_exclusion"], errors="coerce").fillna(0)
        ax.scatter(x, y, s=sizes, alpha=0.75, edgecolor="black", linewidth=0.5)
        label_offsets = {
            "antigen_presentation_mhc_i": (-30, -18),
            "antigen_presentation_mhc_ii": (8, 8),
            "ifn_gamma_response": (-44, -18),
            "cytolytic_t_cell_context": (-34, 10),
            "checkpoint_exhaustion_context": (-58, 18),
            "myeloid_inflammatory_context": (6, -16),
            "tgfb_emt_exclusion_context": (-18, -24),
            "proteasome_antigen_processing": (0, 10),
        }
        for _, row in data.iterrows():
            sx = float(row["tcga_median_anchor_rho"]) if np.isfinite(row["tcga_median_anchor_rho"]) else 0
            sy = float(row["random_control_median_percentile_abs_rho"]) if np.isfinite(row["random_control_median_percentile_abs_rho"]) else 0
            dx, dy = label_offsets.get(safe_str(row["signature"]), (4, 4))
            ax.annotate(
                SIGNATURE_LABELS.get(row["signature"], row["signature"]),
                xy=(sx, sy),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=7,
                ha="left",
                va="bottom",
                bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none", "pad": 1.0},
            )
        ax.set_xlabel("TCGA RNA immune-anchor median rho")
        ax.set_ylabel("CPTAC random-control percentile")
        ax.set_title("Claim map: RNA concordance and CPTAC context")
        ax.set_xlim(0.25, 1.04)
        ax.set_ylim(45, 106)
        fig.tight_layout()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)
    return figures


def write_report(
    report_path: Path,
    summary: dict[str, Any],
    protein_summary: pd.DataFrame,
    inventory: pd.DataFrame,
    coherence: pd.DataFrame,
    motif_summary: pd.DataFrame,
    depmap_exclusion: pd.DataFrame,
    rna_summary: pd.DataFrame,
    claim: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append("# Discovery v0.3 phosphosite and dependency mechanism audit")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(
        "This version adds protein-vs-phosphoprotein stratification, signature-matched phosphosite coherence, exploratory motif classes, DepMap common/global dependency sensitivity, and a TCGA RNA comparator summary. It does not claim validation, biomarker performance, kinase activity, therapeutic vulnerability, or causality."
    )
    lines.append("")
    lines.append("## Key counts")
    lines.append("")
    for key in [
        "protein_vs_phosphoprotein_rows",
        "phosphosite_inventory_rows",
        "phosphosite_coherence_rows",
        "motif_summary_rows",
        "depmap_gene_level_rows",
        "depmap_exclusion_rows",
        "rna_comparator_rows",
        "claim_tier_rows",
    ]:
        lines.append(f"- {key}: {summary.get(key)}")
    lines.append("")
    lines.append("## Main results")
    lines.append("")
    if not rna_summary.empty:
        strong = rna_summary[rna_summary["rna_comparator_class"].astype(str) == "strong_rna_concordant_context"]
        lines.append(
            f"- TCGA RNA comparator: {len(strong)}/{len(rna_summary)} signatures are classified as strong RNA-concordant context signals."
        )
    if not coherence.empty:
        ok = coherence[coherence["status"].astype(str) == "ok"]
        coherent = int(ok.get("coherent_abs_rho_ge_0_35", pd.Series(dtype=int)).sum())
        resid = int(ok.get("residualized_coherent_abs_rho_ge_0_35", pd.Series(dtype=int)).sum())
        lines.append(
            f"- Phosphosite coherence: {len(ok)} evaluable feature-signature rows; {coherent} have abs(rho) >= 0.35 to the parent program score and {resid} retain abs(rho) >= 0.35 after residualized program scoring."
        )
    if not depmap_exclusion.empty:
        proteasome = depmap_exclusion[depmap_exclusion["signature"].astype(str) == "proteasome_antigen_processing"]
        if not proteasome.empty:
            row = proteasome.iloc[0]
            lines.append(
                f"- Proteasome DepMap audit: {int(row['common_or_global_excluded_gene_count'])} genes are common/global dependency flagged; {int(row['remaining_after_common_or_global_exclusion'])} remain after exclusion."
            )
    lines.append("")
    lines.append("## Claim boundaries")
    lines.append("")
    lines.append("- Use `supportive context`, `RNA-concordant immune ecology`, `exploratory phosphosite motif class`, and `common-essential-rich dependency context` wording.")
    lines.append("- Avoid `validation`, `biomarker`, `mechanism`, `kinase activity`, `therapeutic target`, and `tumor-cell-intrinsic immune evasion` unless supported by independent experimental data.")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    for name, path in summary.get("outputs", {}).items():
        lines.append(f"- {name}: `{path}`")
    lines.append("")
    lines.append("## Figures")
    lines.append("")
    for path in summary.get("figures", []):
        lines.append(f"- `{path}`")
    lines.append("")
    lines.append("## Claim tier table")
    lines.append("")
    if not claim.empty:
        cols = ["signature_label", "v0_3_claim_tier", "allowed_claim"]
        lines.append(claim[cols].to_markdown(index=False))
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    start = time.time()
    args = parse_args()
    np.random.seed(args.random_seed)
    root = Path(args.project_root)
    harmonization_dir = resolve_path(args.harmonization_dir, root)
    qc_dir = resolve_path(args.qc_dir, root)
    discovery_v0_1_dir = resolve_path(args.discovery_v0_1_dir, root)
    discovery_v0_2_dir = resolve_path(args.discovery_v0_2_dir, root)
    out_dir = resolve_path(args.out_dir, root)
    fig_dir = resolve_path(args.figure_dir, root)
    report_path = resolve_path(args.report, root)
    summary_path = resolve_path(args.summary, root)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    protein_summary = protein_vs_phosphoprotein_summary(discovery_v0_1_dir, discovery_v0_2_dir, out_dir)
    inventory = load_phosphosite_inventory(harmonization_dir, discovery_v0_1_dir, out_dir)
    coherence = score_phosphosite_coherence(
        root=root,
        qc_dir=qc_dir,
        discovery_v0_1_dir=discovery_v0_1_dir,
        inventory=inventory,
        out_dir=out_dir,
        min_n=args.min_correlation_n,
        min_valid_feature_fraction=args.min_valid_feature_fraction,
    )
    motif_summary = summarize_motifs(coherence, out_dir)
    depmap_gene_level, depmap_exclusion = run_depmap_gene_level_audit(harmonization_dir, qc_dir, out_dir)
    rna_summary = rna_residual_comparator_summary(discovery_v0_2_dir, protein_summary, depmap_exclusion, out_dir)
    claim = build_claim_tier_table(protein_summary, motif_summary, depmap_exclusion, rna_summary, out_dir)
    figures = write_figures(fig_dir, protein_summary, motif_summary, depmap_exclusion, rna_summary)

    summary = {
        "version": args.version,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "project_root": str(root),
        "protein_vs_phosphoprotein_rows": int(len(protein_summary)),
        "phosphosite_inventory_rows": int(len(inventory)),
        "phosphosite_inventory_signature_count": int(inventory["signature"].nunique()) if not inventory.empty else 0,
        "phosphosite_coherence_rows": int(len(coherence)),
        "phosphosite_coherence_ok_rows": int((coherence["status"].astype(str) == "ok").sum()) if not coherence.empty else 0,
        "motif_summary_rows": int(len(motif_summary)),
        "depmap_gene_level_rows": int(len(depmap_gene_level)),
        "depmap_exclusion_rows": int(len(depmap_exclusion)),
        "rna_comparator_rows": int(len(rna_summary)),
        "claim_tier_rows": int(len(claim)),
        "figures": [rel(path, root) for path in figures],
        "outputs": {
            "protein_vs_phosphoprotein": rel(out_dir / "protein_vs_phosphoprotein_stratified_summary_v0_3.csv", root),
            "phosphosite_inventory": rel(out_dir / "signature_phosphosite_feature_inventory_v0_3.csv", root),
            "phosphosite_coherence": rel(out_dir / "signature_phosphosite_coherence_v0_3.csv", root),
            "motif_summary": rel(out_dir / "phosphosite_motif_class_summary_v0_3.csv", root),
            "depmap_gene_level": rel(out_dir / "depmap_signature_gene_level_common_essential_audit_v0_3.csv", root),
            "depmap_common_essential_exclusion": rel(out_dir / "depmap_common_essential_exclusion_summary_v0_3.csv", root),
            "rna_comparator": rel(out_dir / "rna_residual_comparator_summary_v0_3.csv", root),
            "claim_tier": rel(out_dir / "discovery_claim_tier_table_v0_3.csv", root),
        },
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "elapsed_seconds": round(time.time() - start, 2),
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    write_report(
        report_path,
        summary,
        protein_summary,
        inventory,
        coherence,
        motif_summary,
        depmap_exclusion,
        rna_summary,
        claim,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
