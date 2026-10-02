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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(
    _workspace_os.environ.get("ARTICLE002_DATA_ROOT", str(Path.cwd()))
).expanduser().resolve()
VERSION = "v0_1"

GENE_SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9.-]{1,20}$")
PHOSPHOSITE_RE = re.compile(r"^[STY]\d+", re.I)
UUID_D_RE = re.compile(r"^(?P<uuid>[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{6,12})(?:_D\d+)?$", re.I)


SIGNATURES: dict[str, dict[str, Any]] = {
    "antigen_presentation_mhc_i": {
        "description": "MHC-I antigen presentation and peptide loading",
        "provenance": "Curated from canonical antigen-processing genes; v0.1 hypothesis set",
        "genes": [
            "HLA-A",
            "HLA-B",
            "HLA-C",
            "B2M",
            "TAP1",
            "TAP2",
            "TAPBP",
            "NLRC5",
            "PSMB8",
            "PSMB9",
            "PSMB10",
            "ERAP1",
            "ERAP2",
            "CALR",
            "CANX",
            "PDIA3",
        ],
    },
    "antigen_presentation_mhc_ii": {
        "description": "MHC-II antigen presentation and professional antigen-presenting context",
        "provenance": "Curated from canonical MHC-II pathway genes; v0.1 hypothesis set",
        "genes": [
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
        ],
    },
    "ifn_gamma_response": {
        "description": "Interferon response and downstream inflammatory signaling",
        "provenance": "Curated from Hallmark-like IFN response and antigen-response genes; v0.1 hypothesis set",
        "genes": [
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
        ],
    },
    "cytolytic_t_cell_context": {
        "description": "Cytotoxic T/NK-cell context",
        "provenance": "Curated cytolytic lymphocyte marker genes; v0.1 hypothesis set",
        "genes": [
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
        ],
    },
    "checkpoint_exhaustion_context": {
        "description": "Checkpoint and T-cell exhaustion context",
        "provenance": "Curated checkpoint/exhaustion marker genes; v0.1 hypothesis set",
        "genes": [
            "PDCD1",
            "CD274",
            "PDCD1LG2",
            "CTLA4",
            "LAG3",
            "HAVCR2",
            "TIGIT",
            "TOX",
            "ENTPD1",
            "CD80",
            "CD86",
            "ICOS",
            "TNFRSF9",
        ],
    },
    "myeloid_inflammatory_context": {
        "description": "Myeloid/macrophage inflammatory context",
        "provenance": "Curated monocyte/macrophage inflammatory marker genes; v0.1 hypothesis set",
        "genes": [
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
        ],
    },
    "tgfb_emt_exclusion_context": {
        "description": "TGF-beta, stromal, EMT, and exclusion-like context",
        "provenance": "Curated stromal/TGF-beta/EMT marker genes; v0.1 hypothesis set",
        "genes": [
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
        ],
    },
    "proteasome_antigen_processing": {
        "description": "Proteasome and antigen-processing machinery",
        "provenance": "Curated proteasome and antigen-processing machinery genes; v0.1 hypothesis set",
        "genes": [
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
        ],
    },
}

CONTEXT_VARS = [
    "estimate_tumor_purity",
    "estimate_immune_score",
    "estimate_stromal_score",
    "cibersort_t_cell",
    "cibersort_myeloid",
    "cibersort_cytotoxic",
    "xcell_t_cell",
    "xcell_myeloid",
    "xcell_stroma",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="First-pass CPTAC proteomic/phosphoproteomic immune-evasion discovery."
    )
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--version", default=VERSION)
    parser.add_argument("--harmonization-dir", default="04_processed/harmonization_v0_1")
    parser.add_argument("--qc-dir", default="04_processed/first_pass_qc")
    parser.add_argument("--out-dir", default="04_processed/discovery_v0_1")
    parser.add_argument("--figure-dir", default="06_figures/discovery_v0_1")
    parser.add_argument("--report", default="05_reports/discovery_v0_1_report.md")
    parser.add_argument("--summary", default="05_reports/discovery_v0_1_summary.json")
    parser.add_argument("--min-correlation-n", type=int, default=20)
    parser.add_argument("--min-residual-n", type=int, default=25)
    return parser.parse_args()


def safe_str(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value)


def rel(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except Exception:
        return str(path)


def resolve_path(value: str, root: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def read_header(path: Path, delimiter: str = "\t") -> list[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            if row and any(safe_str(x).strip() for x in row):
                return [safe_str(x).strip() for x in row]
    return []


def infer_cptac_case_id(sample_id: str) -> str:
    s = safe_str(sample_id).strip()
    if not s:
        return ""
    if s.upper().startswith(("QC", "JHU-QC")) or "QC" in s.upper():
        return ""
    s = re.sub(r"-(T|N)$", "", s, flags=re.I)
    uuid_match = UUID_D_RE.match(s)
    if uuid_match:
        return uuid_match.group("uuid").lower()
    s = re.sub(r"_D\d+$", "", s, flags=re.I)
    return s


def ensg_stable(value: str) -> str:
    text = safe_str(value).strip()
    if not text:
        return ""
    match = re.search(r"(ENSG\d+)(?:\.\d+)?", text)
    return match.group(1) if match else ""


def parse_feature_gene(raw: str, row_values: list[Any], header: list[str], source: str) -> tuple[str, str, str]:
    raw = safe_str(raw).strip()
    lower = {h.lower(): i for i, h in enumerate(header)}
    gene_identifier = ""
    gene_symbol = ""
    phosphosite = ""
    gene_col = lower.get("gene")
    if gene_col is not None and gene_col < len(row_values):
        gene_identifier = safe_str(row_values[gene_col]).strip()
    pieces = raw.split("|")
    stable = ensg_stable(raw) or ensg_stable(gene_identifier)
    if stable and not gene_identifier:
        gene_identifier = stable
    for token in re.split(r"[|_]", raw):
        if PHOSPHOSITE_RE.match(token):
            phosphosite = token
            break
    if source == "umich":
        for i, candidate in enumerate(pieces):
            prev_token = pieces[i - 1] if i else ""
            if (
                prev_token.endswith("-201")
                and GENE_SYMBOL_RE.match(candidate)
                and not candidate.startswith(("ENSG", "ENSP", "ENST", "OTTHUM"))
            ):
                gene_symbol = candidate
                break
    return gene_symbol, stable, phosphosite


def build_signature_table(
    harmonization_dir: Path, out_path: Path
) -> tuple[pd.DataFrame, dict[str, set[str]], dict[str, set[tuple[str, str]]]]:
    signature_rows = []
    all_symbols = {gene for spec in SIGNATURES.values() for gene in spec["genes"]}
    symbol_to_ensembl: dict[str, set[str]] = defaultdict(set)
    feature_path = harmonization_dir / "cptac_feature_universe_v0_1.csv"
    if feature_path.exists():
        usecols = ["source", "gene_symbol_guess", "gene_identifier_guess"]
        for chunk in pd.read_csv(feature_path, usecols=usecols, chunksize=500_000, low_memory=False):
            chunk["gene_symbol_guess"] = chunk["gene_symbol_guess"].fillna("").astype(str)
            hit = chunk[chunk["gene_symbol_guess"].isin(all_symbols)].copy()
            for _, row in hit.iterrows():
                stable = ensg_stable(row["gene_identifier_guess"])
                if stable:
                    symbol_to_ensembl[row["gene_symbol_guess"]].add(stable)
    signature_symbol_map: dict[str, set[str]] = defaultdict(set)
    signature_ensembl_map: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for signature, spec in SIGNATURES.items():
        for gene in spec["genes"]:
            signature_symbol_map[gene].add(signature)
            ensembl_ids = sorted(symbol_to_ensembl.get(gene, set()))
            for stable in ensembl_ids:
                signature_ensembl_map[stable].add((gene, signature))
            signature_rows.append(
                {
                    "signature": signature,
                    "gene_symbol": gene,
                    "description": spec["description"],
                    "provenance": spec["provenance"],
                    "ensembl_stable_ids_from_cptac": ";".join(ensembl_ids),
                    "has_cptac_ensembl_mapping": int(bool(ensembl_ids)),
                }
            )
    table = pd.DataFrame(signature_rows)
    table.to_csv(out_path, index=False)
    return table, signature_symbol_map, signature_ensembl_map


def zscore_rows(values: np.ndarray) -> np.ndarray:
    means = np.nanmean(values, axis=1, keepdims=True)
    stds = np.nanstd(values, axis=1, ddof=1, keepdims=True)
    stds[(stds == 0) | ~np.isfinite(stds)] = np.nan
    return (values - means) / stds


def score_one_matrix(
    root: Path,
    item: pd.Series,
    sample_map: pd.DataFrame,
    signature_table: pd.DataFrame,
    signature_symbol_map: dict[str, set[str]],
    signature_ensembl_map: dict[str, set[tuple[str, str]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    dataset_id = safe_str(item["dataset_id"])
    source = safe_str(item["source"])
    cancer = safe_str(item["cancer"])
    datatype = safe_str(item["datatype"])
    tumor_normal_hint = safe_str(item["tumor_normal_hint"])
    path = resolve_path(safe_str(item["local_path"]), root)
    header = read_header(path, "\t")
    if not header:
        return [], []
    ann_count = int(item.get("annotation_column_count", 0) or 0)
    if source == "bcm" and ann_count < 1:
        ann_count = 1
    if source == "umich" and ann_count < 5:
        ann_count = 5
    sample_ids = header[ann_count:]
    names = list(range(len(header)))
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        names=names,
        na_values=["NA", "NaN", "nan", ""],
        keep_default_na=True,
        low_memory=False,
    )
    if df.empty or not sample_ids:
        return [], []
    feature_raw = df.iloc[:, 0].astype(str).tolist()
    row_to_signatures: list[set[str]] = []
    row_to_signature_genes: list[dict[str, set[str]]] = []
    row_ensembl: list[str] = []
    for idx, raw in enumerate(feature_raw):
        row_values = df.iloc[idx, :ann_count].tolist()
        symbol, stable, _phosphosite = parse_feature_gene(raw, row_values, header, source)
        signatures = set()
        signature_genes: dict[str, set[str]] = defaultdict(set)
        if symbol:
            for signature in signature_symbol_map.get(symbol, set()):
                signatures.add(signature)
                signature_genes[signature].add(symbol)
        if stable:
            for gene, signature in signature_ensembl_map.get(stable, set()):
                signatures.add(signature)
                signature_genes[signature].add(gene)
        row_to_signatures.append(signatures)
        row_to_signature_genes.append(signature_genes)
        row_ensembl.append(stable)
    numeric = df.iloc[:, ann_count:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    z = zscore_rows(numeric)
    sample_subset = sample_map[sample_map["dataset_id"] == dataset_id].copy()
    case_by_sample = dict(zip(sample_subset["sample_id"].astype(str), sample_subset["inferred_case_id"].astype(str)))
    pattern_by_sample = dict(zip(sample_subset["sample_id"].astype(str), sample_subset["sample_id_pattern"].astype(str)))
    score_rows: list[dict[str, Any]] = []
    recovery_rows: list[dict[str, Any]] = []
    requested_by_signature = signature_table.groupby("signature")["gene_symbol"].nunique().to_dict()
    for signature in SIGNATURES:
        mask = np.array([signature in s for s in row_to_signatures], dtype=bool)
        feature_count = int(mask.sum())
        if feature_count == 0:
            unique_genes = 0
            recovered_symbols: set[str] = set()
            recovered_ensembl: set[str] = set()
        else:
            recovered_symbols = set()
            for i, flag in enumerate(mask):
                if flag:
                    recovered_symbols.update(row_to_signature_genes[i].get(signature, set()))
            recovered_ensembl = {row_ensembl[i] for i, flag in enumerate(mask) if flag and row_ensembl[i]}
            unique_genes = len(recovered_symbols)
        requested = int(requested_by_signature.get(signature, 0))
        recovery_fraction = unique_genes / requested if requested else np.nan
        recovery_rows.append(
            {
                "dataset_id": dataset_id,
                "source": source,
                "cancer": cancer,
                "datatype": datatype,
                "tumor_normal_hint": tumor_normal_hint,
                "signature": signature,
                "requested_gene_count": requested,
                "recovered_gene_or_ensembl_count": unique_genes,
                "recovered_feature_count": feature_count,
                "recovery_fraction": recovery_fraction,
                "recovered_gene_symbols": ";".join(sorted(recovered_symbols)),
                "recovered_ensembl_ids": ";".join(sorted(recovered_ensembl)),
                "scoreable": int(feature_count > 0),
            }
        )
        if feature_count == 0:
            continue
        scores = np.nanmean(z[mask, :], axis=0)
        for col_idx, sample_id in enumerate(sample_ids):
            sid = safe_str(sample_id)
            score = scores[col_idx] if col_idx < len(scores) else np.nan
            score_rows.append(
                {
                    "dataset_id": dataset_id,
                    "source": source,
                    "cancer": cancer,
                    "datatype": datatype,
                    "tumor_normal_hint": tumor_normal_hint,
                    "signature": signature,
                    "sample_id": sid,
                    "inferred_case_id": case_by_sample.get(sid, infer_cptac_case_id(sid)),
                    "sample_id_pattern": pattern_by_sample.get(sid, ""),
                    "signature_score": score,
                    "n_features_used": feature_count,
                    "n_recovered_gene_or_ensembl": unique_genes,
                    "recovery_fraction": recovery_fraction,
                    "is_main_discovery_stratum": int(tumor_normal_hint in {"tumor", "unknown"}),
                }
            )
    return score_rows, recovery_rows


def cibersort_context(df: pd.DataFrame) -> pd.DataFrame:
    sample_col = "Input Sample"
    numeric = df.copy()
    for col in numeric.columns:
        if col != sample_col:
            numeric[col] = pd.to_numeric(numeric[col], errors="coerce")
    t_cols = [c for c in numeric.columns if c.startswith("T cells")]
    myeloid_cols = [c for c in numeric.columns if c.startswith(("Macrophages", "Monocytes", "Dendritic cells"))]
    cytotoxic_cols = [c for c in numeric.columns if c in ["T cells CD8", "NK cells resting", "NK cells activated"]]
    out = pd.DataFrame(
        {
            "sample_id": numeric[sample_col].astype(str),
            "inferred_case_id": numeric[sample_col].map(infer_cptac_case_id),
            "cibersort_t_cell": numeric[t_cols].sum(axis=1, skipna=True) if t_cols else np.nan,
            "cibersort_myeloid": numeric[myeloid_cols].sum(axis=1, skipna=True) if myeloid_cols else np.nan,
            "cibersort_cytotoxic": numeric[cytotoxic_cols].sum(axis=1, skipna=True) if cytotoxic_cols else np.nan,
            "cibersort_absolute_score": numeric.get("Absolute score", np.nan),
        }
    )
    return out


def xcell_context(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    row_names = pd.Series(df.index.astype(str), index=df.index)
    t_rows = row_names[row_names.str.contains("T-cells|Tcm|Tem|Tregs|Th1|Th2", case=False, regex=True)].index
    myeloid_rows = row_names[row_names.str.contains("Macrophages|Monocytes|\\bDC\\b|iDC|aDC|cDC", case=False, regex=True)].index
    stroma_rows = row_names[row_names.str.contains("Fibroblasts|Endothelial|MSC|Pericytes", case=False, regex=True)].index
    out = pd.DataFrame({"sample_id": df.columns.astype(str)})
    out["inferred_case_id"] = out["sample_id"].map(infer_cptac_case_id)
    out["xcell_t_cell"] = df.loc[t_rows].mean(axis=0).values if len(t_rows) else np.nan
    out["xcell_myeloid"] = df.loc[myeloid_rows].mean(axis=0).values if len(myeloid_rows) else np.nan
    out["xcell_stroma"] = df.loc[stroma_rows].mean(axis=0).values if len(stroma_rows) else np.nan
    return out


def build_immune_context(root: Path, cptac_inv: pd.DataFrame) -> pd.DataFrame:
    contexts: list[pd.DataFrame] = []
    for _, item in cptac_inv[cptac_inv["source"] == "washu"].iterrows():
        datatype = safe_str(item["datatype"])
        cancer = safe_str(item["cancer"])
        path = resolve_path(safe_str(item["local_path"]), root)
        if datatype == "cibersort":
            df = pd.read_csv(path, sep="\t")
            ctx = cibersort_context(df)
        elif datatype == "xcell":
            ctx = xcell_context(path)
        elif datatype == "tumor_purity":
            df = pd.read_csv(path, sep="\t")
            ctx = pd.DataFrame(
                {
                    "sample_id": df["Sample_ID"].astype(str),
                    "inferred_case_id": df["Sample_ID"].map(infer_cptac_case_id),
                    "estimate_stromal_score": pd.to_numeric(df.get("StromalScore"), errors="coerce"),
                    "estimate_immune_score": pd.to_numeric(df.get("ImmuneScore"), errors="coerce"),
                    "estimate_score": pd.to_numeric(df.get("ESTIMATEScore"), errors="coerce"),
                    "estimate_tumor_purity": pd.to_numeric(df.get("TumorPurity"), errors="coerce"),
                }
            )
        else:
            continue
        ctx["context_source"] = datatype
        ctx["context_cancer"] = cancer
        contexts.append(ctx)
    if not contexts:
        return pd.DataFrame(columns=["inferred_case_id"])
    stacked = pd.concat(contexts, ignore_index=True, sort=False)
    numeric_cols = [c for c in CONTEXT_VARS + ["estimate_score", "cibersort_absolute_score"] if c in stacked.columns]
    grouped = stacked.groupby("inferred_case_id", dropna=False)
    out = grouped[numeric_cols].mean(numeric_only=True).reset_index()
    out["context_sources"] = grouped["context_source"].agg(lambda x: ";".join(sorted(set(map(str, x))))).values
    out["context_cancers"] = grouped["context_cancer"].agg(lambda x: ";".join(sorted(set(map(str, x))))).values
    out["context_sample_rows"] = grouped.size().values
    return out[out["inferred_case_id"].astype(str) != ""].reset_index(drop=True)


def spearman_table(scores: pd.DataFrame, context: pd.DataFrame, min_n: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = scores.merge(context, on="inferred_case_id", how="left", validate="many_to_one")
    rows = []
    group_cols = ["source", "cancer", "datatype", "tumor_normal_hint", "signature"]
    for keys, group in merged.groupby(group_cols, dropna=False):
        source, cancer, datatype, tumor_normal_hint, signature = keys
        for var in CONTEXT_VARS:
            if var not in group.columns:
                continue
            sub = group[["signature_score", var]].dropna()
            if len(sub) < min_n or sub[var].nunique() < 3 or sub["signature_score"].nunique() < 3:
                continue
            rho, p = spearmanr(sub["signature_score"], sub[var])
            rows.append(
                {
                    "source": source,
                    "cancer": cancer,
                    "datatype": datatype,
                    "tumor_normal_hint": tumor_normal_hint,
                    "signature": signature,
                    "context_variable": var,
                    "n": len(sub),
                    "spearman_rho": rho,
                    "p_value": p,
                    "p_family": "all_signature_context_correlations_v0_1_exploratory",
                }
            )
    corr = pd.DataFrame(rows)
    if not corr.empty:
        corr["fdr_bh_exploratory"] = bh_fdr(corr["p_value"].to_numpy())
    return corr, merged


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


def residualize_scores(merged: pd.DataFrame, min_n: int) -> pd.DataFrame:
    rows = []
    group_cols = ["source", "cancer", "datatype", "tumor_normal_hint", "signature"]
    priority_covars = ["estimate_tumor_purity", "estimate_immune_score", "cibersort_t_cell", "cibersort_myeloid"]
    for keys, group in merged.groupby(group_cols, dropna=False):
        source, cancer, datatype, tumor_normal_hint, signature = keys
        available = []
        for covar in priority_covars:
            if covar in group.columns and group[covar].notna().sum() >= min_n and group[covar].nunique(dropna=True) >= 3:
                available.append(covar)
        available = available[:3]
        if not available:
            continue
        model_df = group[
            [
                "dataset_id",
                "sample_id",
                "inferred_case_id",
                "signature_score",
                "n_features_used",
                "n_recovered_gene_or_ensembl",
            ]
            + available
        ].dropna()
        if len(model_df) < max(min_n, 10 + 5 * len(available)):
            continue
        y = model_df["signature_score"].to_numpy(dtype=float)
        x = model_df[available].to_numpy(dtype=float)
        x = (x - np.nanmean(x, axis=0)) / np.nanstd(x, axis=0, ddof=1)
        x = np.column_stack([np.ones(len(x)), x])
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        residual = y - x @ beta
        for i, (_, row) in enumerate(model_df.iterrows()):
            rows.append(
                {
                    "source": source,
                    "cancer": cancer,
                    "datatype": datatype,
                    "tumor_normal_hint": tumor_normal_hint,
                    "signature": signature,
                    "dataset_id": row["dataset_id"],
                    "sample_id": row["sample_id"],
                    "inferred_case_id": row["inferred_case_id"],
                    "signature_score": row["signature_score"],
                    "residualized_signature_score": residual[i],
                    "residualization_covariates": ";".join(available),
                    "residualization_n": len(model_df),
                    "n_features_used": row["n_features_used"],
                    "n_recovered_gene_or_ensembl": row["n_recovered_gene_or_ensembl"],
                }
            )
    return pd.DataFrame(rows)


def cross_stratum_consistency(recovery: pd.DataFrame, corr: pd.DataFrame, residualized: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for signature, group in recovery.groupby("signature"):
        main = group[group["tumor_normal_hint"].isin(["tumor", "unknown"])]
        scoreable = main[main["scoreable"] == 1]
        row = {
            "signature": signature,
            "scoreable_strata": int(len(scoreable)),
            "scoreable_cancers": int(scoreable["cancer"].nunique()),
            "scoreable_sources": ";".join(sorted(set(scoreable["source"].astype(str)))),
            "scoreable_datatypes": ";".join(sorted(set(scoreable["datatype"].astype(str)))),
            "median_recovery_fraction": float(scoreable["recovery_fraction"].median()) if len(scoreable) else np.nan,
            "median_recovered_features": float(scoreable["recovered_feature_count"].median()) if len(scoreable) else np.nan,
        }
        if not corr.empty:
            sub = corr[(corr["signature"] == signature) & (corr["context_variable"] == "estimate_immune_score")]
            row["immune_score_corr_strata"] = int(len(sub))
            row["immune_score_corr_positive"] = int((sub["spearman_rho"] > 0).sum()) if len(sub) else 0
            row["immune_score_corr_negative"] = int((sub["spearman_rho"] < 0).sum()) if len(sub) else 0
            row["median_immune_score_rho"] = float(sub["spearman_rho"].median()) if len(sub) else np.nan
            subp = corr[(corr["signature"] == signature) & (corr["context_variable"] == "estimate_tumor_purity")]
            row["tumor_purity_corr_strata"] = int(len(subp))
            row["median_tumor_purity_rho"] = float(subp["spearman_rho"].median()) if len(subp) else np.nan
        if not residualized.empty:
            row["residualized_strata"] = int(
                residualized[residualized["signature"] == signature][
                    ["source", "cancer", "datatype", "tumor_normal_hint"]
                ].drop_duplicates().shape[0]
            )
        else:
            row["residualized_strata"] = 0
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["scoreable_strata", "median_recovery_fraction"], ascending=False)


def claim_tier_table(consistency: pd.DataFrame, residualized: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in consistency.iterrows():
        scoreable = int(row.get("scoreable_strata", 0))
        residualized_n = int(row.get("residualized_strata", 0))
        if scoreable >= 20 and residualized_n >= 10:
            tier = "analysis_ready_context_adjusted"
            claim = "Can be carried into v0.2 context-adjusted discovery; still not a validation or mechanism claim."
        elif scoreable >= 10:
            tier = "analysis_ready_unadjusted"
            claim = "Can be described as scoreable/recoverable in CPTAC strata; context adjustment may be incomplete."
        else:
            tier = "low_recovery_or_limited_support"
            claim = "Use only as exploratory or drop from main discovery."
        rows.append(
            {
                "signature": row["signature"],
                "claim_tier": tier,
                "scoreable_strata": scoreable,
                "residualized_strata": residualized_n,
                "allowed_claim": claim,
                "forbidden_claim": "Do not write therapeutic target, validated biomarker, causal immune-evasion mechanism, or independent validation.",
            }
        )
    return pd.DataFrame(rows)


def write_figures(recovery: pd.DataFrame, corr: pd.DataFrame, fig_dir: Path) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    made = []
    main = recovery[recovery["tumor_normal_hint"].isin(["tumor", "unknown"])].copy()
    main["source_datatype"] = main["source"] + "_" + main["datatype"]
    pivot = main.pivot_table(
        index="signature",
        columns="source_datatype",
        values="recovery_fraction",
        aggfunc="median",
    )
    if not pivot.empty:
        fig, ax = plt.subplots(figsize=(8, 4.8))
        data = pivot.fillna(0).to_numpy()
        im = ax.imshow(data, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=8)
        ax.set_title("Signature gene recovery across CPTAC source/datatype strata")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Median recovery fraction")
        fig.tight_layout()
        out = fig_dir / "fig01_signature_recovery_heatmap_v0_1.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        made.append(str(out))
    if not corr.empty:
        subset = corr[corr["context_variable"].isin(["estimate_immune_score", "estimate_tumor_purity"])].copy()
        if not subset.empty:
            plot = subset.groupby(["signature", "context_variable"])["spearman_rho"].median().reset_index()
            signatures = list(plot["signature"].drop_duplicates())
            x = np.arange(len(signatures))
            width = 0.35
            fig, ax = plt.subplots(figsize=(9, 4.8))
            for i, var in enumerate(["estimate_immune_score", "estimate_tumor_purity"]):
                vals = [
                    plot[(plot["signature"] == sig) & (plot["context_variable"] == var)]["spearman_rho"].median()
                    for sig in signatures
                ]
                vals = [0 if pd.isna(v) else v for v in vals]
                ax.bar(x + (i - 0.5) * width, vals, width=width, label=var)
            ax.axhline(0, color="black", linewidth=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(signatures, rotation=35, ha="right", fontsize=8)
            ax.set_ylabel("Median Spearman rho across scoreable strata")
            ax.set_title("Signature-score association with immune/purity context")
            ax.legend(fontsize=8)
            fig.tight_layout()
            out = fig_dir / "fig02_signature_context_correlation_summary_v0_1.png"
            fig.savefig(out, dpi=300)
            plt.close(fig)
            made.append(str(out))
    return made


def write_report(report_path: Path, summary: dict[str, Any], consistency: pd.DataFrame, claim: pd.DataFrame) -> None:
    lines = [
        "# Pan-cancer CPTAC proteomic/phosphoproteomic immune-evasion discovery v0.1",
        "",
        f"Generated: {summary['generated_at']}",
        "",
        "## Scope",
        "",
        "- First-pass discovery only.",
        "- CPTAC BCM and UMich were kept separate.",
        "- Tumor and unknown tumor-normal strata were scored; normal matrices were not used for main discovery.",
        "- WashU CIBERSORT, xCell, and ESTIMATE tumor purity were used only as immune/purity context.",
        "- No survival analysis, clinical outcome model, therapeutic target claim, validation claim, scRNA mapping, or TCGA-CPTAC sample merge was performed.",
        "",
        "## Output Counts",
        "",
    ]
    for key in [
        "signature_gene_set_rows",
        "matrix_strata_scored",
        "signature_recovery_rows",
        "signature_scores_long_rows",
        "immune_context_rows",
        "context_correlation_rows",
        "residualized_score_rows",
        "cross_stratum_consistency_rows",
        "claim_tier_rows",
    ]:
        lines.append(f"- {key}: {summary.get(key)}")
    lines.extend(["", "## Strongest Recoverability Signals", ""])
    if consistency.empty:
        lines.append("- No consistency table generated.")
    else:
        for _, row in consistency.head(8).iterrows():
            lines.append(
                f"- `{row['signature']}`: scoreable strata={int(row['scoreable_strata'])}, "
                f"median recovery={row['median_recovery_fraction']:.3f}, "
                f"residualized strata={int(row.get('residualized_strata', 0))}"
            )
    lines.extend(["", "## Claim Tier Summary", ""])
    if not claim.empty:
        counts = claim["claim_tier"].value_counts().to_dict()
        for tier, n in counts.items():
            lines.append(f"- `{tier}`: {n} signatures")
    lines.extend(
        [
            "",
            "## Required Interpretation Boundary",
            "",
            "- These results can support selection of signatures and strata for v0.2 discovery.",
            "- They cannot yet support a validated immune-evasion state, therapeutic vulnerability, or mechanism claim.",
            "- Any future RNA-protein discordance claim requires a dedicated RNA comparator module and careful platform-level wording.",
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
    out_dir = resolve_path(args.out_dir, root)
    fig_dir = resolve_path(args.figure_dir, root)
    report_path = resolve_path(args.report, root)
    summary_path = resolve_path(args.summary, root)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cptac_inv = pd.read_csv(qc_dir / "cptac_processed_inventory_v0_1.csv")
    sample_map = pd.read_csv(harmonization_dir / "cptac_sample_map_v0_1.csv")

    signature_table, signature_symbol_map, signature_ensembl_map = build_signature_table(
        harmonization_dir, out_dir / "signature_gene_sets_v0_1.csv"
    )
    discovery_mats = cptac_inv[
        cptac_inv["source"].isin(["bcm", "umich"])
        & cptac_inv["datatype"].isin(["proteomics", "phosphoproteomics"])
        & cptac_inv["tumor_normal_hint"].isin(["tumor", "unknown"])
    ].copy()
    all_scores: list[dict[str, Any]] = []
    all_recovery: list[dict[str, Any]] = []
    for _, item in discovery_mats.iterrows():
        scores, recovery = score_one_matrix(
            root,
            item,
            sample_map,
            signature_table,
            signature_symbol_map,
            signature_ensembl_map,
        )
        all_scores.extend(scores)
        all_recovery.extend(recovery)
        print(f"scored {item['dataset_id']}: scores={len(scores)} recovery={len(recovery)}", flush=True)
    recovery_df = pd.DataFrame(all_recovery)
    scores_long = pd.DataFrame(all_scores)
    recovery_df.to_csv(out_dir / "cptac_signature_recovery_v0_1.csv", index=False)
    scores_long.to_csv(out_dir / "cptac_signature_scores_long_v0_1.csv", index=False)
    if not scores_long.empty:
        wide_index = ["dataset_id", "source", "cancer", "datatype", "tumor_normal_hint", "sample_id", "inferred_case_id"]
        scores_wide = scores_long.pivot_table(
            index=wide_index,
            columns="signature",
            values="signature_score",
            aggfunc="mean",
        ).reset_index()
        scores_wide.columns.name = None
    else:
        scores_wide = pd.DataFrame()
    scores_wide.to_csv(out_dir / "cptac_signature_scores_wide_v0_1.csv", index=False)

    context = build_immune_context(root, cptac_inv)
    context.to_csv(out_dir / "cptac_immune_context_scores_v0_1.csv", index=False)
    corr, merged = spearman_table(scores_long, context, args.min_correlation_n)
    corr.to_csv(out_dir / "cptac_signature_context_correlations_v0_1.csv", index=False)
    residualized = residualize_scores(merged, args.min_residual_n)
    residualized.to_csv(out_dir / "cptac_signature_residualized_scores_v0_1.csv", index=False)
    consistency = cross_stratum_consistency(recovery_df, corr, residualized)
    consistency.to_csv(out_dir / "cptac_cross_stratum_consistency_v0_1.csv", index=False)
    claim = claim_tier_table(consistency, residualized)
    claim.to_csv(out_dir / "discovery_claim_tier_table_v0_1.csv", index=False)
    figures = write_figures(recovery_df, corr, fig_dir)

    summary = {
        "version": args.version,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "project_root": str(root),
        "matrix_strata_scored": int(len(discovery_mats)),
        "signature_gene_set_rows": int(len(signature_table)),
        "signature_recovery_rows": int(len(recovery_df)),
        "signature_scores_long_rows": int(len(scores_long)),
        "signature_scores_wide_rows": int(len(scores_wide)),
        "immune_context_rows": int(len(context)),
        "context_correlation_rows": int(len(corr)),
        "residualized_score_rows": int(len(residualized)),
        "cross_stratum_consistency_rows": int(len(consistency)),
        "claim_tier_rows": int(len(claim)),
        "figures": [rel(Path(f), root) for f in figures],
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "elapsed_seconds": round(time.time() - start, 2),
        "outputs": {
            "signature_gene_sets": rel(out_dir / "signature_gene_sets_v0_1.csv", root),
            "signature_recovery": rel(out_dir / "cptac_signature_recovery_v0_1.csv", root),
            "signature_scores_long": rel(out_dir / "cptac_signature_scores_long_v0_1.csv", root),
            "signature_scores_wide": rel(out_dir / "cptac_signature_scores_wide_v0_1.csv", root),
            "immune_context_scores": rel(out_dir / "cptac_immune_context_scores_v0_1.csv", root),
            "context_correlations": rel(out_dir / "cptac_signature_context_correlations_v0_1.csv", root),
            "residualized_scores": rel(out_dir / "cptac_signature_residualized_scores_v0_1.csv", root),
            "cross_stratum_consistency": rel(out_dir / "cptac_cross_stratum_consistency_v0_1.csv", root),
            "claim_tier_table": rel(out_dir / "discovery_claim_tier_table_v0_1.csv", root),
        },
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    write_report(report_path, summary, consistency, claim)
    print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
