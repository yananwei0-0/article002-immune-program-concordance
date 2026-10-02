#!/usr/bin/env python3
"""Memory-bounded real-source H5AD execution for BIO-05 single-cell and spatial support.

The public CELLxGENE files contain large ``uns`` objects and, in several files,
both ``X`` and ``raw.X``.  This implementation reads only required obs/var
fields and scans the selected CSR matrix once in bounded row chunks.  It keeps
the locked donor-compartment and sample-level inference units.
"""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import json
import math
import os
import re
import shutil
import time
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", expand_legacy_paths("${SCI27_TMP_ROOT}/bio05_mpl_cache"))

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROGRAM_LABELS = {
    "antigen_presentation_mhc_i": "MHC-I",
    "antigen_presentation_mhc_ii": "MHC-II",
    "ifn_gamma_response": "IFN-gamma",
    "cytolytic_t_cell_context": "Cytolytic",
    "checkpoint_exhaustion_context": "Checkpoint",
    "myeloid_inflammatory_context": "Myeloid",
    "tgfb_emt_exclusion_context": "TGF/EMT",
    "proteasome_antigen_processing": "Proteasome",
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
    "immune cell", "b cell", "plasma cell", "t cell", "nk", "natural killer",
    "macrophage", "monocyte", "microgl", "dendritic", "myeloid", "mast cell",
    "neutrophil", "granulocyte", "lymphocyte", "megakaryocyte",
}
STROMAL_KEYWORDS = {
    "fibroblast", "myofibroblast", "endothelial", "pericyte", "smooth muscle",
    "mural cell", "vascular associated smooth muscle cell", "mesangial", "podocyte",
    "schwann", "stromal", "astrocyte", "oligodendrocyte",
    "oligodendrocyte precursor", "radial glial", "neuron", "support cell",
}
MALIGNANT_KEYWORDS = {"malignant", "abnormal cell", "tumor", "tumour", "cancer epithelial", "neoplastic"}
DATASET_SPECIFIC_MALIGNANT = {
    "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": {"squamous epithelial cell"},
    "0bebef1a-4607-4584-9070-dacf89a0d635": {"epithelial cell of lung"},
}
DATASET_SPECIFIC_OTHER = {
    "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": {"salivary gland glandular cell"},
}
CANCER_SCORE_COLS = ["Cancer Basal SC", "Cancer Cycling", "Cancer Her2 SC", "Cancer LumA SC", "Cancer LumB SC"]
IMMUNE_SCORE_COLS = [
    "B cells Memory", "B cells Naive", "Cycling T-cells", "Cycling_Myeloid", "DCs",
    "Macrophage", "Monocyte", "NK cells", "NKT cells", "Plasmablasts",
    "T cells CD4+", "T cells CD8+",
]
STROMAL_SCORE_COLS = [
    "CAFs MSC iCAF-like", "CAFs myCAF-like", "Endothelial ACKR1", "Endothelial CXCL12",
    "Endothelial Lymphatic LYVE1", "Endothelial RGS5", "Myoepithelial",
    "PVL Differentiated", "PVL Immature",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["scrna", "spatial", "all"], default="all")
    parser.add_argument("--scrna-queue", required=True)
    parser.add_argument("--spatial-queue", required=True)
    parser.add_argument("--signature-table", required=True)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument(
        "--input-root",
        default=".",
        help="Base directory used to resolve relative remote_path values in the queue TSV files.",
    )
    parser.add_argument("--chunk-rows", type=int, default=2000)
    parser.add_argument("--min-group-cells", type=int, default=20)
    parser.add_argument("--min-signature-genes", type=int, default=3)
    parser.add_argument(
        "--source-cache-dir",
        default="",
        help="Optional internal scratch directory; source bytes are copied and size-verified before H5AD reading.",
    )
    return parser.parse_args()


def decode_array(array: Any) -> np.ndarray:
    values = np.asarray(array)
    if values.dtype.kind in {"S", "O", "U"}:
        return np.asarray([x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values], dtype=object)
    return values


def read_node(node: h5py.Dataset | h5py.Group) -> np.ndarray:
    if isinstance(node, h5py.Dataset):
        return decode_array(node[()])
    encoding = node.attrs.get("encoding-type", b"")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding == "categorical":
        categories = read_node(node["categories"])
        codes = np.asarray(node["codes"][()], dtype=int)
        out = np.empty(len(codes), dtype=object)
        out[:] = None
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    if encoding in {"nullable-integer", "nullable-boolean"}:
        values = np.asarray(node["values"][()])
        mask = np.asarray(node["mask"][()], dtype=bool)
        out = values.astype(object)
        out[mask] = None
        return out
    raise RuntimeError(f"Unsupported H5AD node encoding: {encoding} at {node.name}")


def column_order(group: h5py.Group) -> list[str]:
    raw = group.attrs.get("column-order", [])
    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def read_frame_columns(group: h5py.Group, requested: Iterable[str]) -> pd.DataFrame:
    available = set(column_order(group))
    index_name = group.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    index = read_node(group[str(index_name)])
    data: dict[str, Any] = {}
    for name in requested:
        if name in available and name in group:
            data[name] = read_node(group[name])
    return pd.DataFrame(data, index=pd.Index(index.astype(str), name="_index"))


def load_signature_genes(path: Path) -> dict[str, list[str]]:
    table = pd.read_csv(path)
    return {
        str(program): sorted(set(group["gene_symbol"].dropna().astype(str).str.upper()))
        for program, group in table.groupby("signature", sort=False)
    }


def norm_text(value: Any) -> str:
    return re.sub(r"\s+", " ", "" if value is None else str(value).strip().lower())


def norm_series(values: pd.Series) -> pd.Series:
    """Vectorised equivalent of :func:`norm_text` for large obs frames."""
    return values.fillna("").astype(str).str.strip().str.lower().str.replace(r"\s+", " ", regex=True)


def infer_cancer_frame(disease: pd.Series, tissue: pd.Series, default: str | None) -> pd.Series:
    """Apply the locked disease/tissue mapping without per-row pandas objects."""
    d = norm_series(disease)
    t = norm_series(tissue)
    text = d + " | " + t
    result = pd.Series(pd.NA, index=d.index, dtype="object")
    normal = d.isin(["normal", "healthy"])

    def assign(pattern: str, code: str) -> None:
        mask = ~normal & result.isna() & text.str.contains(pattern, regex=True, na=False)
        result.loc[mask] = code

    assign(r"(?=.*squamous cell carcinoma)(?=.*lung)", "lscc")
    assign(r"clear cell renal|renal cell carcinoma|kidney cancer", "ccrcc")
    assign(r"glioblastoma", "gbm")
    assign(r"head and neck|oropharynx squamous|oral cavity squamous", "hnscc")
    assign(r"lung adenocarcinoma", "luad")
    assign(r"squamous cell lung", "lscc")
    assign(r"breast", "brca")
    assign(r"ovarian|serous adenocarcinoma", "ov")
    assign(r"colon adenocarcinoma|colorectal|cecum adenocarcinoma|rectosigmoid", "coad")
    if default is not None:
        result.loc[~normal & result.isna()] = default
    return result


def classify_compartment_frame(obs: pd.DataFrame, source_id: str) -> tuple[pd.Series, pd.Series]:
    """Vectorised, priority-preserving form of ``classify_compartment``."""
    blank = pd.Series("", index=obs.index, dtype="object")
    cell_type = norm_series(obs.get("cell_type", blank))
    major = norm_series(obs.get("celltype_major", blank))
    author = norm_series(obs.get("author_cell_type", blank))
    cluster = norm_series(obs.get("author_cluster_label", blank))
    normal_call = norm_series(obs.get("normal_cell_call", blank))
    combined = cell_type + " | " + major + " | " + author + " | " + cluster
    compartment = pd.Series("other", index=obs.index, dtype="object")
    rule = pd.Series("unmapped", index=obs.index, dtype="object")
    assigned = pd.Series(False, index=obs.index)

    def apply(mask: pd.Series, label: str, reason: str) -> None:
        take = mask.fillna(False) & ~assigned
        compartment.loc[take] = label
        rule.loc[take] = reason
        assigned.loc[take] = True

    apply(normal_call.eq("cancer"), "malignant", "explicit_normal_cell_call_cancer")
    apply(normal_call.eq("normal") & combined.str.contains("epithelial", regex=False), "other", "explicit_normal_cell_call_normal_epithelial")
    apply(major.eq("cancer epithelial"), "malignant", "explicit_celltype_major_cancer_epithelial")
    apply(major.isin(["t-cells", "myeloid", "plasmablasts", "b-cells"]), "immune", "explicit_celltype_major_immune")
    apply(major.isin(["cafs", "endothelial", "pvl"]), "stromal", "explicit_celltype_major_stromal")
    apply(major.eq("normal epithelial"), "other", "explicit_celltype_major_normal_epithelial")
    specific_other = {norm_text(x) for x in DATASET_SPECIFIC_OTHER.get(source_id, set())}
    specific_malignant = {norm_text(x) for x in DATASET_SPECIFIC_MALIGNANT.get(source_id, set())}
    apply(cell_type.isin(specific_other), "other", "dataset_specific_other")
    apply(cell_type.isin(specific_malignant), "malignant", "dataset_specific_malignant")
    malignant_pattern = "|".join(re.escape(x) for x in sorted(MALIGNANT_KEYWORDS))
    immune_pattern = "|".join(re.escape(x) for x in sorted(IMMUNE_KEYWORDS))
    stromal_pattern = "|".join(re.escape(x) for x in sorted(STROMAL_KEYWORDS))
    apply(combined.str.contains(malignant_pattern, regex=True, na=False), "malignant", "keyword_malignant")
    apply(combined.str.contains(immune_pattern, regex=True, na=False), "immune", "keyword_immune")
    apply(combined.str.contains(stromal_pattern, regex=True, na=False), "stromal", "keyword_stromal")
    fallback_sources = {"dea97145-f712-431c-a223-6b5f565f362a", "4d82bd8e-8827-410b-96ae-2c806799c08c"}
    apply(cell_type.str.contains("epithelial", regex=False) & (source_id in fallback_sources), "malignant", "tumor_epithelial_fallback")
    return compartment, rule


def infer_cancer(disease: Any, tissue: Any, default: str | None) -> str | None:
    d, t = norm_text(disease), norm_text(tissue)
    if d in {"normal", "healthy"}:
        return None
    text = f"{d} | {t}"
    rules = [
        (("clear cell renal", "renal cell carcinoma", "kidney cancer"), "ccrcc"),
        (("glioblastoma",), "gbm"),
        (("head and neck", "oropharynx squamous", "oral cavity squamous"), "hnscc"),
        (("lung adenocarcinoma",), "luad"),
        (("squamous cell lung",), "lscc"),
        (("breast",), "brca"),
        (("ovarian", "serous adenocarcinoma"), "ov"),
        (("colon adenocarcinoma", "colorectal", "cecum adenocarcinoma", "rectosigmoid"), "coad"),
    ]
    if "squamous cell carcinoma" in text and "lung" in text:
        return "lscc"
    for terms, code in rules:
        if any(term in text for term in terms):
            return code
    return default


def classify_compartment(row: pd.Series, source_id: str) -> tuple[str, str]:
    cell_type = norm_text(row.get("cell_type"))
    major = norm_text(row.get("celltype_major"))
    author = norm_text(row.get("author_cell_type"))
    cluster = norm_text(row.get("author_cluster_label"))
    normal_call = norm_text(row.get("normal_cell_call"))
    combined = " | ".join([cell_type, major, author, cluster])
    if normal_call == "cancer":
        return "malignant", "explicit_normal_cell_call_cancer"
    if normal_call == "normal" and "epithelial" in combined:
        return "other", "explicit_normal_cell_call_normal_epithelial"
    if major == "cancer epithelial":
        return "malignant", "explicit_celltype_major_cancer_epithelial"
    if major in {"t-cells", "myeloid", "plasmablasts", "b-cells"}:
        return "immune", "explicit_celltype_major_immune"
    if major in {"cafs", "endothelial", "pvl"}:
        return "stromal", "explicit_celltype_major_stromal"
    if major == "normal epithelial":
        return "other", "explicit_celltype_major_normal_epithelial"
    if cell_type in {norm_text(x) for x in DATASET_SPECIFIC_OTHER.get(source_id, set())}:
        return "other", "dataset_specific_other"
    if cell_type in {norm_text(x) for x in DATASET_SPECIFIC_MALIGNANT.get(source_id, set())}:
        return "malignant", "dataset_specific_malignant"
    if any(term in combined for term in MALIGNANT_KEYWORDS):
        return "malignant", "keyword_malignant"
    if any(term in combined for term in IMMUNE_KEYWORDS):
        return "immune", "keyword_immune"
    if any(term in combined for term in STROMAL_KEYWORDS):
        return "stromal", "keyword_stromal"
    if "epithelial" in cell_type and source_id in {"dea97145-f712-431c-a223-6b5f565f362a", "4d82bd8e-8827-410b-96ae-2c806799c08c"}:
        return "malignant", "tumor_epithelial_fallback"
    return "other", "unmapped"


def pick_field(available: set[str], preferred: list[str], pattern: str) -> str | None:
    for name in preferred:
        if name in available:
            return name
    return next((x for x in available if re.search(pattern, x, flags=re.I)), None)


def gene_symbols(var_group: h5py.Group) -> np.ndarray:
    available = set(column_order(var_group))
    if "feature_name" in available:
        return np.char.upper(read_node(var_group["feature_name"]).astype(str))
    index_name = var_group.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return np.char.upper(read_node(var_group[str(index_name)]).astype(str))


def scan_group_csr(
    matrix: h5py.Group,
    cell_group_codes: np.ndarray,
    target_feature_indices: np.ndarray,
    group_count: int,
    chunk_rows: int,
) -> tuple[np.ndarray, np.ndarray]:
    shape = tuple(int(x) for x in matrix.attrs["shape"])
    if str(matrix.attrs.get("encoding-type", "")).replace("b'", "").replace("'", "") not in {"csr_matrix", "bcsr_matrix"}:
        encoding = matrix.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding != "csr_matrix":
            raise RuntimeError(f"Only CSR H5AD matrices are supported; observed {encoding}")
    print(f"  loading CSR indptr rows={shape[0]}", flush=True)
    indptr = np.asarray(matrix["indptr"][()], dtype=np.int64)
    print(f"  loaded CSR indptr entries={len(indptr)}", flush=True)
    lookup = np.full(shape[1], -1, dtype=np.int32)
    lookup[target_feature_indices] = np.arange(len(target_feature_indices), dtype=np.int32)
    group_total = np.zeros(group_count, dtype=np.float64)
    group_target = np.zeros((group_count, len(target_feature_indices)), dtype=np.float64)
    for row_start in range(0, shape[0], chunk_rows):
        row_end = min(row_start + chunk_rows, shape[0])
        value_start, value_end = int(indptr[row_start]), int(indptr[row_end])
        if value_end <= value_start:
            continue
        indices = np.asarray(matrix["indices"][value_start:value_end], dtype=np.int64)
        values = np.asarray(matrix["data"][value_start:value_end], dtype=np.float64)
        counts = np.diff(indptr[row_start : row_end + 1]).astype(np.int64)
        local_rows = np.repeat(np.arange(row_end - row_start, dtype=np.int32), counts)
        entry_groups = cell_group_codes[row_start:row_end][local_rows]
        eligible = entry_groups >= 0
        if eligible.any():
            group_total += np.bincount(
                entry_groups[eligible], weights=values[eligible], minlength=group_count
            )
            target_positions = lookup[indices]
            hit = eligible & (target_positions >= 0)
            if hit.any():
                flat = entry_groups[hit].astype(np.int64) * len(target_feature_indices) + target_positions[hit]
                group_target += np.bincount(
                    flat, weights=values[hit], minlength=group_count * len(target_feature_indices)
                ).reshape(group_count, len(target_feature_indices))
        if row_start == 0 or row_end == shape[0] or (row_start // chunk_rows) % 50 == 0:
            print(f"  CSR rows {row_end}/{shape[0]}", flush=True)
    return group_total, group_target


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    means = frame.mean(axis=0, skipna=True)
    std = frame.std(axis=0, ddof=0, skipna=True).replace(0, np.nan)
    return (frame - means) / std


def schema_row(source_id: str, title: str, role: str, field: str | None, values: pd.Series | None) -> dict[str, Any]:
    counts = values.astype(str).value_counts().head(10) if values is not None else pd.Series(dtype=int)
    return {
        "source_id": source_id,
        "dataset_title": title,
        "field_role": role,
        "field_name": field or "",
        "nunique": int(values.nunique()) if values is not None else 0,
        "top_values": "; ".join(f"{idx}:{int(val)}" for idx, val in counts.items()),
    }


def process_scrna(row: pd.Series, signature_genes: dict[str, list[str]], chunk_rows: int, min_cells: int, min_genes: int):
    source_id, title, path = str(row.source_id), str(row.title), Path(row.remote_path)
    started = time.time()
    print(f"START scRNA {source_id} {path.name}", flush=True)
    obs_fields = {
        "donor_id", "patient_id", "participant_id", "cell_type", "author_cell_type",
        "celltype_major", "author_cluster_label", "normal_cell_call", "disease", "tissue",
    }
    with h5py.File(path, "r") as handle:
        available = set(column_order(handle["obs"]))
        donor_col = pick_field(available, ["donor_id", "patient_id", "participant_id"], r"donor|patient|participant")
        cell_type_col = pick_field(available, ["cell_type", "author_cell_type"], r"cell.?type|cluster|annotation")
        disease_col = pick_field(available, ["disease"], r"disease|cancer|tumou?r")
        tissue_col = pick_field(available, ["tissue"], r"tissue|organ")
        requested = obs_fields | {x for x in [donor_col, cell_type_col, disease_col, tissue_col] if x}
        obs = read_frame_columns(handle["obs"], requested)
        print(f"  loaded obs fields rows={len(obs)} columns={len(obs.columns)}", flush=True)
        if not donor_col or not cell_type_col:
            schema = [schema_row(source_id, title, role, None, None) for role in ["donor", "cell_type", "disease", "tissue"]]
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(schema)
        obs["donor_id_resolved"] = obs[donor_col].astype(str)
        obs["disease_resolved"] = obs[disease_col].astype(str) if disease_col else ""
        obs["tissue_resolved"] = obs[tissue_col].astype(str) if tissue_col else ""
        obs["cancer_code"] = infer_cancer_frame(obs["disease_resolved"], obs["tissue_resolved"], DATASET_DEFAULT_CANCER.get(source_id))
        obs["broad_compartment"], obs["compartment_rule"] = classify_compartment_frame(obs, source_id)
        valid_donor = ~norm_series(obs["donor_id_resolved"]).isin(["", "nan", "none", "na", "<na>"])
        eligible = obs["cancer_code"].notna() & obs["broad_compartment"].isin(["malignant", "immune", "stromal"]) & valid_donor
        print(f"  mapped obs eligible_cells={int(eligible.sum())}", flush=True)
        obs["group_id"] = obs["cancer_code"].astype(str) + "__" + obs["donor_id_resolved"].astype(str) + "__" + obs["broad_compartment"].astype(str)
        counts = obs[eligible].groupby(["cancer_code", "donor_id_resolved", "broad_compartment"], dropna=False).size().reset_index(name="n_cells")
        rule_summary = obs[eligible].groupby(["cancer_code", "donor_id_resolved", "broad_compartment"])["compartment_rule"].agg(lambda x: ";".join(sorted(set(map(str, x))))).reset_index()
        counts = counts.merge(rule_summary, on=["cancer_code", "donor_id_resolved", "broad_compartment"], validate="one_to_one")
        counts["source_id"], counts["dataset_title"] = source_id, title
        keep = counts[counts["n_cells"] >= min_cells].copy()
        keep["group_id"] = keep["cancer_code"].astype(str) + "__" + keep["donor_id_resolved"].astype(str) + "__" + keep["broad_compartment"].astype(str)
        group_ids = keep["group_id"].tolist()
        group_map = {value: idx for idx, value in enumerate(group_ids)}
        cell_codes = np.full(len(obs), -1, dtype=np.int32)
        eligible_groups = obs["group_id"].isin(group_map)
        cell_codes[eligible_groups.to_numpy()] = obs.loc[eligible_groups, "group_id"].map(group_map).to_numpy(dtype=np.int32)
        matrix_name, var_name = ("raw/X", "raw/var") if "raw" in handle else ("X", "var")
        symbols = gene_symbols(handle[var_name])
        target_symbols = sorted({g for genes in signature_genes.values() for g in genes})
        target_idx = np.where(np.isin(symbols, target_symbols))[0].astype(np.int64)
        print(f"  retained_groups={len(group_ids)} target_features={len(target_idx)} matrix={matrix_name}", flush=True)
        if len(target_idx) == 0 or len(group_ids) == 0:
            return counts, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        group_total, group_target = scan_group_csr(handle[matrix_name], cell_codes, target_idx, len(group_ids), chunk_rows)
        with np.errstate(divide="ignore", invalid="ignore"):
            cpm = np.divide(group_target * 1e6, group_total[:, None], out=np.zeros_like(group_target), where=group_total[:, None] > 0)
            log_cpm = np.log1p(cpm)
        expr = pd.DataFrame(log_cpm, index=group_ids, columns=symbols[target_idx])
        expr = expr.T.groupby(level=0, sort=True).mean().T

    score_frames: list[pd.DataFrame] = []
    for cancer, group_meta in keep.groupby("cancer_code", sort=False):
        ids = group_meta["group_id"].tolist()
        cancer_z = zscore_columns(expr.loc[ids])
        meta = group_meta.set_index("group_id").loc[ids]
        for program, genes in signature_genes.items():
            recovered = [gene for gene in genes if gene in cancer_z.columns]
            if len(recovered) < min_genes:
                continue
            score_frames.append(
                pd.DataFrame(
                    {
                        "group_id": ids,
                        "signature": program,
                        "signature_label": PROGRAM_LABELS[program],
                        "signature_score": cancer_z[recovered].mean(axis=1).to_numpy(dtype=float),
                        "n_genes_used": len(recovered),
                        "source_id": source_id,
                        "dataset_title": title,
                        "cancer_code": cancer,
                        "donor_id": meta["donor_id_resolved"].astype(str).to_numpy(),
                        "compartment": meta["broad_compartment"].astype(str).to_numpy(),
                    }
                )
            )
    scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    pair_rows: list[dict[str, Any]] = []
    if not scores.empty:
        for keys, group in scores.groupby(["source_id", "dataset_title", "cancer_code", "signature", "signature_label"], sort=False):
            sid, dtitle, cancer, program, label = keys
            wide = group.pivot(index="donor_id", columns="compartment", values="signature_score")
            for left, right in [("malignant", "immune"), ("malignant", "stromal"), ("immune", "stromal")]:
                if left not in wide or right not in wide:
                    continue
                delta = (wide[left] - wide[right]).dropna()
                if delta.empty:
                    continue
                pair_rows.append({"source_id": sid, "dataset_title": dtitle, "cancer_code": cancer, "signature": program, "signature_label": label, "contrast": f"{left}_minus_{right}", "paired_donor_n": len(delta), "median_delta": delta.median(), "mean_delta": delta.mean(), "positive_fraction": (delta > 0).mean()})
    schema = pd.DataFrame(
        [
            schema_row(source_id, title, "donor", donor_col, obs[donor_col]),
            schema_row(source_id, title, "cell_type", cell_type_col, obs[cell_type_col]),
            schema_row(source_id, title, "disease", disease_col, obs[disease_col] if disease_col else None),
            schema_row(source_id, title, "tissue", tissue_col, obs[tissue_col] if tissue_col else None),
        ]
    )
    print(f"DONE scRNA {source_id} elapsed={time.time()-started:.1f}s groups={len(keep)} scores={len(scores)}", flush=True)
    return counts, scores, pd.DataFrame(pair_rows), schema


def assign_spot_compartment(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for label, cols in [("cancer_score", CANCER_SCORE_COLS), ("immune_score", IMMUNE_SCORE_COLS), ("stromal_score", STROMAL_SCORE_COLS)]:
        present = [col for col in cols if col in out]
        out[label] = out[present].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1) if present else 0.0
    compartments, rules = [], []
    for row in out.itertuples(index=False):
        values = {"malignant": float(row.cancer_score), "immune": float(row.immune_score), "stromal": float(row.stromal_score)}
        best = max(values, key=values.get)
        classification = norm_text(getattr(row, "Classification", ""))
        if values[best] > 0:
            compartments.append(best); rules.append("spot_score_argmax")
        elif classification == "stroma":
            compartments.append("stromal"); rules.append("classification_stroma")
        elif "lymphocytes" in classification:
            compartments.append("immune"); rules.append("classification_lymphocytes")
        elif "invasive cancer" in classification:
            compartments.append("malignant"); rules.append("classification_invasive_cancer")
        else:
            compartments.append("other"); rules.append("unclassified")
    out["spot_compartment"], out["compartment_rule"] = compartments, rules
    return out


def load_csr(matrix: h5py.Group) -> sparse.csr_matrix:
    return sparse.csr_matrix(
        (np.asarray(matrix["data"][()]), np.asarray(matrix["indices"][()], dtype=np.int32), np.asarray(matrix["indptr"][()], dtype=np.int64)),
        shape=tuple(int(x) for x in matrix.attrs["shape"]),
    )


def process_spatial(row: pd.Series, signature_genes: dict[str, list[str]], min_genes: int):
    source_id, title, path = str(row.source_id), str(row.title), Path(row.remote_path)
    print(f"START spatial {source_id} {path.name}", flush=True)
    fields = set(CANCER_SCORE_COLS + IMMUNE_SCORE_COLS + STROMAL_SCORE_COLS + ["in_tissue", "donor_id", "Classification", "annotation"])
    with h5py.File(path, "r") as handle:
        obs = read_frame_columns(handle["obs"], fields)
        obs["_row_index"] = np.arange(len(obs))
        if "in_tissue" in obs:
            obs = obs[pd.to_numeric(obs["in_tissue"], errors="coerce").fillna(0).astype(int).eq(1)].copy()
        obs = assign_spot_compartment(obs)
        obs = obs[obs["spot_compartment"].isin(["malignant", "immune", "stromal"])].copy()
        matrix = load_csr(handle["X"])[obs["_row_index"].to_numpy(), :]
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
        matrix = sparse.diags(np.divide(1e6, totals, out=np.zeros_like(totals), where=totals > 0)) @ matrix
        matrix.data = np.log1p(matrix.data)
        symbols = gene_symbols(handle["var"])
        spot_ids = obs.index.astype(str)
        sample_id = str(obs["donor_id"].iloc[0]) if "donor_id" in obs else title
        score_frames, summaries = [], []
        for program, genes in signature_genes.items():
            recovered_genes = [gene for gene in genes if gene in set(symbols)]
            if len(recovered_genes) < min_genes:
                continue
            idx = np.where(np.isin(symbols, recovered_genes))[0]
            dense = matrix[:, idx].toarray().astype(float)
            means = np.nanmean(dense, axis=0)
            std = np.nanstd(dense, axis=0, ddof=0)
            std[(std == 0) | ~np.isfinite(std)] = np.nan
            z = (dense - means) / std
            scores = np.nanmean(z, axis=1)
            frame = pd.DataFrame({"spot_id": spot_ids, "sample_id": sample_id, "spot_compartment": obs["spot_compartment"].to_numpy(), "compartment_rule": obs["compartment_rule"].to_numpy(), "Classification": obs.get("Classification", pd.Series(index=obs.index, dtype=object)).to_numpy(), "annotation": obs.get("annotation", pd.Series(index=obs.index, dtype=object)).to_numpy(), "signature": program, "signature_label": PROGRAM_LABELS[program], "signature_score": scores, "n_genes_used": len(recovered_genes), "source_id": source_id, "dataset_title": title})
            score_frames.append(frame)
            med = frame.groupby("spot_compartment")["signature_score"].median()
            summaries.append({"source_id": source_id, "dataset_title": title, "sample_id": sample_id, "signature": program, "signature_label": PROGRAM_LABELS[program], "n_spots": frame["spot_id"].nunique(), "n_genes_used": len(recovered_genes), "median_malignant": med.get("malignant"), "median_immune": med.get("immune"), "median_stromal": med.get("stromal")})
    counts = obs.groupby(["spot_compartment", "compartment_rule"]).size().reset_index(name="n_spots")
    counts["sample_id"], counts["source_id"], counts["dataset_title"] = sample_id, source_id, title
    print(f"DONE spatial {source_id} spots={len(obs)}", flush=True)
    return counts, pd.concat(score_frames, ignore_index=True), pd.DataFrame(summaries)


def write_scrna(
    root: Path,
    queue_path: Path,
    input_root: Path,
    genes: dict[str, list[str]],
    chunk_rows: int,
    min_cells: int,
    min_genes: int,
    source_cache_dir: Path | None,
) -> None:
    out = root / "results/scrna"; out.mkdir(parents=True, exist_ok=True)
    queue = pd.read_csv(queue_path, sep="\t", dtype=str)
    queue = queue[queue["source_type"].eq("cellxgene")].copy()
    queue["remote_path"] = queue["remote_path"].map(
        lambda value: str((input_root / value).resolve()) if not Path(value).is_absolute() else str(Path(value))
    )
    queue["observed_size_bytes"] = queue["remote_path"].map(lambda x: Path(x).stat().st_size if Path(x).exists() else 0)
    queue["is_complete"] = [int(Path(p).is_file() and (not str(e).isdigit() or Path(p).stat().st_size == int(e))) for p, e in zip(queue["remote_path"], queue["expected_size"])]
    queue.to_csv(out / "scrna_queue_audit_v0_7.csv", index=False)
    counts, scores, pairs, schemas = [], [], [], []
    for row in queue[queue["is_complete"].eq(1)].itertuples(index=False):
        source_row = pd.Series(row._asdict())
        if source_cache_dir is not None:
            source_cache_dir.mkdir(parents=True, exist_ok=True)
            original = Path(source_row["remote_path"])
            cached = source_cache_dir / f"{source_row['source_id']}{original.suffix}"
            expected = int(source_row["expected_size"]) if str(source_row["expected_size"]).isdigit() else original.stat().st_size
            if not cached.is_file() or cached.stat().st_size != expected:
                print(f"CACHE scRNA {source_row['source_id']} {original} -> {cached}", flush=True)
                shutil.copy2(original, cached)
            if cached.stat().st_size != expected:
                raise RuntimeError(f"Scratch copy size mismatch for {source_row['source_id']}")
            source_row["remote_path"] = str(cached)
            print(f"CACHE verified {source_row['source_id']} bytes={cached.stat().st_size}", flush=True)
        c, s, p, a = process_scrna(source_row, genes, chunk_rows, min_cells, min_genes)
        if not c.empty: counts.append(c)
        if not s.empty: scores.append(s)
        if not p.empty: pairs.append(p)
        if not a.empty: schemas.append(a)
    count_df = pd.concat(counts, ignore_index=True) if counts else pd.DataFrame()
    score_df = pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()
    pair_df = pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame()
    schema_df = pd.concat(schemas, ignore_index=True) if schemas else pd.DataFrame()
    summary = pd.DataFrame()
    if not score_df.empty:
        summary = score_df.groupby(["source_id", "dataset_title", "cancer_code", "signature", "signature_label", "compartment"])["signature_score"].median().unstack("compartment").reset_index()
        for col in ["malignant", "immune", "stromal"]:
            if col not in summary: summary[col] = np.nan
        summary = summary.rename(columns={"malignant": "median_malignant", "immune": "median_immune", "stromal": "median_stromal"})
    schema_df.to_csv(out / "scrna_schema_audit_v0_7.csv", index=False)
    count_df.to_csv(out / "scrna_group_counts_v0_7.csv", index=False)
    score_df.to_csv(out / "scrna_signature_scores_v0_7.csv", index=False)
    pair_df.to_csv(out / "scrna_pairwise_contrasts_v0_7.csv", index=False)
    summary.to_csv(out / "scrna_dataset_cancer_signature_summary_v0_7.csv", index=False)


def write_spatial(
    root: Path,
    queue_path: Path,
    input_root: Path,
    genes: dict[str, list[str]],
    min_genes: int,
) -> None:
    out = root / "results/spatial"; out.mkdir(parents=True, exist_ok=True)
    queue = pd.read_csv(queue_path, sep="\t", dtype=str)
    queue["remote_path"] = queue["remote_path"].map(
        lambda value: str((input_root / value).resolve()) if not Path(value).is_absolute() else str(Path(value))
    )
    counts, scores, summaries = [], [], []
    for row in queue.itertuples(index=False):
        path = Path(row.remote_path)
        if not path.is_file() or (str(row.expected_size).isdigit() and path.stat().st_size != int(row.expected_size)):
            continue
        c, s, m = process_spatial(pd.Series(row._asdict()), genes, min_genes)
        counts.append(c); scores.append(s); summaries.append(m)
    count_df = pd.concat(counts, ignore_index=True)
    score_df = pd.concat(scores, ignore_index=True)
    summary_df = pd.concat(summaries, ignore_index=True)
    pair_rows = []
    for label, group in summary_df.groupby("signature_label"):
        for left, right in [("malignant", "immune"), ("malignant", "stromal"), ("immune", "stromal")]:
            delta = (group[f"median_{left}"] - group[f"median_{right}"]).dropna()
            if len(delta): pair_rows.append({"signature_label": label, "contrast": f"{left}_minus_{right}", "paired_sample_n": len(delta), "median_delta": delta.median(), "mean_delta": delta.mean(), "positive_fraction": (delta > 0).mean()})
    count_df.to_csv(out / "spatial_brca_spot_compartment_counts_v0_7.csv", index=False)
    score_df.to_csv(out / "spatial_brca_signature_scores_v0_7.csv", index=False)
    summary_df.to_csv(out / "spatial_brca_sample_signature_summary_v0_7.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(out / "spatial_brca_pairwise_contrasts_v0_7.csv", index=False)
    plot = summary_df.pivot(index="signature_label", columns="sample_id", values="median_immune")
    fig, ax = plt.subplots(figsize=(7, 4.5)); im = ax.imshow(plot.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(plot.columns)), plot.columns, rotation=45, ha="right"); ax.set_yticks(range(len(plot.index)), plot.index); fig.colorbar(im, ax=ax, label="Immune-compartment median score"); fig.tight_layout()
    fig_dir = root / "figures/spatial"; fig_dir.mkdir(parents=True, exist_ok=True); fig.savefig(fig_dir / "fig_spatial_brca_sample_medians_v2.png", dpi=300); plt.close(fig)


def main() -> int:
    args = parse_args()
    root = Path(args.analysis_root).resolve()
    input_root = Path(args.input_root).expanduser().resolve()
    genes = load_signature_genes(Path(args.signature_table))
    if args.mode in {"scrna", "all"}:
        cache_dir = Path(args.source_cache_dir).resolve() if args.source_cache_dir else None
        write_scrna(root, Path(args.scrna_queue), input_root, genes, args.chunk_rows, args.min_group_cells, args.min_signature_genes, cache_dir)
    if args.mode in {"spatial", "all"}:
        write_spatial(root, Path(args.spatial_queue), input_root, genes, args.min_signature_genes)
    summary = {"mode": args.mode, "status": "complete", "reader": "h5py_required_fields_plus_bounded_csr_scan"}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
