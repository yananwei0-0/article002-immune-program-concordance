#!/usr/bin/env python3
"""Build the BIO-05 CURRENT_P2 documentation-only closure from retained evidence."""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import ast
import csv
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


PROJECT = "BIO-05"
SOURCE_P1 = "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md"
SOURCE_P1_SHA256 = "7ab1be4929eaa15030fc226b90172c5cbecbe69be42257470562112b83b5d5e2"
P2_PATH = "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P2.md"
FRESH_AUDIT = "current_p2_closure/FRESH_AUDIT_RESPONSE.md"
CONFIG_COPY = "current_p2_closure/PROGRAM_DEFINITION_CURRENT_P2.csv"
ENVIRONMENT = "current_p2_closure/SOFTWARE_ENVIRONMENT_CURRENT_P2.json"
RULES = "current_p2_closure/SUPPORT_MAPPING_RULES_CURRENT_P2.json"
TCGA_SELECTION = "current_p2_closure/TCGA_PRIMARY_SAMPLE_SELECTION_CURRENT_P2.csv"
DERIVATIONS = "current_p2_closure/NUMERIC_DERIVATIONS_CURRENT_P2.json"
NUMERIC_MAP = "current_p2_closure/NUMERIC_MAP_CURRENT_P2.csv"
TCGA_ANCHOR_SELECTION = "current_p2_closure/TCGA_ANCHOR_SELECTION_CURRENT_P2.csv"
CROSS_FREEZE_REVISION_JSON = "cross_freeze_revision/CROSS_FREEZE_REVISION_1.json"
CROSS_FREEZE_REVISION_MD = "cross_freeze_revision/CROSS_FREEZE_REVISION_1.md"
CROSS_FREEZE_AUDIT = Path(
    expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/sci27_current_p2_cross_freeze_20260815/responses/"
    "BIO-05_CURRENT_P2_CROSS_FREEZE_RESPONSE.md")
)
CROSS_FREEZE_AUDIT_SHA256 = "11919a32ea9eb54d90f55b8a556881c0458b45aa00eaf29ca7b681c13951a88f"
PRE_REVISION_P2_SHA256 = "8419628ed8382bbad34ee84322bc66412c77dd3197572d3e822baaa9820fc502"
PRE_REVISION_NUMERIC_MAP_SHA256 = "7fe0ac0f043c291d88685c7376f5f53b8bfd8c714e8f1a0a36d111a0ee74182f"

SNAPSHOTS = {
    "tcga": {
        "basename": "run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py",
        "expected_sha256": "219d3e3ca81377b897172830ce4741fb207beb9d97700b29d4de2858abcaf584",
    },
    "tcga_score_helper": {
        "basename": "run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py",
        "expected_sha256": "79d050d282f361af325b03c248ed8c9044fda6f9dbea439a24c02c0a20944abf",
    },
    "tcga_anchor_helper": {
        "basename": "run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap.py",
        "expected_sha256": "adeee39f677fe0e065d8be399f1b9e9702ce4cbbe7ab7ff0c60c13f53c7b3862",
    },
    "apollo": {
        "basename": "run_bio05_apollo_external_cohort_v0_7.py",
        "expected_sha256": "80c9178104d3350b249b33e02418aa50e40e76e9268bf1a88f4891e057cbf4c5",
    },
    "threeca": {
        "basename": "run_bio05_3ca_source_mapping_v0_8.py",
        "expected_sha256": "209c8d1599f9eb41167f3603bbb68473ec0852f9035b1d215c85914b6bec3939",
    },
}

# Numeric tokens exclude digits embedded in gene symbols and identifiers, while
# retaining dotted software versions, ordinary values, and scientific notation.
NUMERIC_RE = re.compile(
    r"(?<![A-Za-z0-9_])(?:[+-]?\d+(?:\.\d+){2,}|[+-]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][+-]?\d+)?)(?![A-Za-z0-9_])"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", default=".")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(x) for x in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def read_manifest_paths(root: Path) -> list[dict[str, Any]]:
    manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    return list(manifest.get("protocol_source_inputs", [])) + list(manifest.get("executable_assets", []))


def find_manifest_source(root: Path, basename: str, expected_sha256: str) -> Path:
    hits = []
    for row in read_manifest_paths(root):
        path = Path(str(row.get("path", "")))
        if path.name == basename and str(row.get("sha256", "")) == expected_sha256 and path.is_file():
            hits.append(path)
    if not hits:
        raise FileNotFoundError(f"No manifest source matched {basename} / {expected_sha256}")
    # Prefer the data-lake execution copy when duplicate documentary copies exist.
    return sorted(hits, key=lambda x: ("bioinformatics_analysis" not in str(x), str(x)))[0]


def copy_evidence_inputs(root: Path) -> dict[str, dict[str, Any]]:
    out_dir = root / "current_p2_closure" / "executed_code"
    out_dir.mkdir(parents=True, exist_ok=True)
    copied: dict[str, dict[str, Any]] = {}
    for key, spec in SNAPSHOTS.items():
        source = find_manifest_source(root, spec["basename"], spec["expected_sha256"])
        target = out_dir / spec["basename"]
        shutil.copyfile(source, target)
        observed = sha256_file(target)
        if observed != spec["expected_sha256"]:
            raise RuntimeError(f"Copied {key} snapshot hash mismatch: {observed}")
        copied[key] = {
            "relative_path": str(target.relative_to(root)),
            "sha256": observed,
            "size_bytes": target.stat().st_size,
        }

    # The executed 3CA script imported this data-lake helper. It was not separately
    # itemized in SOURCE_MANIFEST, so its identity is retained here without claiming
    # that the original source manifest covered it.
    three_source = find_manifest_source(
        root,
        SNAPSHOTS["threeca"]["basename"],
        SNAPSHOTS["threeca"]["expected_sha256"],
    )
    base_source = three_source.parent / "run_bio05_scrna_source_mapping_v0_7.py"
    if not base_source.is_file():
        raise FileNotFoundError(base_source)
    base_target = out_dir / base_source.name
    shutil.copyfile(base_source, base_target)
    copied["threeca_imported_scrna_helper"] = {
        "relative_path": str(base_target.relative_to(root)),
        "sha256": sha256_file(base_target),
        "size_bytes": base_target.stat().st_size,
        "source_manifest_itemized": False,
    }

    signature = find_manifest_source(
        root,
        "signature_gene_sets_v0_1.csv",
        "b631094186b4f1a69f3e2ee8b948009f10f4b7ad0b644cd1cf446e4152c46ef1",
    )
    config_target = root / CONFIG_COPY
    shutil.copyfile(signature, config_target)
    if sha256_file(config_target) != "b631094186b4f1a69f3e2ee8b948009f10f4b7ad0b644cd1cf446e4152c46ef1":
        raise RuntimeError("Program-definition copy is not byte-identical to the executed configuration")
    return copied


def capture_environment(root: Path) -> dict[str, Any]:
    run_manifest = json.loads((root / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
    runtime = run_manifest["runtime"]
    interpreter = Path(runtime["python_executable"])
    if not interpreter.is_file():
        raise FileNotFoundError(interpreter)
    probe = r'''
import json, platform, sys
import h5py, matplotlib, numpy, pandas, scipy
payload = {
  "python": platform.python_version(),
  "numpy": numpy.__version__,
  "pandas": pandas.__version__,
  "scipy": scipy.__version__,
  "h5py": h5py.__version__,
  "matplotlib": matplotlib.__version__,
}
print(json.dumps(payload, sort_keys=True))
'''
    observed = json.loads(subprocess.check_output([str(interpreter), "-c", probe], text=True))
    if observed["python"] != runtime["python_version"]:
        raise RuntimeError("Recorded and observed runtime Python versions differ")
    site_packages = interpreter.parent.parent / "lib" / f"python{observed['python'].rsplit('.', 1)[0]}" / "site-packages"
    metadata = {}
    for package, version in observed.items():
        if package == "python":
            continue
        candidates = sorted(site_packages.glob(f"{package.replace('-', '_')}-{version}.dist-info/METADATA"))
        if candidates:
            item = candidates[0]
            metadata[package] = {
                "metadata_sha256": sha256_file(item),
                "metadata_mtime_epoch": int(item.stat().st_mtime),
            }
    payload = {
        "project": PROJECT,
        "execution_window_utc": {
            "began": run_manifest["created_at_utc"],
            "completed": run_manifest["completed_at_utc"],
        },
        "recorded_runtime": {
            "python": runtime["python_version"],
            "platform": runtime["platform"],
        },
        "observed_versions_from_recorded_interpreter": observed,
        "package_metadata_evidence": metadata,
        "adjudication_basis": (
            "The interpreter is the executable recorded in RUN_MANIFEST.json. Package metadata files "
            "were installed before the retained execution window; their hashes are frozen here."
        ),
        "functions": {
            "spearman": "scipy.stats.spearmanr",
            "wilcoxon": "scipy.stats.wilcoxon",
            "kruskal": "scipy.stats.kruskal",
            "rank": "scipy.stats.rankdata(method='average')",
            "least_squares": "numpy.linalg.lstsq(rcond=None)",
            "rng": "numpy.random.default_rng",
        },
    }
    write_json(root / ENVIRONMENT, payload)
    return payload


def export_tcga_selection(root: Path) -> pd.DataFrame:
    scores = pd.read_csv(root / "results/tcga/tcga_xena_signature_scores_v0_5.csv")
    cols = [
        "sample_id",
        "patient_id",
        "sample_type_code",
        "sample_type_label",
        "cancer_hint",
        "primary_disease_hint",
    ]
    selected = scores[cols].drop_duplicates().sort_values(["patient_id", "sample_id"]).reset_index(drop=True)
    if len(selected) != 4780 or selected["patient_id"].nunique() != 4780:
        raise RuntimeError("TCGA selected-sample export did not reproduce 4780 one-per-patient rows")
    selected.to_csv(root / TCGA_SELECTION, index=False)
    return selected


def literal_assignment(path: Path, name: str) -> Any:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError(f"Literal assignment {name} not found in {path}")


def export_tcga_anchor_selection(
    root: Path,
    snapshots: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    helper_path = root / snapshots["tcga_anchor_helper"]["relative_path"]
    anchors = literal_assignment(helper_path, "SIGNATURE_ANCHORS")
    labels = literal_assignment(helper_path, "SIGNATURE_LABELS")
    summary_path = root / "results/tcga/tcga_xena_signature_summary_v0_5.csv"
    meta_path = root / "results/tcga/tcga_xena_68immune_meta_summary_v0_5.csv"
    summary = pd.read_csv(summary_path)
    meta = pd.read_csv(meta_path)
    rows: list[dict[str, Any]] = []
    for signature, anchor_names in anchors.items():
        retained = row_one(summary, signature=signature)
        selected = meta[
            meta["signature"].eq(signature) & meta["immune_signature_68"].isin(anchor_names)
        ].copy()
        retained_names = {
            value for value in str(retained["selected_anchor_names"]).split(";") if value
        }
        expected_names = set(anchor_names)
        median_rho = float(selected["median_within_cancer_rho"].median())
        min_positive = float(selected["positive_cancer_fraction"].min())
        membership_match = (
            len(selected) == len(anchor_names)
            and retained_names == expected_names
            and int(retained["selected_anchor_count"]) == len(anchor_names)
        )
        summary_match = (
            np.isclose(median_rho, float(retained["selected_anchor_median_rho"]), rtol=1e-12, atol=1e-14)
            and np.isclose(
                min_positive,
                float(retained["selected_anchor_min_positive_cancer_fraction"]),
                rtol=1e-12,
                atol=1e-14,
            )
        )
        if not membership_match or not summary_match:
            raise RuntimeError(f"TCGA anchor recovery mismatch for {signature}")
        selected_by_name = selected.set_index("immune_signature_68")
        for order, anchor_name in enumerate(anchor_names, start=1):
            anchor_row = selected_by_name.loc[anchor_name]
            rows.append(
                {
                    "signature": signature,
                    "signature_label": labels[signature],
                    "anchor_order": order,
                    "immune_signature_68": anchor_name,
                    "cancer_strata": int(anchor_row["cancer_strata"]),
                    "median_within_cancer_rho": float(anchor_row["median_within_cancer_rho"]),
                    "positive_cancer_fraction": float(anchor_row["positive_cancer_fraction"]),
                    "retained_selected_anchor_count": int(retained["selected_anchor_count"]),
                    "retained_selected_anchor_median_rho": float(retained["selected_anchor_median_rho"]),
                    "retained_selected_anchor_min_positive_cancer_fraction": float(
                        retained["selected_anchor_min_positive_cancer_fraction"]
                    ),
                    "membership_match": membership_match,
                    "summary_match": summary_match,
                }
            )
    artifact = pd.DataFrame(rows)
    if len(artifact) != 17 or artifact["signature"].nunique() != 8:
        raise RuntimeError("Recovered TCGA anchor table was not the expected 17 anchors across 8 programs")
    artifact.to_csv(root / TCGA_ANCHOR_SELECTION, index=False)
    return {
        "source_snapshot": snapshots["tcga_anchor_helper"],
        "artifact": TCGA_ANCHOR_SELECTION,
        "anchor_count": len(artifact),
        "program_count": int(artifact["signature"].nunique()),
        "anchors": anchors,
        "selection": (
            "For each program, retain 68-immune meta-summary rows whose immune_signature_68 is in the fixed "
            "program-specific anchor list; selected_anchor_median_rho is the median of their "
            "median_within_cancer_rho values, and selected_anchor_min_positive_cancer_fraction is the minimum "
            "of their positive_cancer_fraction values."
        ),
        "validation": (
            "All 17 configured anchors were present; membership, selected-anchor counts, median aggregation, "
            "and minimum-positive-fraction aggregation matched the retained signature summary for all 8 programs."
        ),
    }


def write_rules(
    root: Path,
    snapshots: dict[str, dict[str, Any]],
    tcga_anchor_rule: dict[str, Any],
) -> dict[str, Any]:
    rules = {
        "project": PROJECT,
        "program_configuration": {
            "artifact": CONFIG_COPY,
            "sha256": sha256_file(root / CONFIG_COPY),
            "version": "v0_1",
            "cptac_identifier_canonicalization": {
                "gene_symbol": "strip whitespace and convert to uppercase",
                "feature_identifier": "extract ENSG followed by digits and remove an optional dot-version suffix",
                "mapping": "for each program gene, union all matrix rows matching any listed stable ENSG; average standardized matching rows within sample",
                "unmapped_definition_rows": "retained in requested counts but do not enter CPTAC scores",
            },
            "cptac_standardization": {
                "axis": "feature across eligible tumor samples within cancer and matrix",
                "standard_deviation": "sample standard deviation (ddof=1)",
                "zero_or_nonfinite_standard_deviation": "replace standard deviation with missing; resulting feature z-scores are missing",
                "gene_collapse": "mean of all finite standardized feature rows mapped to a gene within sample",
                "program_score": "equal-weight mean of finite gene-level values; require at least 3 nonmissing recovered genes",
                "imputation": "none",
            },
        },
        "tcga": {
            "code_snapshot": snapshots["tcga"],
            "score_helper_snapshot": snapshots["tcga_score_helper"],
            "selection_table": TCGA_SELECTION,
            "selection": {
                "layer": "rna_expression",
                "sample_type_code": 1,
                "eligible_cancers": ["BRCA", "KIRC", "COAD", "GBM", "HNSC", "LUSC", "LUAD", "OV", "PAAD", "UCEC"],
                "tie_break": "sort ascending by patient_id then sample_id; keep the first sample per patient",
            },
            "gene_mapping": "union of program gene-symbol row IDs and DepMap-derived Entrez row IDs that are present in the TCGA matrix",
            "standardization": "within cancer, selected expression rows standardized across samples using sample standard deviation (ddof=1); zero/nonfinite standard deviations become missing",
            "score_aggregation_and_missingness": {
                "matrix_level_recovery_gate": "emit a program only when at least 3 mapped gene symbols are represented by present row IDs",
                "program_score": "unweighted arithmetic mean across finite standardized selected expression rows for each patient (pandas mean with skipna=True)",
                "sample_level_finite_gene_gate": "none beyond at least one finite standardized selected row; an all-missing row set yields a missing score",
                "imputation": "none after the stated row mapping and standardization",
            },
            "anchor_selection": tcga_anchor_rule,
        },
        "apollo": {
            "code_snapshot": snapshots["apollo"],
            "operation_order": {
                "LUAD_RNA": [
                    "outer-union per-sample count-UQ tables by gene",
                    "fill absent gene/sample combinations with zero",
                    "apply log2(count+1)",
                    "uppercase gene symbols and average duplicate rows",
                    "z-standardize each gene across samples with ddof=0; zero/nonfinite standard deviations become missing",
                    "average recovered genes with equal weight",
                ],
                "LUAD_proteomics": [
                    "canonicalize sample identifiers",
                    "uppercase gene symbols and average duplicate rows",
                    "z-standardize each gene across samples with ddof=0; zero/nonfinite standard deviations become missing",
                    "average recovered genes with equal weight",
                ],
                "LUAD_phosphoproteomics": [
                    "retain relative columns and canonicalize sample identifiers",
                    "uppercase Primary_Gene values and average phosphosite rows within gene",
                    "z-standardize each collapsed gene across samples with ddof=0; zero/nonfinite standard deviations become missing",
                    "average recovered genes with equal weight",
                ],
                "OV_RNA_and_proteomics": [
                    "uppercase gene symbols and average duplicate rows",
                    "z-standardize each gene across samples with ddof=0; zero/nonfinite standard deviations become missing",
                    "average recovered genes with equal weight",
                ],
            },
            "score_aggregation_and_missingness": {
                "matrix_level_recovery_gate": "emit a program only when at least 3 program genes are present in the collapsed layer matrix",
                "program_score": "unweighted arithmetic mean of finite standardized recovered-gene values within each sample (pandas mean with skipna=True)",
                "sample_level_finite_gene_gate": "no additional three-finite-gene rule; at least one finite standardized recovered-gene value is sufficient",
                "all_missing_rule": "if all recovered-gene standardized values are missing for a sample, its program score is missing and is removed by complete-case joining",
                "imputation": "none after the explicitly documented LUAD RNA outer-union zero fill",
            },
        },
        "single_cell_cellxgene": {
            "code": {"relative_path": "code/run_h5ad_support_v2.py", "sha256": sha256_file(root / "code/run_h5ad_support_v2.py")},
            "field_selection": {
                "donor": ["donor_id", "patient_id", "participant_id", "regex donor|patient|participant fallback"],
                "cell_type": ["cell_type", "author_cell_type", "regex cell.?type|cluster|annotation fallback"],
                "disease": ["disease", "regex disease|cancer|tumour fallback"],
                "tissue": ["tissue", "regex tissue|organ fallback"],
                "actual_fields": "results/scrna/scrna_schema_audit_v0_7.csv",
            },
            "priority_order": [
                "normal_cell_call=cancer -> malignant",
                "normal_cell_call=normal and epithelial -> other",
                "celltype_major exact cancer/immune/stromal/normal rules",
                "dataset-specific other then malignant exact-value overrides",
                "malignant then immune then stromal keyword rules",
                "dataset-specific epithelial malignant fallback",
                "otherwise other",
            ],
            "dataset_specific_malignant": {
                "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": ["squamous epithelial cell"],
                "0bebef1a-4607-4584-9070-dacf89a0d635": ["epithelial cell of lung"],
            },
            "dataset_specific_other": {
                "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": ["salivary gland glandular cell"]
            },
            "epithelial_fallback_source_literals_in_executed_code": [
                "dea97145-f712-431c-a223-6b5f565f362a",
                "4d82bd8e-8827-410b-96ae-2c806799c08c",
            ],
            "malignant_keywords": ["abnormal cell", "cancer epithelial", "malignant", "neoplastic", "tumor", "tumour"],
            "immune_keywords": [
                "b cell", "dendritic", "granulocyte", "immune cell", "lymphocyte", "macrophage", "mast cell",
                "megakaryocyte", "microgl", "monocyte", "myeloid", "natural killer", "neutrophil", "nk", "plasma cell", "t cell",
            ],
            "stromal_keywords": [
                "astrocyte", "endothelial", "fibroblast", "mesangial", "mural cell", "myofibroblast", "neuron",
                "oligodendrocyte", "oligodendrocyte precursor", "pericyte", "podocyte", "radial glial", "schwann",
                "smooth muscle", "stromal", "support cell", "vascular associated smooth muscle cell",
            ],
            "score_aggregation_and_missingness": {
                "group_gate": "at least 20 cells per source-specific donor-cancer-compartment group",
                "matrix_level_recovery_gate": "emit a program only when at least 3 program genes are present after duplicate-symbol collapse",
                "program_score": "unweighted arithmetic mean of finite within-cancer standardized recovered-gene pseudobulk values (pandas mean with skipna=True)",
                "sample_level_finite_gene_gate": "no additional three-finite-gene rule; an all-missing recovered-gene vector yields a missing donor-compartment score",
                "imputation": "none",
            },
        },
        "single_cell_3ca": {
            "code_snapshot": snapshots["threeca"],
            "imported_helper_snapshot": snapshots["threeca_imported_scrna_helper"],
            "actual_fields": "results/threeca/threeca_schema_audit_v0_8.csv",
            "pdac": {
                "malignant": ["malignant"],
                "immune": ["macrophage", "t_cell", "b_cell", "mast"],
                "stromal": ["fibroblast", "stellate", "endothelial"],
            },
            "ucec": {
                "malignant": ["malignant"],
                "immune": ["macrophage", "t_cell", "b_cell", "mast"],
                "stromal": ["fibroblast", "endothelial", "myocyte"],
            },
            "other_rule": "all unlisted normalized cell_type values -> other",
            "score_aggregation_and_missingness": {
                "group_gate": "at least 20 cells per source-specific donor-cancer-compartment group",
                "matrix_level_recovery_gate": "emit a program only when at least 3 program genes are present after duplicate-symbol collapse",
                "program_score": "unweighted arithmetic mean of finite within-cancer standardized recovered-gene pseudobulk values using the imported v0.7 helper (pandas mean with skipna=True)",
                "sample_level_finite_gene_gate": "no additional three-finite-gene rule; an all-missing recovered-gene vector yields a missing donor-compartment score",
                "imputation": "none",
            },
        },
        "spatial": {
            "code": {"relative_path": "code/run_h5ad_support_v2.py", "sha256": sha256_file(root / "code/run_h5ad_support_v2.py")},
            "score_columns": {
                "malignant": ["Cancer Basal SC", "Cancer Cycling", "Cancer Her2 SC", "Cancer LumA SC", "Cancer LumB SC"],
                "immune": [
                    "B cells Memory", "B cells Naive", "Cycling T-cells", "Cycling_Myeloid", "DCs", "Macrophage",
                    "Monocyte", "NK cells", "NKT cells", "Plasmablasts", "T cells CD4+", "T cells CD8+",
                ],
                "stromal": [
                    "CAFs MSC iCAF-like", "CAFs myCAF-like", "Endothelial ACKR1", "Endothelial CXCL12",
                    "Endothelial Lymphatic LYVE1", "Endothelial RGS5", "Myoepithelial", "PVL Differentiated", "PVL Immature",
                ],
            },
            "assignment": [
                "sum present columns within each of the three score groups, coercing nonnumeric/missing entries to zero",
                "if the largest group score is positive, assign its argmax; ties follow malignant, immune, stromal insertion order",
                "otherwise Classification=stroma -> stromal",
                "otherwise Classification containing lymphocytes -> immune",
                "otherwise Classification containing invasive cancer -> malignant",
                "otherwise other",
            ],
            "score_aggregation_and_missingness": {
                "matrix_level_recovery_gate": "emit a program only when at least 3 program genes are present in the spatial matrix",
                "program_score": "per-spot unweighted numpy.nanmean across recovered genes after within-sample ddof=0 gene standardization",
                "spot_level_finite_gene_gate": "no additional three-finite-gene rule; no finite recovered-gene z value yields a missing spot score",
                "sample_compartment_summary": "pandas median of finite spot scores within sample and compartment; missing spot scores are skipped",
                "imputation": "none",
            },
        },
    }
    write_json(root / RULES, rules)
    return rules


def add_fact(
    facts: dict[str, dict[str, Any]],
    fact_id: str,
    value: Any,
    source_file: str,
    locator: str,
    field: str,
    derivation: str,
) -> None:
    facts[fact_id] = {
        "value": json_safe(value),
        "source_file": source_file,
        "source_locator": locator,
        "source_field": field,
        "derivation": derivation,
    }


def row_one(frame: pd.DataFrame, **filters: Any) -> pd.Series:
    mask = pd.Series(True, index=frame.index)
    for field, value in filters.items():
        mask &= frame[field].eq(value)
    hit = frame[mask]
    if len(hit) != 1:
        raise RuntimeError(f"Expected one row for {filters}; found {len(hit)}")
    return hit.iloc[0]


def compute_facts(root: Path) -> dict[str, dict[str, Any]]:
    facts: dict[str, dict[str, Any]] = {}
    primary_path = "results/cptac/cptac_primary_concordance_registry_v2.csv"
    hetero_path = "results/cptac/cptac_heterogeneity_summary_v2.csv"
    feature_path = "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv"
    threshold_path = "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv"
    loco_path = "results/cptac/cptac_leave_one_cancer_out_v2.csv"
    pooled_path = "results/cptac/cptac_within_cancer_standardized_pooled_sensitivity_v2.csv"
    residual_path = "results/cptac/cptac_residual_context_registry_v2.csv"
    apollo_path = "results/apollo/apollo_concordance_registry_v2.csv"
    scrna_path = "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv"
    spatial_path = "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv"
    tcga_corr_path = "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv"
    tcga_summary_path = "results/tcga/tcga_xena_signature_summary_v0_5.csv"
    primary = pd.read_csv(root / primary_path)
    hetero = pd.read_csv(root / hetero_path)
    feature = pd.read_csv(root / feature_path)
    threshold = pd.read_csv(root / threshold_path)
    loco = pd.read_csv(root / loco_path)
    pooled = pd.read_csv(root / pooled_path)
    residual = pd.read_csv(root / residual_path)
    apollo = pd.read_csv(root / apollo_path)
    scrna = pd.read_csv(root / scrna_path)
    spatial = pd.read_csv(root / spatial_path)
    tcga_corr = pd.read_csv(root / tcga_corr_path)
    tcga_summary = pd.read_csv(root / tcga_summary_path)
    support = json.loads((root / "results/support_layers_summary_v2.json").read_text(encoding="utf-8"))
    env = json.loads((root / ENVIRONMENT).read_text(encoding="utf-8"))
    config = pd.read_csv(root / CONFIG_COPY)

    def f(fid: str, value: Any, source: str, locator: str, field: str, derivation: str) -> None:
        add_fact(facts, fid, value, source, locator, field, derivation)

    f("cancer_count", primary["cancer"].nunique(), primary_path, "all rows", "cancer", "nunique")
    f("program_count", config["signature"].nunique(), CONFIG_COPY, "all rows", "signature", "nunique")
    f("layer_count", primary["comparison_layer"].nunique(), primary_path, "all rows", "comparison_layer", "nunique")
    f("primary_registry_rows", len(primary), primary_path, "all rows", "row count", "count")
    for fid, value, locator in [
        ("min_score_genes", 3, "MIN_GENES = 3; score gate at lines 293-300"),
        ("cptac_ddof", 1, "zscore_rows at lines 207-212"),
        ("primary_min_n", 20, "MIN_N = 20; evaluate_pair at lines 467-496"),
        ("primary_bootstrap_resamples", 20000, "PRIMARY_BOOTSTRAP_RESAMPLES = 20_000"),
        ("primary_permutation_resamples", 999999, "PRIMARY_PERMUTATION_RESAMPLES = 999_999"),
        ("primary_asymptotic_n_cutoff", 30, "PRIMARY_PERMUTATION_N_MAX = 30"),
        ("seed_digest_bytes", 8, "stable_seed uses sha256 digest()[:8]"),
        ("ci_percent", 95, "BCa quantiles 0.025 and 0.975"),
    ]:
        f(fid, value, "code/run_cptac_v2.py", locator, "code literal", "executed code constant")
    f("primary_min_distinct_values", 3, "code/run_cptac_v2.py", "evaluate_pair lines 468-495", "x_distinct/y_distinct gate", "executed operation-specific code literal")
    f("residual_min_distinct_values", 3, "code/run_cptac_v2.py", "residual model gate lines 699-709", "rna_score/omic_score nunique gate", "executed operation-specific code literal")
    f("primary_q_cutoff", 0.05, "code/run_cptac_v2.py", "results summary line 1076", "q_value_bh threshold", "executed primary result-classification literal")
    f("bca_alpha_low", 0.025, "code/run_cptac_v2.py", "spearman_bca_ci line 406", "alpha", "code literal")
    f("bca_alpha_high", 0.975, "code/run_cptac_v2.py", "spearman_bca_ci line 406", "alpha", "code literal")
    f("i2_lower_bound", 0, "code/run_cptac_v2.py", "build_heterogeneity line 607", "max lower bound", "code literal")
    f("tcga_primary_type_code", 1, "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", "load_tcga_meta lines 127-134", "sample_type_code", "code literal")
    f("tcga_min_correlation_n", 25, "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", "MIN_CORRELATION_N", "code literal", "executed code constant")
    f("tcga_min_subtype_n", 20, "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", "MIN_SUBTYPE_N line 65; subtype gate lines 327-350", "subtype group n gate", "executed operation-specific code literal")
    f("tcga_min_program_genes", 3, "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", "compute_signature_scores lines 181-209", "mapped_gene_symbol_count gate", "executed operation-specific code literal")
    f("tcga_ddof", 1, "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py", "zscore_rows lines 404-408", "ddof", "executed imported helper operation")
    f("tcga_q_cutoff", 0.05, "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", "meta-summary line 304", "fdr_bh threshold", "executed TCGA support-summary literal")
    f("tcga_anchor_rule_rows", 17, TCGA_ANCHOR_SELECTION, "all rows", "row count", "count of fixed program-anchor links recovered from executed helper")
    f("apollo_ddof", 0, "current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py", "zscore_columns lines 85-94", "ddof", "code literal")
    f("apollo_log_pseudocount", 1, "current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py", "load_apollo_luad_rna line 189", "pseudocount", "code literal")
    f("apollo_min_program_genes", 3, "current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py", "MIN_SIGNATURE_GENES line 24; score_signature_matrix lines 113-125", "matrix-level recovered-gene gate", "executed operation-specific code literal")
    f("support_ddof", 0, "code/run_h5ad_support_v2.py", "zscore_columns lines 357-360 and spatial lines 535-539", "ddof", "code literal")
    f("apollo_min_n", 10, "code/postprocess_support_v2.py", "apollo_registry lines 119-128", "matched n gate", "code literal")
    f("apollo_min_distinct_values", 3, "code/postprocess_support_v2.py", "apollo_registry line 121", "rna_score/omic_score nunique gate", "executed operation-specific code literal")
    f("apollo_family_rows", 8, apollo_path, "group by fdr_family", "row count", "each of 3 families has 8 rows")
    f("apollo_family_count", apollo["fdr_family"].nunique(), apollo_path, "all rows", "fdr_family", "nunique")
    f("apollo_registered_rows", len(apollo), apollo_path, "all rows", "row count", "count")
    f("apollo_fisher_z_constant", 1.9599639845, "code/postprocess_support_v2.py", "spearman_ci lines 82-87", "normal critical value", "code literal")
    f("apollo_ci_percent", 95, "code/postprocess_support_v2.py", "spearman_ci lines 82-87", "two-sided nominal confidence level", "coverage implied by the executed 1.9599639845 normal critical value")
    f("fisher_z_n_offset", 3, "code/postprocess_support_v2.py", "spearman_ci line 86", "n offset", "code literal")
    f("wilcoxon_auto_exact_max_n", 50, ENVIRONMENT, "SciPy 1.18.0 scipy.stats.wilcoxon method='auto'", "len(d)", "exact when no ties/zeros and len(d)<=50")
    f("wilcoxon_auto_tied_permutation_max_n", 13, ENVIRONMENT, "SciPy 1.18.0 scipy.stats.wilcoxon method='auto'", "len(d)", "exhaustive sign permutations for ties/zeros and len(d)<=13")
    f("wilcoxon_all_zero_p", 1, "code/postprocess_support_v2.py", "signed_rank lines 246-249", "P value", "manual all-zero special case")
    f("scrna_cellxgene_min_group_cells", 20, "code/run_h5ad_support_v2.py", "parse_args line 94; group gate lines 404-409", "min_group_cells", "executed streaming H5AD code literal")
    f("scrna_cellxgene_min_program_genes", 3, "code/run_h5ad_support_v2.py", "parse_args line 95; program gate lines 434-444", "min_signature_genes", "executed streaming H5AD code literal")
    f("scrna_3ca_min_group_cells", 20, "current_p2_closure/executed_code/run_bio05_3ca_source_mapping_v0_8.py", "MIN_GROUP_CELLS line 23; group gate lines 234-244", "min_group_cells", "executed 3CA code literal")
    f("scrna_3ca_min_program_genes", 3, "current_p2_closure/executed_code/run_bio05_3ca_source_mapping_v0_8.py", "MIN_SIGNATURE_GENES line 24; program gate lines 273-277", "min_signature_genes", "executed 3CA code literal")
    f("spatial_min_program_genes", 3, "code/run_h5ad_support_v2.py", "parse_args line 95; spatial program gate lines 529-539", "min_signature_genes", "executed spatial code literal")
    f("scrna_min_paired_n", 6, "code/postprocess_support_v2.py", "signed_rank default min_n", "minimum paired donors", "code literal")
    f("scrna_family_rows", 24, "code/postprocess_support_v2.py", "8 programs x 3 contrasts per source-cancer", "family size", "deterministic product")
    f("scrna_sensitivity_n", 10, "code/postprocess_support_v2.py", "n_ge_10_sensitivity_window", "paired donor threshold", "code literal")
    f("scrna_support_q_cutoff", 0.05, "code/postprocess_support_v2.py", "within-source support rules lines 287-295", "q_value_bh threshold", "executed source-bounded support literal")
    f("spatial_sample_n", spatial["paired_sample_n"].dropna().unique()[0], spatial_path, "all rows", "paired_sample_n", "single unique value")
    for package, version in env["observed_versions_from_recorded_interpreter"].items():
        f(f"software_{package}", version, ENVIRONMENT, "observed_versions_from_recorded_interpreter", package, "identity")

    evaluable = primary[primary["status"].eq("evaluable")]
    prot = evaluable[evaluable["comparison_layer"].eq("proteomics")]
    phospho = evaluable[evaluable["comparison_layer"].eq("phosphoproteomics")]
    f("primary_unique_cases", json.loads((root / "RESULTS_SUMMARY.json").read_text())["primary"]["unique_cases_any_matched_layer"], "RESULTS_SUMMARY.json", ".primary.unique_cases_any_matched_layer", "value", "JSON field")
    f("primary_evaluable_rows", len(evaluable), primary_path, "status=evaluable", "row count", "count")
    f("primary_not_evaluable_rows", len(primary) - len(evaluable), primary_path, "status!=evaluable", "row count", "count")
    f("matched_min", pd.read_csv(root / "tables/TABLE1_CPTAC_CASE_FLOW_V2.csv")["maximum_matched_unique_cases"].min(), "tables/TABLE1_CPTAC_CASE_FLOW_V2.csv", "all cancer-layer rows", "maximum_matched_unique_cases", "minimum")
    f("matched_max", pd.read_csv(root / "tables/TABLE1_CPTAC_CASE_FLOW_V2.csv")["maximum_matched_unique_cases"].max(), "tables/TABLE1_CPTAC_CASE_FLOW_V2.csv", "all cancer-layer rows", "maximum_matched_unique_cases", "maximum")
    for prefix, frame in [("protein", prot), ("phosphoprotein", phospho)]:
        f(f"{prefix}_evaluable", len(frame), primary_path, f"comparison_layer={frame['comparison_layer'].iloc[0]}, status=evaluable", "row count", "count")
        f(f"{prefix}_q_lt_005", int((frame["q_value_bh"] < 0.05).sum()), primary_path, f"{prefix} evaluable q<0.05", "q_value_bh", "count")
        f(f"{prefix}_rho_min", frame["spearman_rho"].min(), primary_path, f"{prefix} evaluable", "spearman_rho", "minimum")
        f(f"{prefix}_rho_max", frame["spearman_rho"].max(), primary_path, f"{prefix} evaluable", "spearman_rho", "maximum")
        f(f"{prefix}_rho_median", frame["spearman_rho"].median(), primary_path, f"{prefix} evaluable", "spearman_rho", "median")
    nonsig = row_one(primary, cancer="pdac", program_label="Proteasome", comparison_layer="proteomics")
    for field, fid in [("complete_pair_n", "nonsig_n"), ("spearman_rho", "nonsig_rho"), ("p_value", "nonsig_p"), ("q_value_bh", "nonsig_q")]:
        f(fid, nonsig[field], primary_path, "cancer=pdac; program_label=Proteasome; comparison_layer=proteomics", field, "direct cell")
    f("positive_fraction_all", hetero["positive_cancer_fraction"].dropna().min(), hetero_path, "rows with evaluable cancers", "positive_cancer_fraction", "minimum")
    f("i2_min", hetero["descriptive_i2"].dropna().min(), hetero_path, "rows with evaluable cancers", "descriptive_i2", "minimum")
    f("i2_max", hetero["descriptive_i2"].dropna().max(), hetero_path, "rows with evaluable cancers", "descriptive_i2", "maximum")
    f("zero_evaluable_program_layers", 16 - len(hetero), "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv", "evaluable_cancer_n=0", "row count", "count")

    # Publication table cells are all individually addressable.
    table = pd.read_csv(root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv")
    table_fields = [
        "registered_cancer_rows", "evaluable_cancer_n", "total_complete_pairs_across_strata",
        "median_within_cancer_rho", "positive_cancer_fraction", "descriptive_i2",
    ]
    for ridx, row in table.iterrows():
        key = f"table_{ridx+1:02d}"
        locator = f"comparison_layer={row.comparison_layer}; program_label={row.program_label}"
        for field in table_fields:
            value = row[field]
            if pd.isna(value):
                continue
            f(f"{key}_{field}", value, "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv", locator, field, "direct cell")

    joined = primary[["cancer", "program", "comparison_layer", "spearman_rho"]].merge(
        feature[["cancer", "program", "comparison_layer", "spearman_rho"]],
        on=["cancer", "program", "comparison_layer"], suffixes=("_gene", "_feature"), validate="one_to_one",
    ).dropna()
    phospho_joined = joined[joined["comparison_layer"].eq("phosphoproteomics")]
    f("sensitivity_joined", len(joined), feature_path, "finite one-to-one join to primary", "row count", "count")
    f("sensitivity_phospho_joined", len(phospho_joined), feature_path, "phosphoproteomics finite join", "row count", "count")
    delta = (phospho_joined["spearman_rho_gene"] - phospho_joined["spearman_rho_feature"]).abs()
    f("sensitivity_phospho_median_delta", delta.median(), feature_path, "phosphoproteomics finite join", "absolute rho difference", "median")
    f("sensitivity_phospho_max_delta", delta.max(), feature_path, "phosphoproteomics finite join", "absolute rho difference", "maximum")
    f("threshold_evaluable", int(threshold["status"].eq("evaluable").sum()), threshold_path, "status=evaluable", "row count", "count")
    f("threshold_supported", int((threshold["q_value_bh"] < 0.05).sum()), threshold_path, "q_value_bh<0.05", "row count", "count")
    f("loco_rows", len(loco), loco_path, "all rows", "row count", "count")
    f("loco_median_rho_min", loco["leave_one_cancer_out_median_rho"].min(), loco_path, "all rows", "leave_one_cancer_out_median_rho", "minimum")
    f("loco_median_rho_max", loco["leave_one_cancer_out_median_rho"].max(), loco_path, "all rows", "leave_one_cancer_out_median_rho", "maximum")
    f("loco_positive_fraction_min", loco["leave_one_cancer_out_positive_fraction"].min(), loco_path, "all rows", "leave_one_cancer_out_positive_fraction", "minimum")
    pooled_finite = pooled[pooled["within_cancer_standardized_pooled_spearman_rho"].notna()].copy()
    f("pooled_registered_rows", len(pooled), pooled_path, "all rows", "row count", "count")
    f("pooled_finite_rows", len(pooled_finite), pooled_path, "within_cancer_standardized_pooled_spearman_rho finite", "row count", "count")
    f("pooled_positive_rows", int((pooled_finite["within_cancer_standardized_pooled_spearman_rho"] > 0).sum()), pooled_path, "finite rho > 0", "row count", "count")
    f("pooled_rho_min", pooled_finite["within_cancer_standardized_pooled_spearman_rho"].min(), pooled_path, "finite rows", "within_cancer_standardized_pooled_spearman_rho", "minimum")
    f("pooled_rho_max", pooled_finite["within_cancer_standardized_pooled_spearman_rho"].max(), pooled_path, "finite rows", "within_cancer_standardized_pooled_spearman_rho", "maximum")
    influence = primary.loc[primary["max_abs_leave_one_out_rho_delta"].idxmax()]
    for field, fid in [
        ("max_abs_leave_one_out_rho_delta", "max_influence_delta"), ("complete_pair_n", "max_influence_n"),
        ("spearman_rho", "max_influence_rho"), ("ci_lower", "max_influence_ci_low"),
        ("ci_upper", "max_influence_ci_high"), ("q_value_bh", "max_influence_q"),
    ]:
        f(fid, influence[field], primary_path, "row at maximum max_abs_leave_one_out_rho_delta", field, "direct cell")

    f("residual_strata", int(pd.read_csv(root / "results/cptac/cptac_residual_model_diagnostics_v2.csv")["status"].eq("evaluable").sum()), "results/cptac/cptac_residual_model_diagnostics_v2.csv", "status=evaluable", "row count", "count")
    f("residual_registry", len(residual), residual_path, "all rows", "row count", "count")
    f("residual_linear_finite", residual["linear_p_value"].notna().sum(), residual_path, "linear_p_value finite", "row count", "count")
    f("residual_rank_finite", residual["rank_p_value"].notna().sum(), residual_path, "rank_p_value finite", "row count", "count")
    f("residual_linear_supported", int((residual["linear_q_value_bh"] < 0.05).sum()), residual_path, "linear_q_value_bh<0.05", "row count", "count")
    f("residual_rank_supported", int((residual["rank_q_value_bh"] < 0.05).sum()), residual_path, "rank_q_value_bh<0.05", "row count", "count")
    both = (residual["linear_q_value_bh"] < 0.05) & (residual["rank_q_value_bh"] < 0.05)
    f("residual_both_supported", int(both.sum()), residual_path, "linear and rank q<0.05", "row count", "count")
    for variable, fid in [
        ("cibersort_myeloid", "residual_myeloid_supported"), ("cibersort_t_cell", "residual_tcell_supported"),
        ("estimate_immune_score", "residual_estimate_supported"), ("estimate_tumor_purity", "residual_purity_supported"),
    ]:
        f(fid, int(((residual["context_variable"] == variable) & (residual["linear_q_value_bh"] < 0.05)).sum()), residual_path, f"context_variable={variable}; linear_q_value_bh<0.05", "row count", "count")

    f("tcga_patients", support["tcga"]["primary_tumor_patients"], "results/support_layers_summary_v2.json", ".tcga.primary_tumor_patients", "value", "JSON field")
    f("tcga_samples", support["tcga"]["scored_samples"], "results/support_layers_summary_v2.json", ".tcga.scored_samples", "value", "JSON field")
    f("tcga_corr_rows", len(tcga_corr), tcga_corr_path, "all rows", "row count", "count")
    f("tcga_immune_signature_count", tcga_corr["immune_signature_68"].nunique(), tcga_corr_path, "all rows", "immune_signature_68", "nunique")
    f("tcga_supported", int((tcga_corr["fdr_bh"] < 0.05).sum()), tcga_corr_path, "fdr_bh<0.05", "row count", "count")
    nonprot = tcga_summary[~tcga_summary["signature_label"].eq("Proteasome")]
    f("tcga_anchor_min", nonprot["selected_anchor_median_rho"].min(), tcga_summary_path, "signature_label!=Proteasome", "selected_anchor_median_rho", "minimum")
    f("tcga_anchor_max", nonprot["selected_anchor_median_rho"].max(), tcga_summary_path, "signature_label!=Proteasome", "selected_anchor_median_rho", "maximum")
    f("tcga_anchor_positive", nonprot["selected_anchor_min_positive_cancer_fraction"].min(), tcga_summary_path, "signature_label!=Proteasome", "selected_anchor_min_positive_cancer_fraction", "minimum")
    f("tcga_proteasome", row_one(tcga_summary, signature_label="Proteasome")["selected_anchor_median_rho"], tcga_summary_path, "signature_label=Proteasome", "selected_anchor_median_rho", "direct cell")
    subtype_tests = pd.read_csv(root / "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv")
    f("tcga_subtype_tests", len(subtype_tests), "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv", "all rows", "row count", "count")
    subtype = pd.read_csv(root / "results/tcga/tcga_xena_immune_subtype_summary_v0_5.csv")
    f("tcga_c5_n", subtype[subtype["immune_subtype_short"].eq("Immune C5")]["n"].max(), "results/tcga/tcga_xena_immune_subtype_summary_v0_5.csv", "immune_subtype_short=Immune C5", "n", "maximum")

    for status, fid in [("evaluable", "apollo_evaluable"), ("supported", "apollo_supported"), ("opposite", "apollo_opposite"), ("null_result", "apollo_null"), ("not_evaluable", "apollo_not_evaluable")]:
        if status == "evaluable":
            value = int(apollo["status"].eq(status).sum())
            field = "status"
        else:
            value = int(apollo["evidence_status"].eq(status).sum())
            field = "evidence_status"
        f(fid, value, apollo_path, f"{field}={status}", "row count", "count")
    for cohort, layer, prefix in [
        ("APOLLO-LUAD", "proteomics", "apollo_luad_ifng"),
        ("APOLLO-OV", "proteomics", "apollo_ov_ifng"),
        ("APOLLO-LUAD", "phosphoproteomics", "apollo_luad_ifng_opposite"),
    ]:
        row = row_one(apollo, cohort_label=cohort, comparison_layer=layer, program="ifn_gamma_response")
        locator = f"cohort_label={cohort}; comparison_layer={layer}; program=ifn_gamma_response"
        for field, suffix in [("matched_sample_n", "n"), ("spearman_rho", "rho"), ("ci_lower", "ci_low"), ("ci_upper", "ci_high"), ("q_value_bh", "q")]:
            f(f"{prefix}_{suffix}", row[field], apollo_path, locator, field, "direct cell")

    f("scrna_donors", support["single_cell"]["unique_donor_units"], "results/support_layers_summary_v2.json", ".single_cell.unique_donor_units", "value", "JSON field")
    f("scrna_sources", support["single_cell"]["unique_source_ids"], "results/support_layers_summary_v2.json", ".single_cell.unique_source_ids", "value", "JSON field")
    f("scrna_pooled_rows", 24, "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv", "all rows", "row count", "count")
    f("scrna_evaluable", int(scrna["status"].eq("evaluable").sum()), scrna_path, "status=evaluable", "row count", "count")
    support_mask = scrna["within_source_support_window"].astype(bool)
    f("scrna_supported", int(support_mask.sum()), scrna_path, "within_source_support_window=True", "row count", "count")
    f("scrna_positive", int((support_mask & (scrna["median_delta"] > 0)).sum()), scrna_path, "support=True; median_delta>0", "row count", "count")
    f("scrna_negative", int((support_mask & (scrna["median_delta"] < 0)).sum()), scrna_path, "support=True; median_delta<0", "row count", "count")
    sensitivity_mask = scrna["n_ge_10_sensitivity_window"].astype(bool)
    f("scrna_sensitivity_supported", int(sensitivity_mask.sum()), scrna_path, "n_ge_10_sensitivity_window=True", "row count", "count")
    f("scrna_sensitivity_positive", int((sensitivity_mask & (scrna["median_delta"] > 0)).sum()), scrna_path, "n_ge_10_sensitivity_window=True; median_delta>0", "row count", "count")
    f("scrna_sensitivity_negative", int((sensitivity_mask & (scrna["median_delta"] < 0)).sum()), scrna_path, "n_ge_10_sensitivity_window=True; median_delta<0", "row count", "count")
    score_table = pd.read_csv(root / "results/scrna_combined/public_scrna_and_3ca_donor_compartment_scores_v2.csv")
    score_table["donor_unit"] = score_table["source_id"].astype(str) + "::" + score_table["cancer_code"].astype(str) + "::" + score_table["donor_id"].astype(str)
    method_counts = {"exact": 0, "asymptotic": 0, "permutation": 0, "all_zero": 0}
    for (_, _, _), group in score_table.groupby(["source_id", "cancer_code", "signature"], sort=False):
        wide = group.pivot(index="donor_unit", columns="compartment", values="signature_score")
        for left, right in [("malignant", "immune"), ("malignant", "stromal"), ("immune", "stromal")]:
            if left not in wide or right not in wide:
                continue
            delta_values = (wide[left] - wide[right]).dropna().to_numpy(dtype=float)
            n = len(delta_values)
            if n < 6:
                continue
            if np.allclose(delta_values, 0):
                method_counts["all_zero"] += 1
                continue
            has_zero = bool(np.any(delta_values == 0))
            has_ties = len(np.unique(np.abs(delta_values))) < n
            if n > 50:
                method_counts["asymptotic"] += 1
            elif not (has_zero or has_ties):
                method_counts["exact"] += 1
            elif n <= 13:
                method_counts["permutation"] += 1
            else:
                method_counts["asymptotic"] += 1
    for key, value in method_counts.items():
        f(f"wilcoxon_actual_{key}", value, "results/scrna_combined/public_scrna_and_3ca_donor_compartment_scores_v2.csv", "all source-cancer-program paired contrasts with n>=6", "resolved method", "SciPy 1.18.0 method='auto' decision rule")
    f("wilcoxon_actual_total", sum(method_counts.values()), scrna_path, "status=evaluable", "row count", "sum resolved methods")
    gbm_source = "999f2a15-3d7e-440b-96ae-2c806799c08c"
    for program, prefix in [("myeloid_inflammatory_context", "scrna_gbm_myeloid"), ("checkpoint_exhaustion_context", "scrna_gbm_checkpoint"), ("proteasome_antigen_processing", "scrna_gbm_proteasome")]:
        row = row_one(scrna, source_id=gbm_source, cancer_code="gbm", program=program, contrast="malignant_minus_immune")
        locator = f"source_id={gbm_source}; cancer_code=gbm; program={program}; contrast=malignant_minus_immune"
        for field, suffix in [("paired_donor_n", "n"), ("median_delta", "delta"), ("q_value_bh", "q")]:
            f(f"{prefix}_{suffix}", row[field], scrna_path, locator, field, "direct cell")
    f("spatial_rows", len(spatial), spatial_path, "all rows", "row count", "count")
    for label, contrast, prefix in [("MHC-I", "immune_minus_stromal", "spatial_mhci"), ("Cytolytic", "malignant_minus_immune", "spatial_cytolytic")]:
        row = row_one(spatial, signature_label=label, contrast=contrast)
        locator = f"signature_label={label}; contrast={contrast}"
        f(f"{prefix}_n", row["paired_sample_n"], spatial_path, locator, "paired_sample_n", "direct cell")
        f(f"{prefix}_delta", row["median_delta"], spatial_path, locator, "median_delta", "direct cell")

    return facts


def display(value: Any, rule: str) -> str:
    if rule == "int":
        return str(int(round(float(value))))
    if rule == "comma_int":
        return f"{int(round(float(value))):,}"
    if rule == "fixed3":
        return f"{float(value):.3f}"
    if rule == "fixed2":
        return f"{float(value):.2f}"
    if rule == "fixed10":
        return f"{float(value):.10f}"
    if rule == "sci3":
        return f"{float(value):.2e}"
    if rule == "identity":
        return str(value)
    raise ValueError(rule)


def numeric_tokens(text: str) -> list[str]:
    tokens = []
    for match in NUMERIC_RE.finditer(text):
        # Exclude the numeric tail of alphanumeric specimen/project identifiers
        # such as C3L-01043; signed numeric values remain included.
        if match.start() > 0 and text[match.start() - 1] == "-":
            prefix = re.search(r"([A-Za-z0-9_]+)-$", text[: match.start()])
            if prefix and re.search(r"[A-Za-z]", prefix.group(1)) and re.search(r"\d", prefix.group(1)):
                continue
        tokens.append(match.group(0))
    return tokens


class Document:
    def __init__(self, facts: dict[str, dict[str, Any]]):
        self.lines: list[str] = []
        self.maps: list[dict[str, Any]] = []
        self.facts = facts

    def add(self, text: str = "", mappings: Iterable[tuple[str, str]] = ()) -> None:
        self.lines.append(text)
        line_no = len(self.lines)
        tokens = numeric_tokens(text)
        specs = list(mappings)
        if len(tokens) != len(specs):
            raise RuntimeError(f"Line {line_no} numeric-map mismatch: tokens={tokens}, specs={specs}, text={text}")
        for ordinal, (token, (fact_id, rule)) in enumerate(zip(tokens, specs), start=1):
            if fact_id not in self.facts:
                raise KeyError(fact_id)
            expected = display(self.facts[fact_id]["value"], rule)
            if token != expected:
                raise RuntimeError(f"Line {line_no} token {token!r} != fact {fact_id} display {expected!r}")
            fact = self.facts[fact_id]
            self.maps.append(
                {
                    "occurrence_id": f"L{line_no:04d}_N{ordinal:02d}",
                    "manuscript_line": line_no,
                    "token_ordinal": ordinal,
                    "numeric_token": token,
                    "numeric_id": fact_id,
                    "source_value": fact["value"],
                    "display_rule": rule,
                    "source_file": fact["source_file"],
                    "source_locator": fact["source_locator"],
                    "source_field": fact["source_field"],
                    "derivation": fact["derivation"],
                    "verification_status": "TRACEABLE_MACHINE_CHECK",
                }
            )


def build_manuscript(root: Path, facts: dict[str, dict[str, Any]]) -> tuple[str, list[dict[str, Any]]]:
    d = Document(facts)
    d.add("# Methods and Results")
    d.add()
    d.add("## Methods")
    d.add()
    d.add("### Study Design and Evidence Boundary")
    d.add()
    d.add(
        "We performed a cross-sectional proteogenomic concordance analysis in 10 cancer types represented in processed Clinical Proteomic Tumor Analysis Consortium data. The independent unit was the patient or case. The study had no longitudinal time origin, clinical endpoint, treatment exposure, or post-baseline covariate. It estimated within-cancer concordance between RNA and protein or phosphoprotein immune-ecology program scores; it was not designed to establish prognosis, prediction, diagnosis, treatment response, causality, mediation, malignant-cell-intrinsic origin, post-transcriptional mechanism, or therapeutic vulnerability.",
        [("cancer_count", "int")],
    )
    d.add()
    d.add("### Sources, Case Mapping, and Inference Units")
    d.add()
    d.add(
        "The primary layer used processed RNA, proteome, phosphoproteome, ESTIMATE, xCell, and CIBERSORT matrices for breast cancer, clear-cell renal cell carcinoma, colon adenocarcinoma, glioblastoma, head and neck squamous carcinoma, lung squamous carcinoma, lung adenocarcinoma, ovarian cancer, pancreatic ductal adenocarcinoma, and endometrial cancer. Files were identified by provider, cancer, data type, filename, byte count, and cryptographic digest; the source-provenance manifest was documentary and postdated result generation. No file-level public URL was available in that manifest for the core CPTAC matrices.")
    d.add()
    d.add(
        "RNA aliquots labelled tumor were retained; unsuffixed RNA identifiers were considered tumor-like only in cohorts whose processed RNA matrix used unsuffixed tumor identifiers. Adjacent, normal, and quality-control identifiers were excluded. Terminal tumor, adjacent, normal, and technical suffixes were removed to form a normalized case identifier. When multiple eligible aliquots mapped to one case, an explicit tumor suffix was preferred and remaining ties were resolved by lexical sample identifier. RNA and each omic layer were joined one-to-one by normalized case, with uniqueness required within cancer, program, layer, and case."
    )
    d.add()
    d.add("### Program Definition, Identifier Mapping, and Scoring")
    d.add()
    d.add(
        "The immutable configuration contained 8 programs and is distributed as a byte-identical configuration table with its version and SHA256 identity. The programs and complete gene memberships were: MHC-I antigen presentation (HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2, TAPBP, NLRC5, PSMB8, PSMB9, PSMB10, ERAP1, ERAP2, CALR, CANX, PDIA3); MHC-II antigen presentation (HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DMA, HLA-DMB, CD74, CIITA, LGMN, CTSS); interferon-gamma response (IFNG, STAT1, IRF1, IRF9, JAK1, JAK2, CXCL9, CXCL10, CXCL11, GBP1, GBP2, GBP5, IDO1, ISG15, IFIT1, IFIT2, IFIT3, MX1, OAS1, OAS2, OAS3); cytolytic T-cell context (CD8A, CD8B, GZMA, GZMB, GZMH, GZMK, PRF1, NKG7, GNLY, CTSW, CST7, CCL5, CXCR3, TRAC, TRBC1, TRBC2); checkpoint/exhaustion context (PDCD1, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT, TOX, ENTPD1, CD80, CD86, ICOS, TNFRSF9); myeloid inflammatory context (LST1, LYZ, TYROBP, AIF1, FCGR3A, FCGR1A, CSF1R, ITGAM, CD68, CD163, MSR1, C1QA, C1QB, C1QC, IL1B, TNF, CXCL8, S100A8, S100A9); TGF/EMT exclusion context (TGFB1, TGFBR1, TGFBR2, SMAD2, SMAD3, VIM, FN1, COL1A1, COL1A2, COL3A1, ACTA2, TAGLN, SERPINE1, THBS1, ITGA5, ITGB1, LOXL2, ZEB1, SNAI1, SNAI2); and proteasome antigen processing (PSMB8, PSMB9, PSMB10, PSME1, PSME2, PSME3, PSMA1, PSMA2, PSMA3, PSMA4, PSMA5, PSMA6, PSMA7, PSMB1, PSMB2, PSMB3, PSMB4, PSMB5, PSMB6, PSMB7, UBC, UBB).",
        [("program_count", "int")],
    )
    d.add()
    d.add(
        "For CPTAC, gene symbols were stripped and uppercased. Stable Ensembl identifiers were extracted as ENSG followed by digits with an optional dot-version removed; every matrix feature row matching any stable Ensembl identifier listed for a program gene was retained. Within each cancer and matrix, each feature was standardized across eligible tumors using the sample standard deviation (ddof=1); a zero or nonfinite standard deviation yielded missing standardized values. Standardized feature rows mapped to the same gene, including phosphosites, were averaged within sample before equal-weight averaging across recovered genes. A sample score required at least 3 finite gene-level values. The same gene-equal rule was used for RNA, and no primary score was imputed.",
        [("cptac_ddof", "int"), ("min_score_genes", "int")],
    )
    d.add()
    d.add("### Primary Estimand, Multiplicity, and Heterogeneity")
    d.add()
    d.add(
        "The unchanged primary family was 10 cancers x 8 programs x 2 comparison layers, yielding 160 rows. The estimand was the within-cancer patient-level Spearman correlation between matched RNA and omic scores. A row required at least 20 complete cases and at least 3 distinct values in each score.",
        [("cancer_count", "int"), ("program_count", "int"), ("layer_count", "int"), ("primary_registry_rows", "int"), ("primary_min_n", "int"), ("primary_min_distinct_values", "int")],
    )
    d.add(
        "Spearman correlations used SciPy. Confidence intervals used paired-case BCa bootstrap with 20,000 resamples and tail probabilities 0.025 and 0.975. Bias correction used the proportion of bootstrap estimates below the observed rho with half weight for equality; acceleration used leave-one-case-out jackknife estimates. For rows with more than 30 complete cases and no ties, the two-sided Spearman P value was SciPy's asymptotic rank-correlation test. Rows with no more than 30 cases or any score ties used 999,999 random permutations of the matched omic ranks, counted absolute permuted rho at least as large as observed rho, and applied the plus-one correction.",
        [("primary_bootstrap_resamples", "comma_int"), ("bca_alpha_low", "fixed3"), ("bca_alpha_high", "fixed3"), ("primary_asymptotic_n_cutoff", "int"), ("primary_asymptotic_n_cutoff", "int"), ("primary_permutation_resamples", "comma_int")],
    )
    d.add(
        "Bootstrap and permutation generators used NumPy's default random generator. The seed was the little-endian integer represented by the first 8 bytes of the SHA256 digest of the project identifier, analysis component, cancer, program, and layer joined with double colons; the analysis component was primary-ci for intervals and primary-p for P values. Benjamini-Hochberg adjustment was applied once across all finite primary P values while retaining the complete 160-row registry.",
        [("seed_digest_bytes", "int"), ("primary_registry_rows", "int")],
    )
    d.add(
        "Descriptive heterogeneity was calculated on Fisher-z-transformed cancer-specific rho values with inverse-variance weights equal to complete-pair count minus three. The fixed weighted mean z was used in Cochran Q; I2 was max(0, [Q minus the number of contributing cancers plus one]/Q), with I2 set to zero when Q was zero. These summaries were descriptive; no pooled pan-cancer coefficient was a primary estimand.",
        [("i2_lower_bound", "int")],
    )
    d.add()
    d.add("### Sensitivity and RNA-Adjusted Context Analyses")
    d.add()
    d.add(
        "The feature-row sensitivity averaged standardized feature rows directly. A second sensitivity lowered the per-sample scoring requirement to one recovered gene without changing the primary registry. Leave-one-cancer-out medians and within-cancer-standardized pooled correlations were descriptive only."
    )
    d.add(
        "Within each primary stratum meeting the analysis gate, ordinary least squares fit omic score as intercept plus beta times RNA score; linear residual was observed minus fitted omic score. A rank sensitivity replaced both scores with average ranks before the same regression. Residuals were joined by case to tumor purity, ESTIMATE immune score, CIBERSORT T-cell signal, and CIBERSORT myeloid signal. Each association required at least 20 complete pairs and at least 3 distinct values in both variables. Two-sided asymptotic Spearman tests were adjusted in separate Benjamini-Hochberg families for linear and rank residuals.",
        [("primary_min_n", "int"), ("residual_min_distinct_values", "int")],
    )
    d.add()
    d.add("### Public Supporting Layers")
    d.add()
    d.add(
        "For TCGA, RNA-expression rows with primary-solid-tumor sample-type code 1 in the 10 target cancers were sorted by patient and sample identifier, and the lexically first eligible sample was retained per patient; the complete selected-sample table is distributed with the analysis. Candidate expression-row identifiers were the union of program gene symbols and their DepMap-derived Entrez identifiers. A program was emitted only when at least 3 mapped gene symbols were represented by present rows. Within cancer, selected rows were standardized across patients with ddof=1; zero or nonfinite standard deviations became missing. Each patient score was the unweighted mean of finite standardized selected rows, with no additional per-patient finite-row minimum and no imputation. Within-cancer program-to-immune-signature correlations required n>=25 and used the two-sided SciPy Spearman asymptotic P value: average ranks, Pearson correlation of ranks, and the Student t transformation with n minus two degrees of freedom. Benjamini-Hochberg adjustment was applied once across the full 8 x 68 x 10 matrix.",
        [("tcga_primary_type_code", "int"), ("cancer_count", "int"), ("tcga_min_program_genes", "int"), ("tcga_ddof", "int"), ("tcga_min_correlation_n", "int"), ("program_count", "int"), ("tcga_immune_signature_count", "int"), ("cancer_count", "int")],
    )
    d.add(
        "The executed fixed anchor map contained 17 program-anchor links: MHC-I to MHC1_21978456 and IFNG_score_21050467; MHC-II to MHC2_21978456 and TcClassII_score; IFN-gamma to IFNG_score_21050467, IFN_21978456, and Module3_IFN_score; Cytolytic to CD8A and Tcell_21978456; Checkpoint to PD1_PDL1_score and Tcell_21978456; Myeloid to CD68 and CSF1_response; TGF/EMT to TGFB_score_21050467 and TGFB_PCA_17349583; and Proteasome to MHC1_21978456 and IFN_21978456. For each anchor, cancer-specific correlations were reduced to the median rho and positive-cancer fraction; a program's selected-anchor median was the median of its anchor-level median rhos, and its selected-anchor positive fraction was the minimum anchor-level positive-cancer fraction. Immune-subtype Kruskal-Wallis tests retained groups with at least 20 patients and were adjusted across the 8 program tests separately.",
        [("tcga_anchor_rule_rows", "int"), ("tcga_min_subtype_n", "int"), ("tcga_subtype_tests", "int")],
    )
    d.add(
        "For APOLLO LUAD RNA, sample count tables were outer-unioned, absent gene/sample combinations were set to zero, log2(count+1) was applied, and duplicated uppercase gene symbols were averaged. For LUAD proteomics, LUAD phosphoproteomics, ovarian RNA, and ovarian proteomics, duplicated uppercase gene symbols or phosphosite-to-gene rows were likewise averaged before standardization. Each collapsed gene was then standardized across samples with ddof=0; zero or nonfinite standard deviations became missing. A program was emitted only when at least 3 genes were present in that layer matrix. Its sample score was the unweighted mean of finite standardized recovered-gene values, using default missing-value skipping; there was no additional per-sample three-finite-gene rule, an all-missing vector yielded a missing score, and no score was imputed. RNA and omic scores were joined one-to-one by public sample identifier, and each row required n>=10 and at least 3 distinct values per score. The registry contained 3 cohort-layer families of 8 programs. P values used the same two-sided SciPy Spearman asymptotic algorithm; confidence intervals were tanh[atanh(rho) +/- 1.9599639845/sqrt(n-3)]. Benjamini-Hochberg adjustment was applied within each cohort-layer family.",
        [("apollo_log_pseudocount", "int"), ("apollo_ddof", "int"), ("apollo_min_program_genes", "int"), ("apollo_min_n", "int"), ("apollo_min_distinct_values", "int"), ("apollo_family_count", "int"), ("apollo_family_rows", "int"), ("apollo_fisher_z_constant", "fixed10"), ("fisher_z_n_offset", "int")],
    )
    d.add(
        "Single-cell sources were reduced to source-specific donor-by-compartment pseudobulk profiles. Exact field-selection priorities, cancer mappings, compartment keyword dictionaries, dataset-specific overrides, and the executed priority order are distributed in the mapping-rule artifact. The streaming H5AD sources required at least 20 cells per donor-compartment group and at least 3 recovered program genes; the 3CA sources independently used a 20-cell group gate and a 3-gene program-recovery gate. Counts were aggregated, converted to counts per million, log1p transformed, and duplicate uppercase gene symbols were averaged. Genes were standardized within cancer with ddof=0 and zero/nonfinite variance set to missing. Each donor-compartment program score was the unweighted mean of finite standardized recovered-gene values with default missing-value skipping; there was no additional per-group three-finite-gene rule, an all-missing vector yielded a missing score, and no score was imputed. Malignant-minus-immune, malignant-minus-stromal, and immune-minus-stromal differences were paired within donor, source, cancer, and program.",
        [("scrna_cellxgene_min_group_cells", "int"), ("scrna_cellxgene_min_program_genes", "int"), ("scrna_3ca_min_group_cells", "int"), ("scrna_3ca_min_program_genes", "int"), ("support_ddof", "int")],
    )
    d.add(
        "Two-sided paired Wilcoxon tests used SciPy with zero_method=wilcox, method=auto, and no continuity correction; all-zero difference vectors were assigned P=1. Under the executed SciPy version, auto selected exact calculation for n<=50 without ties or zeros, exhaustive sign permutations for tied or zero-containing n<=13, and the asymptotic normal method otherwise. Among 256 tested rows, 224 resolved to exact and 32 to asymptotic; none used the tied small-sample permutation or all-zero special case. Benjamini-Hochberg adjustment was applied within each source-cancer family of 24 program-by-contrast rows. Support required paired donor n>=6 and q<0.05; n>=10 was retained only as a sensitivity.",
        [("wilcoxon_all_zero_p", "int"), ("wilcoxon_auto_exact_max_n", "int"), ("wilcoxon_auto_tied_permutation_max_n", "int"), ("wilcoxon_actual_total", "int"), ("wilcoxon_actual_exact", "int"), ("wilcoxon_actual_asymptotic", "int"), ("scrna_family_rows", "int"), ("scrna_min_paired_n", "int"), ("scrna_support_q_cutoff", "fixed2"), ("scrna_sensitivity_n", "int")],
    )
    d.add(
        "Spatial breast-cancer spots were assigned by the argmax of summed malignant, immune, and stromal source-score columns when the largest sum was positive. If all sums were nonpositive, source Classification values mapped stroma to stromal, values containing lymphocytes to immune, and values containing invasive cancer to malignant; remaining spots were excluded as other. Within each sample, counts were converted to counts per million and log1p transformed. A program was emitted when at least 3 genes were present; each gene was standardized across eligible spots with ddof=0 and zero/nonfinite variance set to missing. Each spot score was the unweighted mean of its finite standardized recovered-gene values; no additional per-spot three-finite-gene rule was imposed, an all-missing vector yielded a missing score, and no score was imputed. Sample-compartment medians skipped missing spot scores. Contrasts used 4 paired sample medians; no spot-level P value was calculated.",
        [("spatial_min_program_genes", "int"), ("support_ddof", "int"), ("spatial_sample_n", "int")],
    )
    d.add()
    d.add("### Software")
    d.add()
    d.add(
        "The retained execution environment used Python 3.14.6, NumPy 2.5.1, pandas 3.0.5, SciPy 1.18.0, h5py 3.16.0, and Matplotlib 3.11.1. Core statistical calls were scipy.stats.spearmanr, scipy.stats.wilcoxon, scipy.stats.kruskal, scipy.stats.rankdata with average ties, numpy.linalg.lstsq with rcond=None, and numpy.random.default_rng.",
        [("software_python", "identity"), ("software_numpy", "identity"), ("software_pandas", "identity"), ("software_scipy", "identity"), ("software_h5py", "identity"), ("software_matplotlib", "identity")],
    )
    d.add()
    d.add("## Results")
    d.add()
    d.add("### CPTAC Case Flow and Primary Concordance")
    d.add()
    d.add(
        "The one-to-one matched table contained 1023 unique cases across either protein layer. The 160-row registry retained 129 rows meeting analysis criteria and 31 that did not; all 31 had fewer than 20 complete score pairs under the 3-gene rule. There were no duplicated cancer-program-layer-case keys. Maximum matched case counts across cancer-layer combinations ranged from 82 to 121.",
        [("primary_unique_cases", "int"), ("primary_registry_rows", "int"), ("primary_evaluable_rows", "int"), ("primary_not_evaluable_rows", "int"), ("primary_not_evaluable_rows", "int"), ("primary_min_n", "int"), ("min_score_genes", "int"), ("matched_min", "int"), ("matched_max", "int")],
    )
    d.add(
        "Of 77 RNA-protein rows meeting analysis criteria, 76 had q<0.05; rho ranged from 0.012 to 0.950, with median 0.755. All 52 RNA-phosphoprotein rows meeting analysis criteria had q<0.05; rho ranged from 0.206 to 0.842, with median 0.566. The sole primary row not passing Benjamini-Hochberg adjustment was PDAC Proteasome RNA-protein concordance (n=105, rho=0.012, P=0.902, q=0.902).",
        [("protein_evaluable", "int"), ("protein_q_lt_005", "int"), ("primary_q_cutoff", "fixed2"), ("protein_rho_min", "fixed3"), ("protein_rho_max", "fixed3"), ("protein_rho_median", "fixed3"), ("phosphoprotein_evaluable", "int"), ("primary_q_cutoff", "fixed2"), ("phosphoprotein_rho_min", "fixed3"), ("phosphoprotein_rho_max", "fixed3"), ("phosphoprotein_rho_median", "fixed3"), ("nonsig_n", "int"), ("nonsig_rho", "fixed3"), ("nonsig_p", "fixed3"), ("nonsig_q", "fixed3")],
    )
    d.add(
        "All summaries with contributing cancers had a positive-cancer fraction of 1.000, while descriptive I2 fractions ranged from 0.352 to 0.878. Two registered phosphoproteomic program-layer combinations had no cancer meeting the complete-pair gate: MHC-II and Cytolytic.",
        [("positive_fraction_all", "fixed3"), ("i2_min", "fixed3"), ("i2_max", "fixed3")],
    )
    d.add()
    d.add("| Layer | Program | Registered cancers | Cancers meeting criteria | Complete pairs | Median rho | Positive fraction | I2 fraction |")
    d.add("|:--|:--|--:|--:|--:|--:|--:|--:|")
    table = pd.read_csv(root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv")
    for ridx, row in table.iterrows():
        key = f"table_{ridx+1:02d}"
        parts = [str(row.comparison_layer), str(row.program_label)]
        mappings: list[tuple[str, str]] = []
        for field, rule in [
            ("registered_cancer_rows", "int"), ("evaluable_cancer_n", "int"),
            ("total_complete_pairs_across_strata", "int"), ("median_within_cancer_rho", "fixed3"),
            ("positive_cancer_fraction", "fixed3"), ("descriptive_i2", "fixed3"),
        ]:
            value = row[field]
            if pd.isna(value):
                parts.append("NA")
            else:
                fid = f"{key}_{field}"
                parts.append(display(facts[fid]["value"], rule))
                mappings.append((fid, rule))
        d.add("| " + " | ".join(parts) + " |", mappings)
    d.add()
    d.add("### Scoring and Influence Sensitivities")
    d.add()
    d.add(
        "Feature-row and gene-equal estimates had the same sign in all 129 jointly available rows. Protein estimates were identical because those matrices contained one recovered feature per mapped gene. Among 52 phosphoprotein rows, the median absolute rho difference was 0.019 and the maximum was 0.265. The one-gene threshold sensitivity yielded 158 testable rows, of which 157 passed its separate Benjamini-Hochberg family, without changing the primary registry. The largest leave-one-case change in primary rho was 0.118 in the GBM Checkpoint RNA-phosphoprotein stratum after omitting case C3L-01043; the full-stratum result remained positive and supported (n=26, rho=0.506, 95% CI 0.062 to 0.784, q=0.010).",
        [("sensitivity_joined", "int"), ("sensitivity_phospho_joined", "int"), ("sensitivity_phospho_median_delta", "fixed3"), ("sensitivity_phospho_max_delta", "fixed3"), ("threshold_evaluable", "int"), ("threshold_supported", "int"), ("max_influence_delta", "fixed3"), ("max_influence_n", "int"), ("max_influence_rho", "fixed3"), ("ci_percent", "int"), ("max_influence_ci_low", "fixed3"), ("max_influence_ci_high", "fixed3"), ("max_influence_q", "fixed3")],
    )
    d.add(
        "Across 129 leave-one-cancer-out summaries, median rho values ranged from 0.317 to 0.858 and the minimum retained positive-cancer fraction was 1.000. The within-cancer-standardized pooled sensitivity had finite coefficients for 14 of 16 registered program-layer combinations; all 14 were positive, ranging from 0.319 to 0.866. Phosphoproteomic MHC-II and Cytolytic remained non-estimable because no cancer met the primary complete-pair gate.",
        [("loco_rows", "int"), ("loco_median_rho_min", "fixed3"), ("loco_median_rho_max", "fixed3"), ("loco_positive_fraction_min", "fixed3"), ("pooled_finite_rows", "int"), ("pooled_registered_rows", "int"), ("pooled_positive_rows", "int"), ("pooled_rho_min", "fixed3"), ("pooled_rho_max", "fixed3")],
    )
    d.add()
    d.add("### RNA-Adjusted Residual Context")
    d.add()
    d.add(
        "Linear residual models were available in 129 cancer-program-layer strata. The residual-context registry contained 640 rows; 516 had finite linear-residual tests and 516 had finite rank-residual tests. 35 rows passed the linear-residual family and 52 passed the rank sensitivity family; all 35 passing both had concordant directions. Linear-residual supported rows were distributed across CIBERSORT myeloid (11), CIBERSORT T-cell (5), ESTIMATE immune score (10), and tumor purity (9). These associations remained secondary descriptions of immune and purity context after RNA adjustment.",
        [("residual_strata", "int"), ("residual_registry", "int"), ("residual_linear_finite", "int"), ("residual_rank_finite", "int"), ("residual_linear_supported", "int"), ("residual_rank_supported", "int"), ("residual_both_supported", "int"), ("residual_myeloid_supported", "int"), ("residual_tcell_supported", "int"), ("residual_estimate_supported", "int"), ("residual_purity_supported", "int")],
    )
    d.add()
    d.add("### RNA-Only TCGA Sensitivity")
    d.add()
    d.add(
        "The TCGA layer included 4780 primary-tumor patients and 4780 selected samples. The full evidence surface contained 5440 within-cancer correlations, of which 4850 had q<0.05 in the single family. Across the seven non-proteasome programs, canonical-anchor median correlations ranged from 0.838 to 0.963 and the minimum positive-cancer fraction was 1.000. The proteasome comparator had anchor median rho 0.337. Immune-subtype testing generated 8 tests in a separate family; Immune C5 contained only 4 patients and remained descriptive. These findings support RNA-level immune-landscape alignment only.",
        [("tcga_patients", "int"), ("tcga_samples", "int"), ("tcga_corr_rows", "int"), ("tcga_supported", "int"), ("tcga_q_cutoff", "fixed2"), ("tcga_anchor_min", "fixed3"), ("tcga_anchor_max", "fixed3"), ("tcga_anchor_positive", "fixed3"), ("tcga_proteasome", "fixed3"), ("tcga_subtype_tests", "int"), ("tcga_c5_n", "int")],
    )
    d.add()
    d.add("### APOLLO Same-Cancer Support")
    d.add()
    d.add(
        "APOLLO contributed 87 LUAD sample units with RNA, protein, and phosphoprotein measurements and 70 ovarian sample units with RNA and protein measurements. Across 24 registered rows, 17 met analysis criteria: 11 were supported, 1 had a significant opposite direction, and 5 were null; 7 did not meet the score-recovery or n>=10 gate. Representative supported rows were LUAD IFN-gamma proteomics (rho=0.894, 95% CI 0.842 to 0.929, q=1.42e-30) and ovarian IFN-gamma proteomics (rho=0.907, 95% CI 0.853 to 0.941, q=2.16e-26). The opposite row was retained: LUAD IFN-gamma phosphoproteomics (rho=-0.281, 95% CI -0.464 to -0.075, q=0.042). This heterogeneous mixture supports bounded same-cancer proteomic context, not donor-proven independent replication.",
        [("apollo_luad_ifng_n", "int"), ("apollo_ov_ifng_n", "int"), ("apollo_registered_rows", "int"), ("apollo_evaluable", "int"), ("apollo_supported", "int"), ("apollo_opposite", "int"), ("apollo_null", "int"), ("apollo_not_evaluable", "int"), ("apollo_min_n", "int"), ("apollo_luad_ifng_rho", "fixed3"), ("apollo_ci_percent", "int"), ("apollo_luad_ifng_ci_low", "fixed3"), ("apollo_luad_ifng_ci_high", "fixed3"), ("apollo_luad_ifng_q", "sci3"), ("apollo_ov_ifng_rho", "fixed3"), ("apollo_ci_percent", "int"), ("apollo_ov_ifng_ci_low", "fixed3"), ("apollo_ov_ifng_ci_high", "fixed3"), ("apollo_ov_ifng_q", "sci3"), ("apollo_luad_ifng_opposite_rho", "fixed3"), ("apollo_ci_percent", "int"), ("apollo_luad_ifng_opposite_ci_low", "fixed3"), ("apollo_luad_ifng_opposite_ci_high", "fixed3"), ("apollo_luad_ifng_opposite_q", "fixed3")],
    )
    d.add()
    d.add("### Donor-Aware Single-Cell and Sample-Aware Spatial Triangulation")
    d.add()
    d.add(
        "The combined single-cell and 3CA layer contained 345 source-specific donor units from 10 public sources. No cross-repository participant crosswalk established global biological independence, so the 24 pooled program-by-compartment rows were descriptive and had no pooled P value or support claim. Within source-cancer families, 256 rows met analysis criteria and 116 met the separately adjusted n>=6 support rule; 42 supported rows had positive median differences and 74 had negative median differences. Representative source-bounded rows were GBM Myeloid malignant-minus-immune (paired n=68, median difference=-1.661, q=5.24e-12), GBM Checkpoint malignant-minus-immune (paired n=68, median difference=-0.822, q=5.24e-12), and GBM Proteasome malignant-minus-immune (paired n=68, median difference=0.484, q=5.24e-12). These contrasts locate signal across broad compartments within sources but do not establish malignant-cell origin or mechanism.",
        [("scrna_donors", "int"), ("scrna_sources", "int"), ("scrna_pooled_rows", "int"), ("scrna_evaluable", "int"), ("scrna_supported", "int"), ("scrna_min_paired_n", "int"), ("scrna_positive", "int"), ("scrna_negative", "int"), ("scrna_gbm_myeloid_n", "int"), ("scrna_gbm_myeloid_delta", "fixed3"), ("scrna_gbm_myeloid_q", "sci3"), ("scrna_gbm_checkpoint_n", "int"), ("scrna_gbm_checkpoint_delta", "fixed3"), ("scrna_gbm_checkpoint_q", "sci3"), ("scrna_gbm_proteasome_n", "int"), ("scrna_gbm_proteasome_delta", "fixed3"), ("scrna_gbm_proteasome_q", "sci3")],
    )
    d.add(
        "Under the paired-donor n>=10 sensitivity, 85 source-cancer-program-contrast rows remained supported: 31 had positive and 54 had negative median differences.",
        [("scrna_sensitivity_n", "int"), ("scrna_sensitivity_supported", "int"), ("scrna_sensitivity_positive", "int"), ("scrna_sensitivity_negative", "int")],
    )
    d.add(
        "All 4 breast Visium samples yielded sample-level program and compartment summaries. The pairwise table contained 24 program-contrast rows, each with paired sample n=4. Examples were immune-minus-stromal MHC-I median difference 0.299 and malignant-minus-immune Cytolytic median difference -0.348. These were same-collection sample-level contrasts without spot-level P values or an independent-validation claim.",
        [("spatial_sample_n", "int"), ("spatial_rows", "int"), ("spatial_mhci_n", "int"), ("spatial_mhci_delta", "fixed3"), ("spatial_cytolytic_delta", "fixed3")],
    )
    text = "\n".join(d.lines) + "\n"
    return text, d.maps


def attach_source_hashes(root: Path, rows: list[dict[str, Any]]) -> None:
    cache: dict[str, str] = {}
    for row in rows:
        source = str(row["source_file"])
        if source not in cache:
            path = root / source
            if not path.is_file():
                raise FileNotFoundError(path)
            cache[source] = sha256_file(path)
        row["source_sha256"] = cache[source]


def write_numeric_map(root: Path, rows: list[dict[str, Any]]) -> None:
    attach_source_hashes(root, rows)
    fields = [
        "occurrence_id", "manuscript_line", "token_ordinal", "numeric_token", "numeric_id", "source_value",
        "display_rule", "source_file", "source_sha256", "source_locator", "source_field", "derivation", "verification_status",
    ]
    with (root / NUMERIC_MAP).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def adjudication_records(root: Path, snapshots: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "finding_id": "FRESH-P1-01",
            "severity": "P1",
            "allegation": "Exact program gene membership, identifier mapping/collapse, and zero-variance handling were not supplied.",
            "disposition": "ACCEPTED",
            "exact_evidence": [
                {"file": CONFIG_COPY, "locator": "139 data rows; fields signature, gene_symbol, ensembl_stable_ids_from_cptac", "fact": "Byte-identical executed program configuration; SHA256 b631094186b4f1a69f3e2ee8b948009f10f4b7ad0b644cd1cf446e4152c46ef1."},
                {"file": "code/run_cptac_v2.py", "locator": "lines 119-121, 207-234, 237-317", "fact": "Stable ENSG canonicalization, ddof=1 feature standardization, zero/nonfinite-variance-to-missing rule, within-gene collapse, equal-gene scoring, and 3-gene gate."},
                {"file": "results/cptac/cptac_gene_recovery_v2.csv", "locator": "all cancer/layer/program rows; recovered_genes and recovered_feature_row_count", "fact": "Executed per-matrix recovery surface."},
            ],
            "action": "Exported the immutable configuration, froze its hash, and added exact canonicalization, collapse, missingness, zero-variance, and scoring rules to Methods and mapping metadata.",
            "rerun_status": "NOT_REQUIRED_DOCUMENTATION_ONLY",
            "affected_artifacts": [CONFIG_COPY, RULES, P2_PATH],
            "resolution": "The executed scoring configuration is now locally reconstructible; numerical estimates were not changed.",
        },
        {
            "finding_id": "FRESH-P1-02",
            "severity": "P1",
            "allegation": "TCGA selection, single-cell/spatial mapping, and APOLLO operation order were incomplete or ambiguous.",
            "disposition": "ACCEPTED",
            "exact_evidence": [
                {"file": snapshots["tcga"]["relative_path"], "locator": "lines 117-134", "fact": "Primary tumor filter, cancer filter, patient/sample sort, and first-sample selection."},
                {"file": TCGA_SELECTION, "locator": "4780 unique patient rows", "fact": "Exact selected TCGA sample membership recovered from the retained score table."},
                {"file": "code/run_h5ad_support_v2.py", "locator": "lines 34-83, 169-232, 479-500", "fact": "Dataset defaults, keyword/override priority, and spatial argmax/classification fallback."},
                {"file": snapshots["threeca"]["relative_path"], "locator": "lines 38-71 and 133-142", "fact": "Exact 3CA dataset filters and compartment dictionaries."},
                {"file": snapshots["apollo"]["relative_path"], "locator": "lines 85-136 and 171-228", "fact": "Every APOLLO layer collapses duplicates before z-standardization; LUAD RNA performs union/fill/log transform before collapse."},
                {"file": RULES, "locator": "tcga, apollo, single_cell_cellxgene, single_cell_3ca, spatial", "fact": "Machine-readable consolidation of exact executed rules."},
            ],
            "action": "Exported selection membership, source dictionaries, spatial fallback rules, and per-layer APOLLO order; corrected Methods wording without changing execution.",
            "rerun_status": "NOT_REQUIRED_DOCUMENTATION_ONLY",
            "affected_artifacts": [TCGA_SELECTION, RULES, P2_PATH, "current_p2_closure/executed_code"],
            "resolution": "All alleged supporting-layer choices are explicit and traceable to retained executed code or row-level outputs.",
        },
        {
            "finding_id": "FRESH-P1-03",
            "severity": "P1",
            "allegation": "Software/function identities, seeds, P-value algorithms, Wilcoxon mode, and Q/I2 formula were insufficiently specified.",
            "disposition": "ACCEPTED",
            "exact_evidence": [
                {"file": ENVIRONMENT, "locator": "recorded_runtime, observed_versions_from_recorded_interpreter, functions", "fact": "Python and package identities recovered from the interpreter recorded in RUN_MANIFEST."},
                {"file": "code/run_cptac_v2.py", "locator": "lines 330-439, 453-535, 595-624", "fact": "BH, SHA256 seed derivation, BCa, Monte Carlo permutation, SciPy asymptotic branch, and Fisher-z Q/I2 weights/formula."},
                {"file": snapshots["tcga"]["relative_path"], "locator": "lines 251-288", "fact": "TCGA uses scipy.stats.spearmanr and a single BH family."},
                {"file": "code/postprocess_support_v2.py", "locator": "lines 82-87, 119-154, 241-249", "fact": "APOLLO Fisher-z interval/SciPy Spearman and single-cell scipy.stats.wilcoxon method=auto, zero_method=wilcox."},
                {"file": DERIVATIONS, "locator": "wilcoxon_actual_*", "fact": "Resolved-method audit from retained donor-compartment scores: 224 exact, 32 asymptotic, 0 permutation, 0 all-zero."},
            ],
            "action": "Added exact runtime/function identities, seed construction, Spearman algorithms, Wilcoxon auto-resolution, and Q/I2 formula; retained a method-resolution audit.",
            "rerun_status": "NO_SCIENTIFIC_RERUN; DETERMINISTIC_METHOD_CLASSIFICATION_EXECUTED",
            "affected_artifacts": [ENVIRONMENT, DERIVATIONS, P2_PATH],
            "resolution": "Bit-level algorithm choices are documented; no endpoint, family, estimate, P value, or q value changed.",
        },
        {
            "finding_id": "FRESH-P1-04",
            "severity": "P1",
            "allegation": "The numeric map did not permit value-by-value verification of all manuscript numbers.",
            "disposition": "ACCEPTED",
            "exact_evidence": [
                {"file": "current_p1_closure/NUMERIC_MAP.csv", "locator": "75 rows", "fact": "The prior map omitted multiple displayed Table-2, APOLLO, single-cell, and spatial values."},
                {"file": NUMERIC_MAP, "locator": "one row per manuscript numeric-token occurrence", "fact": "Every occurrence has line/ordinal, numeric ID, source path/hash, row/field locator, derivation, and display rule."},
                {"file": "current_p2_closure/run_current_p2_gate.py", "locator": "numeric inventory and fact comparison checks", "fact": "The executable gate rejects missing, duplicate, stale, or source-inconsistent mappings."},
            ],
            "action": "Replaced selective mapping with occurrence-level traceability and deterministic source-value verification.",
            "rerun_status": "NOT_REQUIRED_DOCUMENTATION_ONLY; NUMERIC_VALIDATION_EXECUTED",
            "affected_artifacts": [NUMERIC_MAP, DERIVATIONS, "current_p2_closure/NUMERIC_GATE_CURRENT_P2.json"],
            "resolution": "The complete manuscript numeric surface is machine-covered without changing reported values.",
        },
    ]


def write_adjudication(root: Path, records: list[dict[str, Any]]) -> None:
    payload = {
        "project": PROJECT,
        "governing_review": FRESH_AUDIT,
        "declared_inventory": {"P0": 0, "P1": 4, "P2": 0},
        "adjudicated_records": records,
        "inventory_complete": len(records) == 4 and all(x["severity"] == "P1" for x in records),
    }
    write_json(root / "current_p2_closure/CURRENT_P2_ADJUDICATION.json", payload)
    lines = [
        "# BIO-05 CURRENT_P2 adjudication",
        "",
        "Governing inventory: P0=0, P1=4, P2=0. All four alleged findings were adjudicated from local evidence.",
        "",
        "| Finding | Severity | Disposition | Evidence-based resolution | Rerun status |",
        "|:--|:--|:--|:--|:--|",
    ]
    for row in records:
        evidence = "; ".join(f"{x['file']} ({x['locator']})" for x in row["exact_evidence"])
        lines.append(f"| {row['finding_id']} | {row['severity']} | {row['disposition']} | {row['resolution']} Evidence: {evidence}. | {row['rerun_status']} |")
    lines.extend(
        [
            "",
            "No numerical contradiction was found. No scientific analysis layer was rerun; only immutable configuration/code export, deterministic method classification, manuscript generation, and validation were executed.",
            "",
        ]
    )
    (root / "current_p2_closure/CURRENT_P2_ADJUDICATION.md").write_text("\n".join(lines), encoding="utf-8")


def cross_freeze_adjudications() -> list[dict[str, Any]]:
    return [
        {
            "finding_id": "CROSS-FREEZE-P1-01",
            "severity": "P1",
            "finding": "Several numeric-map rows were value-correct but semantically bound to unrelated CPTAC constants or interval methods.",
            "disposition": "ACCEPTED_AND_CLOSED",
            "evidence": [
                {
                    "path": "code/run_cptac_v2.py",
                    "locator": "lines 468-495 and 699-709",
                    "fact": "The primary and residual distinct-value gates are separate executed operations from the MIN_GENES score-recovery rule.",
                },
                {
                    "path": "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py",
                    "locator": "zscore_rows lines 404-408",
                    "fact": "The TCGA-imported standardization helper uses ddof=1.",
                },
                {
                    "path": "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py",
                    "locator": "MIN_SUBTYPE_N line 65 and subtype gate lines 327-350",
                    "fact": "The TCGA subtype group threshold is n>=20.",
                },
                {
                    "path": "current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py",
                    "locator": "MIN_SIGNATURE_GENES line 24 and score_signature_matrix lines 105-136",
                    "fact": "The APOLLO matrix-level program recovery gate and finite-value averaging are defined in the APOLLO execution.",
                },
                {
                    "path": "code/postprocess_support_v2.py",
                    "locator": "lines 82-87, 119-154, and 287-295",
                    "fact": "APOLLO intervals are Fisher-z intervals; APOLLO distinct-value and q thresholds and single-cell support thresholds are analysis-specific.",
                },
                {
                    "path": "current_p2_closure/NUMERIC_MAP_CURRENT_P2.csv",
                    "locator": "numeric_id, source_file, source_locator for every revised occurrence",
                    "fact": "The regenerated map uses operation-specific numeric IDs and code/result authorities.",
                },
            ],
            "commands": [
                "rg -n 'nunique|ddof|MIN_SUBTYPE_N|MIN_SIGNATURE_GENES|spearman_ci|q_value_bh' code current_p2_closure/executed_code",
                "python3 current_p2_closure/build_current_p2_outputs.py --analysis-root .",
                "python3 current_p2_closure/run_current_p2_gate.py --analysis-root .",
            ],
            "changes": [
                "Replaced score-recovery IDs used for primary/residual distinct-value occurrences with operation-specific IDs.",
                "Rebound TCGA ddof and subtype thresholds to TCGA execution snapshots.",
                "Rebound APOLLO and single-cell recovery/distinct/q thresholds to their executed modules.",
                "Rebound all three APOLLO 95% interval labels to the Fisher-z implementation rather than CPTAC BCa constants.",
            ],
            "rerun": "NO_SCIENTIFIC_RERUN; deterministic document/map regeneration only",
            "resolution": "All named semantic misbindings are corrected and covered by an explicit semantic-binding gate.",
        },
        {
            "finding_id": "CROSS-FREEZE-P1-02",
            "severity": "P1",
            "finding": "The TCGA canonical-anchor rule and supporting-layer program-score aggregation/missingness rules were not fully reproducible from Methods.",
            "disposition": "ACCEPTED_AND_CLOSED",
            "evidence": [
                {
                    "path": "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap.py",
                    "locator": "SIGNATURE_ANCHORS lines 79-88",
                    "fact": "The executed helper fixes the complete program-to-anchor mapping.",
                },
                {
                    "path": "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py",
                    "locator": "meta-summary lines 290-308 and signature-summary lines 385-407",
                    "fact": "Anchor selection uses membership in the fixed map; summary fields are the median of anchor-level median rhos and minimum anchor-level positive fractions.",
                },
                {
                    "path": "current_p2_closure/TCGA_ANCHOR_SELECTION_CURRENT_P2.csv",
                    "locator": "17 anchor rows across 8 programs; membership_match and summary_match",
                    "fact": "All fixed anchors and both aggregation functions match retained TCGA summaries.",
                },
                {
                    "path": "current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py",
                    "locator": "zscore_columns and score_signature_matrix lines 85-136",
                    "fact": "APOLLO uses matrix-level recovery followed by equal-weight finite-value means with default missing-value skipping.",
                },
                {
                    "path": "code/run_h5ad_support_v2.py",
                    "locator": "single-cell lines 404-445 and spatial lines 510-543",
                    "fact": "Streaming single-cell and spatial program-score construction and all-missing behavior are explicit.",
                },
                {
                    "path": "current_p2_closure/executed_code/run_bio05_3ca_source_mapping_v0_8.py",
                    "locator": "lines 234-284",
                    "fact": "The 3CA group/recovery gates and equal-weight finite-value aggregation match the single-cell Methods clarification.",
                },
                {
                    "path": "current_p2_closure/SUPPORT_MAPPING_RULES_CURRENT_P2.json",
                    "locator": "tcga.anchor_selection and score_aggregation_and_missingness for tcga/apollo/single_cell_cellxgene/single_cell_3ca/spatial",
                    "fact": "The rules are consolidated in a deterministic machine-readable artifact.",
                },
            ],
            "commands": [
                "python3 current_p2_closure/build_current_p2_outputs.py --analysis-root .",
                "python3 current_p2_closure/run_current_p2_gate.py --analysis-root .",
            ],
            "changes": [
                "Exported the two TCGA imported helper snapshots named by the source manifest.",
                "Exported and cross-checked the 17-link TCGA anchor-selection table.",
                "Added the exact TCGA anchor identities and aggregation rule to Methods.",
                "Added matrix-level recovery, equal-weight finite-value aggregation, sample/group/spot missingness, and no-imputation rules for TCGA, APOLLO, single-cell, 3CA, and spatial layers.",
            ],
            "rerun": "NO_SCIENTIFIC_RERUN; retained aggregate/result cross-checks only",
            "resolution": "The reported supporting-layer results are now reproducible from explicit executed selection, aggregation, and missingness rules.",
        },
        {
            "finding_id": "CROSS-FREEZE-P2-01",
            "severity": "P2",
            "finding": "Named leave-one-cancer-out, within-cancer-standardized pooled, and paired-donor n>=10 sensitivities lacked Results outcomes.",
            "disposition": "ACCEPTED_AND_CLOSED_BY_REPORTING_RETAINED_RESULTS",
            "evidence": [
                {
                    "path": "results/cptac/cptac_leave_one_cancer_out_v2.csv",
                    "locator": "129 rows; leave_one_cancer_out_median_rho and leave_one_cancer_out_positive_fraction",
                    "fact": "Retained leave-one-cancer-out medians range 0.3169652855543113-0.8575531770762768 and the minimum positive fraction is 1.0.",
                },
                {
                    "path": "results/cptac/cptac_within_cancer_standardized_pooled_sensitivity_v2.csv",
                    "locator": "16 program-layer rows; within_cancer_standardized_pooled_spearman_rho",
                    "fact": "Fourteen coefficients are finite and positive, ranging 0.318787008844282-0.866281400913768; phosphoproteomic MHC-II and Cytolytic are non-estimable.",
                },
                {
                    "path": "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv",
                    "locator": "n_ge_10_sensitivity_window=True",
                    "fact": "The n>=10 sensitivity retains 85 supported rows: 31 positive and 54 negative median differences.",
                },
            ],
            "commands": [
                "python3 current_p2_closure/build_current_p2_outputs.py --analysis-root .",
                "python3 current_p2_closure/run_current_p2_gate.py --analysis-root .",
            ],
            "changes": [
                "Added minimal Results statements for all three already-executed sensitivity components.",
                "Added occurrence-level numeric-map entries for every newly reported sensitivity value.",
            ],
            "rerun": "NO_SCIENTIFIC_RERUN; direct summaries of retained result tables only",
            "resolution": "No named sensitivity component remains outcome-free in CURRENT_P2.",
        },
    ]


def render_cross_freeze_revision_md(payload: dict[str, Any]) -> str:
    lines = [
        "# BIO-05 cross-freeze revision 1",
        "",
        f"- Status: **{payload['status']}**",
        f"- Governing review: `{payload['governing_review']['path']}`",
        f"- Governing review SHA256: `{payload['governing_review']['sha256']}`",
        f"- Immutable CURRENT_P1: `{payload['immutable_current_p1']['path']}`; unchanged `{str(payload['immutable_current_p1']['unchanged']).lower()}`.",
        f"- Declared findings: P0={payload['declared_inventory']['P0']}, P1={payload['declared_inventory']['P1']}, P2={payload['declared_inventory']['P2']}.",
        f"- Unresolved after adjudication: P0={len(payload['unresolved']['P0'])}, P1={len(payload['unresolved']['P1'])}, P2={len(payload['unresolved']['P2'])}.",
        "- Scientific rerun: no.",
        "- Existing retained results changed: no.",
        "- Existing reported numerical estimates changed: no.",
        "- Retained sensitivity outcomes newly reported: yes.",
        "- Primary conclusion changed: no.",
        "",
        "## Finding-by-finding adjudication",
        "",
    ]
    for item in payload["adjudications"]:
        lines.extend(
            [
                f"### {item['finding_id']} — {item['severity']} — {item['disposition']}",
                "",
                item["finding"],
                "",
                f"Resolution: {item['resolution']}",
                "",
                "Evidence:",
                "",
            ]
        )
        for evidence in item["evidence"]:
            lines.append(f"- `{evidence['path']}` ({evidence['locator']}): {evidence['fact']}")
        lines.extend(["", "Changes:", ""])
        lines.extend(f"- {change}" for change in item["changes"])
        lines.extend(["", f"Rerun status: {item['rerun']}", ""])
    lines.extend(["## Commands", ""])
    lines.extend(f"- `{command}`" for command in payload["commands"])
    lines.extend(
        [
            "",
            "## Result and conclusion status",
            "",
            payload["result_conclusion_change_status"]["detail"],
            "",
            "## Gate and deterministic-hash status",
            "",
            f"- Numeric/semantic gate: {payload['gate_validation']['numeric_and_semantic_gate']}.",
            f"- Forbidden-language gate: {payload['gate_validation']['forbidden_language_gate']}.",
            f"- Two-pass deterministic hashes: {payload['gate_validation']['two_pass_deterministic_hashes']}.",
            "",
        ]
    )
    return "\n".join(lines)


def write_cross_freeze_revision(root: Path) -> None:
    if not CROSS_FREEZE_AUDIT.is_file() or sha256_file(CROSS_FREEZE_AUDIT) != CROSS_FREEZE_AUDIT_SHA256:
        raise RuntimeError("Cross-freeze governing response is missing or hash-mismatched")
    payload = {
        "schema_version": "1.0",
        "project": PROJECT,
        "revision_id": "CROSS_FREEZE_REVISION_1",
        "governing_review": {"path": str(CROSS_FREEZE_AUDIT), "sha256": CROSS_FREEZE_AUDIT_SHA256},
        "declared_inventory": {"P0": 0, "P1": 2, "P2": 1},
        "immutable_current_p1": {
            "path": SOURCE_P1,
            "expected_sha256": SOURCE_P1_SHA256,
            "observed_sha256": sha256_file(root / SOURCE_P1),
            "unchanged": sha256_file(root / SOURCE_P1) == SOURCE_P1_SHA256,
        },
        "before_revision": {
            "current_p2_sha256": PRE_REVISION_P2_SHA256,
            "numeric_map_sha256": PRE_REVISION_NUMERIC_MAP_SHA256,
        },
        "after_revision": {
            "current_p2_path": P2_PATH,
            "current_p2_sha256": sha256_file(root / P2_PATH),
            "numeric_map_path": NUMERIC_MAP,
            "numeric_map_sha256": sha256_file(root / NUMERIC_MAP),
        },
        "adjudications": cross_freeze_adjudications(),
        "changes": {
            "modified": [
                P2_PATH,
                "current_p2_closure/build_current_p2_outputs.py",
                "current_p2_closure/run_current_p2_gate.py",
                NUMERIC_MAP,
                DERIVATIONS,
                RULES,
                "current_p2_closure/CURRENT_P2_CLOSURE.json",
                "current_p2_closure/CURRENT_P2_CLOSURE.md",
                "current_p2_closure/COMMANDS_RUN.txt",
            ],
            "generated": [
                TCGA_ANCHOR_SELECTION,
                "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py",
                "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap.py",
                CROSS_FREEZE_REVISION_JSON,
                CROSS_FREEZE_REVISION_MD,
                "current_p2_closure/TWO_PASS_DETERMINISTIC_HASHES_CURRENT_P2.json",
            ],
            "out_of_scope_files_modified": [],
        },
        "commands": command_list(),
        "unresolved": {"P0": [], "P1": [], "P2": []},
        "result_conclusion_change_status": {
            "scientific_reanalysis_performed": False,
            "retained_result_values_changed": False,
            "existing_reported_numeric_estimates_changed": False,
            "retained_sensitivity_outcomes_added_to_results": True,
            "primary_conclusion_changed": False,
            "detail": (
                "No retained estimate, interval, P value, q value, endpoint, threshold, multiplicity family, or primary conclusion changed. "
                "CURRENT_P2 now reports already-retained leave-one-cancer-out, within-cancer-standardized pooled, and paired-donor n>=10 sensitivity outcomes."
            ),
        },
        "gate_validation": {
            "numeric_and_semantic_gate": "PENDING",
            "forbidden_language_gate": "PENDING",
            "two_pass_deterministic_hashes": "PENDING",
            "gate_path": "current_p2_closure/NUMERIC_GATE_CURRENT_P2.json",
            "hash_proof_path": "current_p2_closure/TWO_PASS_DETERMINISTIC_HASHES_CURRENT_P2.json",
        },
        "status": "PENDING_GATE",
    }
    write_json(root / CROSS_FREEZE_REVISION_JSON, payload)
    (root / CROSS_FREEZE_REVISION_MD).write_text(render_cross_freeze_revision_md(payload), encoding="utf-8")


def command_list() -> list[str]:
    return [
        expand_legacy_paths("shasum -a 256 manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md ${SCI27_RESEARCH_ROOT}/outputs/sci27_current_p2_cross_freeze_20260815/responses/BIO-05_CURRENT_P2_CROSS_FREEZE_RESPONSE.md"),
        "python3 -m py_compile current_p2_closure/build_current_p2_outputs.py current_p2_closure/run_current_p2_gate.py",
        "python3 current_p2_closure/build_current_p2_outputs.py --analysis-root .",
        "python3 current_p2_closure/run_current_p2_gate.py --analysis-root .",
        "python3 current_p2_closure/build_current_p2_outputs.py --analysis-root .",
        "python3 current_p2_closure/run_current_p2_gate.py --analysis-root .",
        "shasum -a 256 manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md manuscript/BIO-05_METHODS_RESULTS_CURRENT_P2.md current_p2_closure/NUMERIC_MAP_CURRENT_P2.csv current_p2_closure/NUMERIC_GATE_CURRENT_P2.json current_p2_closure/TWO_PASS_DETERMINISTIC_HASHES_CURRENT_P2.json cross_freeze_revision/CROSS_FREEZE_REVISION_1.json cross_freeze_revision/CROSS_FREEZE_REVISION_1.md",
    ]


def write_commands(root: Path) -> None:
    lines = [
        "# BIO-05 CURRENT_P2 closure command ledger",
        "# Commands use the cloned project root as the working directory.",
        "",
        *command_list(),
        "",
        "# The repeated build/gate sequence is the deterministic regeneration check.",
        "# No command reruns CPTAC, TCGA, APOLLO, single-cell, 3CA, or spatial scientific analyses.",
    ]
    (root / "current_p2_closure/COMMANDS_RUN.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_initial_closure(root: Path) -> None:
    p2_hash = sha256_file(root / P2_PATH)
    payload = {
        "project": PROJECT,
        "source_current_p1": {"path": SOURCE_P1, "sha256": SOURCE_P1_SHA256, "observed_sha256": sha256_file(root / SOURCE_P1), "unchanged": sha256_file(root / SOURCE_P1) == SOURCE_P1_SHA256},
        "current_p2": {"path": P2_PATH, "sha256": p2_hash},
        "commands": command_list(),
        "result_change_summary": {
            "scientific_reanalysis_performed": False,
            "numerical_results_changed": False,
            "existing_reported_numeric_estimates_changed": False,
            "retained_sensitivity_outcomes_added_to_results": True,
            "primary_conclusion_changed": False,
            "detail": "Documentation/configuration closure only. Already-retained sensitivity outcomes were added to Results; no retained estimate, confidence interval, P value, q value, endpoint, threshold, multiplicity family, or primary conclusion was changed.",
        },
        "unresolved_p0": [],
        "unresolved_p1": [],
        "unresolved_p2": [],
        "deferred_p2": [],
        "numeric_gate": {"status": "PENDING", "path": "current_p2_closure/NUMERIC_GATE_CURRENT_P2.json"},
        "forbidden_language_gate": {"status": "PENDING", "path": "current_p2_closure/NUMERIC_GATE_CURRENT_P2.json"},
        "deterministic_hash_proof": {"status": "PENDING", "path": "current_p2_closure/TWO_PASS_DETERMINISTIC_HASHES_CURRENT_P2.json"},
        "cross_freeze_revision": {"path": CROSS_FREEZE_REVISION_JSON, "status": "PENDING_GATE"},
        "claim_boundary": "Cancer-specific patient/case-level cross-sectional CPTAC RNA-protein and RNA-phosphoprotein concordance with bounded RNA-only TCGA, same-cancer APOLLO, source-bounded donor-compartment single-cell, and four-sample spatial context; no prognosis, prediction, diagnosis, treatment response, causality, mediation, mechanism, malignant-cell-intrinsic origin, therapeutic vulnerability, homogeneous pan-cancer effect, or independent-validation claim.",
        "status": "PENDING_GATE",
    }
    write_json(root / "current_p2_closure/CURRENT_P2_CLOSURE.json", payload)
    md = [
        "# BIO-05 CURRENT_P2 closure",
        "",
        f"- Source CURRENT_P1: `{SOURCE_P1}`; SHA256 `{SOURCE_P1_SHA256}`; unchanged: `{str(payload['source_current_p1']['unchanged']).lower()}`.",
        f"- CURRENT_P2: `{P2_PATH}`; SHA256 `{p2_hash}`.",
        "- Scientific reanalysis: no.",
        "- Numerical result change: no.",
        "- Existing reported numerical estimate change: no.",
        "- Retained sensitivity outcomes added to Results: yes.",
        "- Primary conclusion change: no.",
        "- Unresolved P0: none.",
        "- Unresolved P1: none.",
        "- Unresolved P2: none.",
        "- Deferred P2: none.",
        "- Numeric gate: PENDING.",
        "- Forbidden-language gate: PENDING.",
        "- Two-pass deterministic hash proof: PENDING.",
        "- Final status: PENDING_GATE.",
        "",
        "The cross-freeze revision corrected semantic numeric bindings, recovered the TCGA anchor rule, specified supporting-layer aggregation/missingness, and reported retained sensitivity outcomes. The executable numeric, semantic, forbidden-language, and deterministic-hash gates determine the final status.",
        "",
    ]
    (root / "current_p2_closure/CURRENT_P2_CLOSURE.md").write_text("\n".join(md), encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = Path(args.analysis_root).resolve()
    if sha256_file(root / SOURCE_P1) != SOURCE_P1_SHA256:
        raise RuntimeError("Immutable CURRENT_P1 hash mismatch before generation")
    snapshots = copy_evidence_inputs(root)
    capture_environment(root)
    export_tcga_selection(root)
    tcga_anchor_rule = export_tcga_anchor_selection(root, snapshots)
    write_rules(root, snapshots, tcga_anchor_rule)
    facts = compute_facts(root)
    write_json(root / DERIVATIONS, {"project": PROJECT, "facts": facts})
    manuscript, numeric_rows = build_manuscript(root, facts)
    (root / P2_PATH).write_text(manuscript, encoding="utf-8")
    write_numeric_map(root, numeric_rows)
    records = adjudication_records(root, snapshots)
    write_adjudication(root, records)
    write_commands(root)
    write_cross_freeze_revision(root)
    write_initial_closure(root)
    print(json.dumps({"project": PROJECT, "status": "BUILT_PENDING_GATE", "p2_sha256": sha256_file(root / P2_PATH), "numeric_occurrences": len(numeric_rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
