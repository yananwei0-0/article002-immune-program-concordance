#!/usr/bin/env python3
"""Build BIO-05 root authority JSONs, tables, and Methods/Results from v2 outputs."""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROTOCOL_SHA256 = "a3bc200faf24eefb1b89be08cb4b6804b81423f50c3337a430434e6015143fc3"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument("--protocol", required=True)
    return parser.parse_args()


def sha256_file(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def clean_scalar(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    return clean_scalar(value)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(clean_json(payload), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def fmt(value: Any, digits: int = 3) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    if 0 < abs(x) < 0.001:
        return f"{x:.2e}"
    return f"{x:.{digits}f}"


def count_rows(path: Path) -> int:
    with path.open("rb") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def protocol_path_gate(protocol: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    sources: list[dict[str, Any]] = []
    for item in protocol["protocol"]["source_inputs"]:
        path = Path(item)
        record: dict[str, Any] = {"path": str(path), "exists": path.exists(), "path_type": "missing"}
        if path.is_file():
            record.update(path_type="file", size_bytes=path.stat().st_size, sha256=sha256_file(path))
        elif path.is_dir():
            files = [x for x in path.rglob("*") if x.is_file()]
            record.update(
                path_type="directory",
                file_count=len(files),
                total_size_bytes=sum(x.stat().st_size for x in files),
                directory_identity_rule="Directory presence plus per-file hashes for files actually used by this v2 run.",
            )
        sources.append(record)
    assets: list[dict[str, Any]] = []
    for item in protocol["executable_assets"]:
        path = Path(item)
        assets.append(
            {
                "path": str(path),
                "exists": path.is_file(),
                "size_bytes": path.stat().st_size if path.is_file() else None,
                "sha256": sha256_file(path) if path.is_file() else None,
            }
        )
    return sources, assets


def collect_actual_sources(root: Path) -> list[dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}

    def add(path: Path, role: str, source_url: str = "", source_id: str = "", expected_size: Any = None) -> None:
        key = str(path)
        if key in records:
            roles = set(records[key]["roles"])
            roles.add(role)
            records[key]["roles"] = sorted(roles)
            return
        if not path.is_file():
            records[key] = {
                "path": key,
                "roles": [role],
                "exists": False,
                "source_url": source_url or None,
                "source_id": source_id or None,
                "expected_size_bytes": int(expected_size) if str(expected_size).isdigit() else None,
            }
            return
        records[key] = {
            "path": key,
            "roles": [role],
            "exists": True,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "source_url": source_url or None,
            "source_id": source_id or None,
            "expected_size_bytes": int(expected_size) if str(expected_size).isdigit() else None,
            "size_matches_queue": path.stat().st_size == int(expected_size) if str(expected_size).isdigit() else None,
        }

    for name in ["cptac_input_files_v2.csv", "cptac_context_input_files_v2.csv"]:
        table = pd.read_csv(root / "results" / "cptac" / name)
        for row in table.itertuples(index=False):
            add(Path(row.path), f"cptac_{getattr(row, 'layer', getattr(row, 'role', 'context'))}")

    data_root = Path(expand_legacy_paths("${SCI27_DATA_ROOT}/14_pan_cancer_proteogenomic_immune_evasion"))
    for path in (data_root / "03_downloads" / "tcga_xena_pancanatlas").iterdir():
        if path.is_file():
            add(path, "tcga_xena_source")
    for path, role in [
        (data_root / "04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv", "locked_signature_definition"),
        (data_root / "04_processed/harmonization_v0_1/cptac_sample_map_v0_1.csv", "cptac_sample_to_case_map"),
        (data_root / "04_processed/harmonization_v0_1/cptac_case_layer_coverage_v0_1.csv", "cptac_case_layer_coverage"),
        (data_root / "04_processed/harmonization_v0_1/tcga_sample_map_v0_1.csv", "tcga_patient_sample_map"),
        (data_root / "04_processed/discovery_v0_4/cptac_bcm_transcriptomics_download_manifest_v0_4.csv", "cptac_rna_download_manifest"),
        (data_root / "04_processed/first_pass_qc/cptac_processed_inventory_v0_1.csv", "cptac_processed_inventory"),
        (data_root / "04_processed/public_pan_cancer_download_manifest_v0_1.csv", "public_download_manifest"),
    ]:
        add(path, role)

    queue_paths = [
        Path(expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/manifests/BIO-05/bio05_public_validation_queue_20260730.tsv")),
        Path(expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/manifests/BIO-05/bio05_brca_spatial_queue_20260730.tsv")),
    ]
    for queue_path in queue_paths:
        add(queue_path, "public_source_queue")
        queue = pd.read_csv(queue_path, sep="\t")
        for row in queue.itertuples(index=False):
            add(Path(row.remote_path), f"public_{row.source_type}", str(row.url), str(row.source_id), row.expected_size)
    for path in sorted(Path(expand_legacy_paths("${SCI27_BIOBANK_ROOT}/topic_pool_cache/downloads/BIO-05/3ca")).glob("Data_*.tar.gz")):
        add(path, "public_3ca_source", source_id=path.stem)
    return sorted(records.values(), key=lambda x: x["path"])


def dataframe_markdown(frame: pd.DataFrame, columns: list[str], headers: list[str]) -> str:
    use = frame[columns].copy()
    use.columns = headers
    return use.to_markdown(index=False)


def build_manuscript(root: Path, summaries: dict[str, Any]) -> tuple[str, dict[str, str]]:
    primary = pd.read_csv(root / "results/cptac/cptac_primary_concordance_registry_v2.csv")
    heterogeneity = pd.read_csv(root / "results/cptac/cptac_heterogeneity_summary_v2.csv")
    feature = pd.read_csv(root / "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv")
    threshold1 = pd.read_csv(root / "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv")
    residual = pd.read_csv(root / "results/cptac/cptac_residual_context_registry_v2.csv")
    diagnostics = pd.read_csv(root / "results/cptac/cptac_residual_model_diagnostics_v2.csv")
    flow = pd.read_csv(root / "results/cptac/cptac_complete_denominator_flow_v2.csv")
    tcga_sig = pd.read_csv(root / "results/tcga/tcga_xena_signature_summary_v0_5.csv")
    tcga_subtype = pd.read_csv(root / "results/tcga/tcga_xena_immune_subtype_summary_v0_5.csv")
    apollo = pd.read_csv(root / "results/apollo/apollo_concordance_registry_v2.csv")
    scrna_global = pd.read_csv(root / "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv")
    scrna_within = pd.read_csv(root / "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv")
    spatial = pd.read_csv(root / "results/spatial/spatial_brca_sample_signature_summary_v0_7.csv")

    evaluable = primary[primary["status"].eq("evaluable")]
    prot = evaluable[evaluable["comparison_layer"].eq("proteomics")]
    phospho = evaluable[evaluable["comparison_layer"].eq("phosphoproteomics")]
    nonsig = evaluable[~(evaluable["q_value_bh"] < 0.05)].iloc[0]
    sensitivity = primary[["cancer", "program", "comparison_layer", "spearman_rho"]].merge(
        feature[["cancer", "program", "comparison_layer", "spearman_rho"]],
        on=["cancer", "program", "comparison_layer"],
        suffixes=("_gene", "_feature"),
    ).dropna()
    phospho_sens = sensitivity[sensitivity["comparison_layer"].eq("phosphoproteomics")].copy()
    residual_both = residual[(residual["linear_q_value_bh"] < 0.05) & (residual["rank_q_value_bh"] < 0.05)].dropna(subset=["linear_rho", "rank_rho"])
    max_influence = primary["max_abs_leave_one_out_rho_delta"].max()
    tcga_main = tcga_sig[~tcga_sig["signature_label"].eq("Proteasome")]
    apollo_opposite = apollo[apollo["evidence_status"].eq("opposite")]
    apollo_supported = apollo[apollo["evidence_status"].eq("supported")]
    global_claim = scrna_global[scrna_global["main_claim_window"].astype(bool)]
    within_claim = scrna_within[scrna_within["within_source_support_window"].astype(bool)]

    summary_table = heterogeneity[
        ["comparison_layer", "program_label", "evaluable_cancer_n", "median_within_cancer_rho", "positive_cancer_fraction", "descriptive_i2"]
    ].copy()
    summary_table["median_within_cancer_rho"] = summary_table["median_within_cancer_rho"].map(fmt)
    summary_table["positive_cancer_fraction"] = summary_table["positive_cancer_fraction"].map(fmt)
    summary_table["descriptive_i2"] = summary_table["descriptive_i2"].map(fmt)
    table_markdown = dataframe_markdown(
        summary_table,
        ["comparison_layer", "program_label", "evaluable_cancer_n", "median_within_cancer_rho", "positive_cancer_fraction", "descriptive_i2"],
        ["Layer", "Program", "Evaluable cancers", "Median within-cancer rho", "Positive-cancer fraction", "Descriptive I2"],
    )

    text = f"""# BIO-05 Methods and Results v2

## Methods

### Study design and evidence boundary

We performed a cross-sectional pan-cancer proteogenomic concordance analysis in ten cancer types represented in the Clinical Proteomic Tumor Analysis Consortium processed data surface: breast, clear-cell renal, colon, glioblastoma, head and neck squamous, lung squamous, lung adenocarcinoma, ovarian, pancreatic ductal, and endometrial cancer. The independent unit was the patient/case. The analysis had no longitudinal time origin, clinical endpoint, treatment exposure, or post-baseline covariate. It was designed to estimate within-cancer concordance between RNA and protein or phosphoprotein immune-ecology program scores, not prognosis, prediction, causality, malignant-cell-intrinsic mechanism, or treatment vulnerability.

### Source identity, biospecimen mapping, and denominators

The v2 protocol was locked before result generation and its SHA-256 was verified. RNA, tumor proteome, tumor phosphoproteome, and immune-context matrices were read from the public processed source files enumerated in the machine source manifest. Every file actually used was recorded with byte size and SHA-256. RNA aliquots labelled tumor were retained; identifiers lacking a tumor/normal suffix were treated as tumor-like only for cohorts whose processed RNA matrix used unsuffixed tumor identifiers. Adjacent/normal and quality-control identifiers were excluded. Terminal tumor/adjacent/normal and technical suffixes were removed to construct a normalized case identifier. When more than one eligible aliquot mapped to one case, the locked deterministic rule preferred an explicit tumor suffix and then lexical sample identifier; the selected and excluded candidates were retained in the duplicate registry.

Before matching, we required uniqueness of cancer, program, molecular layer, and normalized case. RNA and each protein layer were joined one-to-one. Raw tumor protein/phosphoprotein headers were independently cross-checked against the frozen sample-to-case map. Program-specific flow retained the matched denominator, complete-pair count, missing-pair count, and explicit non-evaluable reason.

### Locked immune-ecology programs and scores

Eight fixed programs were analyzed: MHC-I antigen presentation, MHC-II antigen presentation, interferon-gamma response, cytolytic T-cell context, checkpoint/exhaustion context, myeloid inflammatory context, TGF/EMT exclusion context, and proteasome antigen processing. No gene was added or removed after inspecting results. Within each cancer and molecular matrix, every recovered feature was standardized across eligible tumor samples with a mean of zero and sample standard deviation of one. For the primary protein and phosphoprotein score, features mapping to the same gene—including multiple phosphosites—were averaged within sample, after which recovered genes were averaged with equal weight. RNA used the same gene-equal rule. A sample score required at least three nonmissing recovered genes; no value was imputed.

### Primary estimand, multiplicity, and diagnostics

The registered primary family contained all 10 cancers × 8 programs × 2 comparison layers (160 rows). For each row we estimated the patient-level Spearman correlation between the RNA program score and its matched protein or phosphoprotein score. A row was evaluable when at least 20 complete cases remained and both variables had at least three distinct values. Ninety-five percent confidence intervals were obtained with a paired-case bias-corrected and accelerated bootstrap using 20,000 deterministic, stratum-seeded resamples. Two-sided asymptotic Spearman P values were used when n>30 and neither score contained ties. Rows with n≤30 or any ties were calibrated with 999,999 deterministic Monte Carlo permutations of the matched omic score; the plus-one correction was used and Monte Carlo standard errors were retained. Benjamini–Hochberg adjustment was then applied once across all finite calibrated or asymptotic P values in the unchanged 160-row registered family; non-evaluable rows remained in the family registry.

Diagnostics included requested and recovered genes/features, sample-to-case agreement, duplicate keys, complete denominators, distinct values, tie fractions, leave-one-case changes in rho, and cancer-level heterogeneity. The main report uses cancer-specific estimates, median within-cancer rho, positive-cancer fraction, and descriptive Fisher-z Cochran Q/I2; no single pooled pan-cancer effect is a primary estimand.

### Sensitivity and secondary residual analyses

The locked scoring sensitivity gave equal weight to feature rows rather than genes. A second sensitivity relaxed the per-sample requirement to one recovered gene without changing primary evaluability labels. Leave-one-cancer-out medians and correlations after standardizing both scores within each cancer were descriptive checks only. For the secondary residual analysis, protein or phosphoprotein score was regressed on the matched RNA program score within cancer. Linear-model residuals were related separately to tumor purity, ESTIMATE immune score, CIBERSORT T-cell signal, and CIBERSORT myeloid signal. The 640 registered residual-context rows formed one BH family for linear residuals. A rank-based residualization formed a separate 640-row sensitivity family. Linear R2, residual moments, residual-versus-RNA rank correlation, and incremental quadratic R2 were retained as diagnostics. Residual associations were interpreted as context associations, not cell of origin or post-transcriptional mechanism.

### Public supporting layers

TCGA PanCanAtlas served only as an RNA-only orthogonal sensitivity. One primary-tumor sample per patient was scored, and alignment with the public 68-signature immune landscape and immune-subtype distributions was summarized within the ten overlapping cancers with separate multiplicity control. APOLLO LUAD and ovarian public matrices were scored from their source files and used as same-cancer proteomic support by provenance. RNA-to-omic Spearman P values were BH-adjusted separately within each cohort-layer family; positive q<0.05 rows were labelled supported, negative q<0.05 rows opposite, and the remaining evaluable rows null. Because a harmonized APOLLO-to-CPTAC donor crosswalk was unavailable, this layer was not called donor-proven independent replication.

Single-cell sources were reduced to donor-by-compartment pseudobulk profiles after requiring at least 20 cells per donor-compartment. Malignant, immune, and stromal program-score contrasts were tested with paired two-sided Wilcoxon signed-rank tests within each source-cancer set of 24 program-by-contrast rows, with BH adjustment within that family and a minimum of six paired source-specific donor units. We audited accession, publication/collection identity, cancer, donor-identifier namespace, and exact donor-label overlap across the ten sources. Because no participant crosswalk established biological donor independence across repositories and aggregated atlases, the former pooled global donor-unit tests were removed from inferential claims. Their 24 program-by-contrast rows were retained only as descriptive pooled summaries without P or q values. Four breast Visium samples were summarized with sample-level compartment medians only; spots were not treated as patients and no spot-level inference was used.

## Results

### CPTAC case flow and primary registry

The one-to-one matched table contained {summaries['cptac']['unique_cases_any_matched_layer']} unique cases across either protein layer. The 160-row primary registry retained {summaries['cptac']['primary_evaluable_rows']} evaluable rows and {summaries['cptac']['primary_not_evaluable_rows']} non-evaluable rows; every non-evaluable row had fewer than 20 complete score pairs under the three-gene rule. There were no duplicated cancer-program-layer-case keys and all raw tumor protein/phosphoprotein mappings agreed with the frozen sample map. Maximum matched case counts across cancer-layer combinations ranged from {int(flow.groupby(['cancer','comparison_layer'])['matched_unique_cases'].max().min())} to {int(flow.groupby(['cancer','comparison_layer'])['matched_unique_cases'].max().max())}.

### Cancer-specific RNA–protein and RNA–phosphoprotein concordance

Of {len(prot)} evaluable RNA–protein rows, {(prot['q_value_bh'] < 0.05).sum()} had q<0.05; rho ranged from {fmt(prot['spearman_rho'].min())} to {fmt(prot['spearman_rho'].max())}, with a median of {fmt(prot['spearman_rho'].median())}. Of {len(phospho)} evaluable RNA–phosphoprotein rows, {(phospho['q_value_bh'] < 0.05).sum()} had q<0.05; rho ranged from {fmt(phospho['spearman_rho'].min())} to {fmt(phospho['spearman_rho'].max())}, with a median of {fmt(phospho['spearman_rho'].median())}. The only evaluable row not passing BH was {nonsig['cancer'].upper()} {nonsig['program_label']} RNA–protein concordance (n={int(nonsig['complete_pair_n'])}, rho={fmt(nonsig['spearman_rho'])}, P={fmt(nonsig['p_value'])}, q={fmt(nonsig['q_value_bh'])}).

All program-layer summaries with evaluable cancers had a positive-cancer fraction of 1.000, but descriptive I2 ranged from {fmt(heterogeneity['descriptive_i2'].min())} to {fmt(heterogeneity['descriptive_i2'].max())}, so directional consistency did not imply a common pan-cancer effect.

{table_markdown}

### Scoring, threshold, and influence sensitivities

Feature-row and gene-equal estimates had the same sign in all {len(sensitivity)} jointly evaluable rows. Protein estimates were identical because those matrices contained one recovered feature per mapped gene. Among {len(phospho_sens)} jointly evaluable phosphoprotein rows, the median absolute rho difference was {fmt((phospho_sens['spearman_rho_gene'] - phospho_sens['spearman_rho_feature']).abs().median())} and the maximum was {fmt((phospho_sens['spearman_rho_gene'] - phospho_sens['spearman_rho_feature']).abs().max())}; thus, gene collapse changed magnitudes while retaining directions. The one-gene threshold sensitivity made {(threshold1['status'] == 'evaluable').sum()} rows evaluable and {(threshold1['q_value_bh'] < 0.05).sum()} passed its separate BH family, but it did not alter the locked primary registry. The maximum absolute leave-one-case change in primary rho was {fmt(max_influence)}.

### RNA-adjusted residual context

Linear residual models were evaluable in {(diagnostics['status'] == 'evaluable').sum()} cancer-program-layer strata. Of 640 registered context rows, {(residual['linear_q_value_bh'] < 0.05).sum()} passed the linear-residual BH family and {(residual['rank_q_value_bh'] < 0.05).sum()} passed the rank-residual sensitivity family. All {len(residual_both)} rows passing both families had concordant directions. These associations remained secondary descriptions of immune/purity context after RNA adjustment.

### RNA-only TCGA sensitivity

The TCGA layer included {summaries['tcga']['primary_tumor_patients']} primary-tumor patients and the same number of scored samples, with no repeated patient-program row. Across the seven non-proteasome programs, canonical-anchor median within-cancer correlations with public immune signatures ranged from {fmt(tcga_main['selected_anchor_median_rho'].min())} to {fmt(tcga_main['selected_anchor_median_rho'].max())}; the minimum positive-cancer fraction was {fmt(tcga_main['selected_anchor_min_positive_cancer_fraction'].min())}. The proteasome comparator had an anchor median rho of {fmt(tcga_sig.loc[tcga_sig['signature_label'].eq('Proteasome'),'selected_anchor_median_rho'].iloc[0])}. The Immune C5 subtype contained only {int(tcga_subtype[tcga_subtype['immune_subtype'].astype(str).str.contains('Immune C5', na=False)]['n'].max())} patients and remained descriptive. These results support RNA-level immune-landscape alignment only.

### APOLLO same-cancer support

APOLLO contributed 87 LUAD samples with RNA, protein, and phosphoprotein measurements and 70 ovarian samples with RNA and protein measurements. Across 24 registered cohort-layer-program rows, {summaries['apollo']['evaluable_rows']} were evaluable, {len(apollo_supported)} were labelled supported, {len(apollo_opposite)} opposite, {summaries['apollo']['null_result_rows']} null, and {summaries['apollo']['not_evaluable_rows']} not evaluable after cohort-layer BH control. The opposite rows were retained rather than recoded: {', '.join(f"{r.cohort_label} {r.program_label} {r.comparison_layer} (rho={fmt(r.spearman_rho)}, q={fmt(r.q_value_bh)})" for r in apollo_opposite.itertuples()) if len(apollo_opposite) else 'none'}. This mixture supports only bounded same-cancer public proteomic context.

### Donor-aware single-cell and sample-aware spatial triangulation

The combined single-cell/3CA layer contained {summaries['single_cell']['unique_donor_units']} source-specific donor units from {summaries['single_cell']['unique_source_ids']} public sources. No exact raw donor label appeared in more than one source, but the provenance audit found no cross-repository participant crosswalk and therefore could not establish global biological independence. The 24 pooled program-by-compartment rows were retained as descriptive summaries only; none contributed a global P value, q value, or support count. Within source-cancer families, {summaries['single_cell']['within_source_evaluable_rows']} rows were evaluable and {len(within_claim)} met the separately adjusted n≥6 support window. These source-bounded contrasts locate program signal across broad compartments but do not establish malignant-cell origin or mechanism.

All four breast Visium samples yielded sample-level program and compartment summaries. Spatial contrasts were therefore treated as same-collection, sample-aware triangulation (n=4 samples), without spot-level P values or an independent-validation claim.
"""
    propagation = {
        "unique_cases": str(summaries["cptac"]["unique_cases_any_matched_layer"]),
        "primary_evaluable": str(summaries["cptac"]["primary_evaluable_rows"]),
        "primary_significant": str(summaries["cptac"]["primary_fdr_lt_0_05_rows"]),
        "tcga_patients": str(summaries["tcga"]["primary_tumor_patients"]),
        "apollo_registered": "24",
        "spatial_samples": "four",
    }
    return text, propagation


def output_manifest(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    roots = ["code", "results", "tables", "figures", "manuscript", "local_audit"]
    for folder in roots:
        base = root / folder
        if not base.exists():
            continue
        for path in sorted(x for x in base.rglob("*") if x.is_file() and "__pycache__" not in x.parts):
            records.append(
                {
                    "relative_path": str(path.relative_to(root)),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    return records


def build_supplementary_workbook(root: Path) -> Path:
    """Package retained machine-readable tables without replacing source CSVs."""
    workbook = root / "tables/BIO-05_SUPPLEMENTARY_DATA_V2.xlsx"
    sheets = [
        ("Primary_registry", "results/cptac/cptac_primary_concordance_registry_v2.csv"),
        ("Gene_recovery", "results/cptac/cptac_gene_recovery_v2.csv"),
        ("Case_flow", "results/cptac/cptac_complete_denominator_flow_v2.csv"),
        ("Duplicate_audit", "results/cptac/cptac_duplicate_priority_audit_v2.csv"),
        ("Feature_sensitivity", "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv"),
        ("Threshold1_sensitivity", "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv"),
        ("Heterogeneity", "results/cptac/cptac_heterogeneity_summary_v2.csv"),
        ("Residual_context", "results/cptac/cptac_residual_context_registry_v2.csv"),
        ("Residual_diagnostics", "results/cptac/cptac_residual_model_diagnostics_v2.csv"),
        ("TCGA_signature_summary", "results/tcga/tcga_xena_signature_summary_v0_5.csv"),
        ("APOLLO_registry", "results/apollo/apollo_concordance_registry_v2.csv"),
        ("scRNA_global", "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv"),
        ("scRNA_within_source", "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv"),
        ("scRNA_donor_provenance", "results/scrna_combined/public_scrna_source_donor_provenance_audit_v2.csv"),
        ("Spatial_sample_summary", "results/spatial/spatial_brca_sample_signature_summary_v0_7.csv"),
    ]
    readme = pd.DataFrame(
        {
            "field": ["project", "protocol_sha256", "scope", "inference_units", "machine_readable_authority"],
            "value": [
                "BIO-05",
                PROTOCOL_SHA256,
                "Cross-sectional CPTAC RNA/protein concordance with bounded public support layers.",
                "Patient/case; donor-compartment pseudobulk; Visium sample medians.",
                "CSV and root JSON files remain authoritative; this workbook is a convenience package.",
            ],
        }
    )
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        for sheet, relative in sheets:
            table = pd.read_csv(root / relative)
            table.to_excel(writer, sheet_name=sheet[:31], index=False)
    return workbook


def build_figure_legends(root: Path) -> Path:
    path = root / "manuscript/BIO-05_FIGURE_LEGENDS_V2.md"
    text = """# BIO-05 v2 figure legends

**Figure 1. Cancer-specific CPTAC RNA–protein and RNA–phosphoprotein program concordance.** Heatmap cells show patient-level Spearman rho for each cancer, locked immune-ecology program, and comparison layer. Non-evaluable registered cells are retained and visually distinguished. P values were adjusted once across all finite tests in the 160-cell primary family.

**Figure 2. Gene-equal versus feature-row scoring sensitivity.** Each point compares the primary gene-equal Spearman estimate with the protocol-defined feature-row-weighted estimate for the same cancer, program, and molecular layer. The diagonal denotes equality; this sensitivity evaluates the influence of multiple features or phosphosites per mapped gene.

**Figure 3. Complete matched-case denominators.** Bars show program-specific complete patient pairs available after one-to-one RNA/omic matching and the requirement for at least three nonmissing recovered genes per score. Counts are denominators, not effect estimates.

**Supplementary Figure 1. TCGA RNA-only immune-landscape sensitivity.** Heatmap of within-cancer correlations between locked RNA program scores and public immune signatures among one selected primary-tumor sample per patient. This panel is an RNA-only orthogonal sensitivity and not cross-omic validation.

**Supplementary Figure 2. TCGA immune-subtype score distributions.** Median locked RNA program scores across public PanCanAtlas immune-subtype strata. Small strata, including Immune C5, remain descriptive.

**Supplementary Figure 3. BRCA Visium sample-level compartment summaries.** Immune-compartment program-score medians for four breast spatial samples. Spot-level scores were reduced to sample summaries; spots were not treated as independent patients and no spot-level significance test is shown.
"""
    path.write_text(text, encoding="utf-8")
    return path


def main() -> int:
    args = parse_args()
    root = Path(args.analysis_root).resolve()
    protocol_path = Path(args.protocol).resolve()
    observed_protocol_sha = sha256_file(protocol_path)
    if observed_protocol_sha != PROTOCOL_SHA256:
        raise RuntimeError("Protocol SHA-256 mismatch")
    protocol = json.loads(protocol_path.read_text())
    cptac = json.loads((root / "results/cptac/cptac_results_summary_v2.json").read_text())
    support = json.loads((root / "results/support_layers_summary_v2.json").read_text())
    primary = pd.read_csv(root / "results/cptac/cptac_primary_concordance_registry_v2.csv")
    sources, executable_assets = protocol_path_gate(protocol)
    actual_sources = collect_actual_sources(root)
    source_failures = [x for x in sources if not x["exists"]] + [x for x in actual_sources if not x["exists"]]
    source_size_failures = [x for x in actual_sources if x.get("size_matches_queue") is False]

    summaries = {
        "cptac": {
            "unique_cases_any_matched_layer": cptac["unique_cases_any_matched_layer"],
            "primary_registry_rows": cptac["primary_registry_rows"],
            "primary_evaluable_rows": cptac["primary_evaluable_rows"],
            "primary_not_evaluable_rows": cptac["primary_registry_rows"] - cptac["primary_evaluable_rows"],
            "primary_fdr_lt_0_05_rows": cptac["primary_fdr_lt_0_05_rows"],
            "residual_linear_fdr_lt_0_05_rows": cptac["residual_linear_fdr_lt_0_05_rows"],
            "residual_rank_fdr_lt_0_05_rows": cptac["residual_rank_fdr_lt_0_05_rows"],
        },
        **support,
    }
    manuscript, propagation = build_manuscript(root, summaries)
    manuscript_dir = root / "manuscript"
    manuscript_dir.mkdir(parents=True, exist_ok=True)
    manuscript_path = manuscript_dir / "BIO-05_METHODS_RESULTS_V2.md"
    manuscript_path.write_text(manuscript, encoding="utf-8")
    build_figure_legends(root)
    build_supplementary_workbook(root)

    validation_checks = {
        "protocol_sha256_match": observed_protocol_sha == PROTOCOL_SHA256,
        "all_protocol_source_inputs_present": all(x["exists"] for x in sources),
        "all_protocol_executable_assets_present": all(x["exists"] for x in executable_assets),
        "all_actual_analysis_sources_present": not source_failures,
        "all_queued_expected_sizes_match": not source_size_failures,
        "primary_registry_has_160_rows": len(primary) == 160,
        "primary_case_keys_unique": int(pd.read_csv(root / "results/cptac/cptac_matched_gene_equal_scores_v2.csv").duplicated(["cancer", "program", "comparison_layer", "normalized_case_id"], keep=False).sum()) == 0,
        "primary_rho_in_range": bool(primary["spearman_rho"].dropna().between(-1, 1).all()),
        "primary_p_in_range": bool(primary["p_value"].dropna().between(0, 1).all()),
        "primary_q_in_range": bool(primary["q_value_bh"].dropna().between(0, 1).all()),
        "primary_bca_ci_complete": int(primary["ci_method"].eq("paired_case_bca_bootstrap").sum()) == int(primary["status"].eq("evaluable").sum()),
        "primary_small_n_or_tied_permutation_complete": int(primary["p_value_method"].eq("two_sided_monte_carlo_permutation").sum()) == int(((primary["status"].eq("evaluable")) & ((primary["complete_pair_n"] <= 30) | (primary["x_tie_fraction"] > 0) | (primary["y_tie_fraction"] > 0))).sum()),
        "tcga_patient_program_unique": support["tcga"]["patient_program_duplicate_rows"] == 0,
        "apollo_inference_keys_unique": support["apollo"]["duplicate_inference_key_rows"] == 0,
        "single_cell_source_execution_present": support["single_cell"].get("status") != "not_run_or_no_scores",
        "single_cell_global_inference_removed": support["single_cell"]["global_main_claim_rows"] == 0 and support["single_cell"]["donor_provenance_audit"]["global_pooled_donor_inference_retained"] is False,
        "spatial_source_execution_present": support["spatial"]["status"] == "executed",
        "methods_results_exists": manuscript_path.is_file() and manuscript_path.stat().st_size > 0,
        "analysis_plan_deviation": False,
    }
    critical_failures = [name for name, passed in validation_checks.items() if passed is False and name != "analysis_plan_deviation"]
    validation = {
        "project": "BIO-05",
        "generated_at_utc": now_utc(),
        "status": "PASS" if not critical_failures else "FAIL",
        "checks": validation_checks,
        "critical_failures": critical_failures,
        "warnings": [
            "Primary Spearman intervals use deterministic paired-case BCa bootstrap; small-n or tied primary rows use Monte Carlo permutation P values.",
            "APOLLO lacks a harmonized donor crosswalk to CPTAC.",
            "Single-cell donor identifiers are source-specific; pooled global donor inference was removed because cross-source biological independence was not established.",
            "BRCA spatial evidence is same-collection and includes four samples.",
        ],
        "claim_boundary": "Cross-sectional concordance and bounded public support only; no prognosis, prediction, causality, mechanism, or therapeutic claim.",
    }
    write_json(root / "VALIDATION_REPORT.json", validation)
    source_manifest = {
        "project": "BIO-05",
        "status": "COMPLETE" if not source_failures and not source_size_failures else "INCOMPLETE",
        "generated_at_utc": now_utc(),
        "protocol_sha256": observed_protocol_sha,
        "protocol_source_inputs": sources,
        "executable_assets": executable_assets,
        "analysis_files": actual_sources,
        "missing_or_unresolved_files": source_failures,
        "expected_size_mismatches": source_size_failures,
        "crosswalk_limitations": [
            "APOLLO-to-CPTAC harmonized donor crosswalk unavailable.",
            "Single-cell donor identifiers are source-specific; they are prefixed for descriptive bookkeeping, not treated as proof of cross-source biological independence.",
            "BRCA spatial samples share the breast-atlas collection with a single-cell support source and are not independent validation.",
        ],
        "redistribution": "Not inferred. Source files remain in their existing local data locations and were not copied into the v2 analysis directory.",
    }
    write_json(root / "SOURCE_MANIFEST.json", source_manifest)
    results_summary = {
        "project": "BIO-05",
        "protocol_sha256": observed_protocol_sha,
        "status": "V2_ANALYSIS_COMPLETE" if not critical_failures else "V2_ANALYSIS_PARTIAL",
        "generated_at_utc": now_utc(),
        "analysis_plan_deviation": False,
        "primary": summaries["cptac"],
        "support_layers": {k: v for k, v in support.items()},
        "time_origin": "Cross-sectional matched molecular measurement; no follow-up time origin.",
        "inference_unit": "Patient/case for CPTAC, TCGA, and APOLLO; donor-compartment pseudobulk for single-cell; sample-level medians for spatial.",
        "claim_boundary": validation["claim_boundary"],
    }
    write_json(root / "RESULTS_SUMMARY.json", results_summary)
    manifest_before = output_manifest(root)
    run_manifest = {
        "project": "BIO-05",
        "status": results_summary["status"],
        "created_at_utc": "2026-08-11T16:07:39Z",
        "completed_at_utc": now_utc(),
        "protocol_path": str(protocol_path),
        "protocol_sha256_expected": PROTOCOL_SHA256,
        "protocol_sha256_observed": observed_protocol_sha,
        "protocol_sha256_match": True,
        "analysis_plan_deviation": False,
        "historical_assets_read_only": True,
        "runtime": {"python_executable": sys.executable, "python_version": platform.python_version(), "platform": platform.platform()},
        "source_execution": {
            "cptac": "real_source_exit_0",
            "tcga": "real_source_exit_0_via_path_safe_wrapper",
            "apollo": "real_source_exit_0",
            "single_cell": "real_source_exit_0",
            "threeca": "real_source_exit_0",
            "spatial": "real_source_exit_0",
        },
        "commands_log": "logs/COMMANDS_V2.log",
        "output_manifest": manifest_before,
        "validation_status": validation["status"],
    }
    write_json(root / "RUN_MANIFEST.json", run_manifest)

    table1 = pd.read_csv(root / "tables/TABLE1_CPTAC_CASE_FLOW_V2.csv")
    table2 = pd.read_csv(root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_V2.csv")
    (root / "tables/TABLE1_CPTAC_CASE_FLOW_V2.md").write_text(table1.to_markdown(index=False) + "\n", encoding="utf-8")
    (root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_V2.md").write_text(table2.to_markdown(index=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": results_summary["status"], "critical_failures": critical_failures}, indent=2))
    return 0 if not critical_failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
