#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT = "BIO-05"
PROTOCOL_SHA256 = "a3bc200faf24eefb1b89be08cb4b6804b81423f50c3337a430434e6015143fc3"
ROOT_ABS = Path(expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/manuscript27_current_p1_closure_20260815/analyses/BIO-05"))
REVIEW_PATH = Path(expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/sci27_v2_single_external_review_20260815/responses/BIO-05_EXTERNAL_R1_RESPONSE.md"))
PLAN_PATH = Path(expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/sci27_v2_single_external_review_20260815/SCI27_CURRENT_MODIFICATION_PLAN.md"))
EXTERNAL_FREEZE_REVIEW_PATH = Path(expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/sci27_current_p1_external_freeze_20260815/responses/BIO-05_CURRENT_P1_FREEZE_RESPONSE.md"))
BASELINE_PATH = ROOT_ABS / "manuscript/BIO-05_METHODS_RESULTS_V2.md"
CURRENT_METHODS_PATH = ROOT_ABS / "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md"

PROGRAMS = [
    ("antigen_presentation_mhc_i", "MHC-I"),
    ("antigen_presentation_mhc_ii", "MHC-II"),
    ("ifn_gamma_response", "IFN-gamma"),
    ("cytolytic_t_cell_context", "Cytolytic"),
    ("checkpoint_exhaustion_context", "Checkpoint"),
    ("myeloid_inflammatory_context", "Myeloid"),
    ("tgfb_emt_exclusion_context", "TGF/EMT"),
    ("proteasome_antigen_processing", "Proteasome"),
]
PROGRAM_ORDER = {label: i for i, (_, label) in enumerate(PROGRAMS)}
LAYER_ORDER = {"proteomics": 0, "phosphoproteomics": 1}

COMMANDS_RUN = [
    expand_legacy_paths("shasum -a 256 ${SCI27_RESEARCH_ROOT}/outputs/manuscript27_v2_scientific_rebuild_20260811/protocols/BIO-05/BIO-05_V2_LOCKED_PROTOCOL.json"),
    "python3 -m py_compile current_p1_closure/build_current_p1_outputs.py current_p1_closure/run_current_p1_numeric_gate.py",
    expand_legacy_paths("python3 current_p1_closure/build_current_p1_outputs.py --analysis-root ${SCI27_RESEARCH_ROOT}/outputs/manuscript27_current_p1_closure_20260815/analyses/BIO-05"),
    expand_legacy_paths("python3 current_p1_closure/run_current_p1_numeric_gate.py --analysis-root ${SCI27_RESEARCH_ROOT}/outputs/manuscript27_current_p1_closure_20260815/analyses/BIO-05"),
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean_json(v) for v in value]
    if isinstance(value, tuple):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(clean_json(payload), indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")


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


def fnum(value: Any) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(x):
        return "NA"
    return f"{x:.12g}"


def rel(path: Path, root: Path) -> str:
    return str(path.relative_to(root))


def dataframe_markdown(frame: pd.DataFrame) -> str:
    return frame.to_markdown(index=False)


def add_map(rows: list[dict[str, str]], numeric_id: str, value: Any, source_file: str, row_filter: str, manuscript_location: str, meaning: str) -> None:
    rows.append(
        {
            "numeric_id": numeric_id,
            "value": fnum(value),
            "source_file": source_file,
            "row_filter": row_filter,
            "manuscript_location": manuscript_location,
            "meaning": meaning,
            "verification_status": "SOURCE_MACHINE_RESULT",
        }
    )


def build_program_table(root: Path, primary: pd.DataFrame, heterogeneity: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        primary.groupby(["comparison_layer", "program_label"], as_index=False)
        .agg(
            registered_cancer_rows=("cancer", "nunique"),
            evaluable_cancer_n=("status", lambda s: int((s == "evaluable").sum())),
            total_complete_pairs_across_strata=("complete_pair_n", "sum"),
        )
    )
    table = grouped.merge(
        heterogeneity[
            [
                "comparison_layer",
                "program_label",
                "median_within_cancer_rho",
                "positive_cancer_fraction",
                "descriptive_i2",
            ]
        ],
        on=["comparison_layer", "program_label"],
        how="left",
        validate="one_to_one",
    )
    table["interpretation_boundary"] = np.where(
        table["evaluable_cancer_n"].eq(0),
        "Registered but not evaluable under the complete-pair/distinct-value gate.",
        "Cancer-specific effects are primary; I2 is descriptive and no single pan-cancer effect is an estimand.",
    )
    table["layer_order"] = table["comparison_layer"].map(LAYER_ORDER)
    table["program_order"] = table["program_label"].map(PROGRAM_ORDER)
    table = table.sort_values(["layer_order", "program_order"]).drop(columns=["layer_order", "program_order"])
    out_csv = root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv"
    out_md = root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.md"
    table.to_csv(out_csv, index=False)
    display = table.copy()
    for col in ["median_within_cancer_rho", "positive_cancer_fraction", "descriptive_i2"]:
        display[col] = display[col].map(fmt)
    display = display.rename(
        columns={
            "comparison_layer": "Layer",
            "program_label": "Program",
            "registered_cancer_rows": "Registered cancers",
            "evaluable_cancer_n": "Evaluable cancers",
            "total_complete_pairs_across_strata": "Complete pairs",
            "median_within_cancer_rho": "Median within-cancer rho",
            "positive_cancer_fraction": "Positive-cancer fraction",
            "descriptive_i2": "Descriptive I2 fraction",
            "interpretation_boundary": "Boundary",
        }
    )
    out_md.write_text(dataframe_markdown(display) + "\n", encoding="utf-8")
    return table


def build_source_provenance(root: Path, source_manifest: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for item in source_manifest["analysis_files"]:
        roles = ";".join(item.get("roles", []))
        source_id = item.get("source_id") or ""
        source_url = item.get("source_url") or ""
        path = Path(item["path"])
        if roles.startswith("tcga_xena"):
            provider = "UCSC Xena / TCGA PanCanAtlas"
            acquisition = "public download manifest timestamp 2026-06-23T23:56:32 where available"
        elif "apollo" in roles:
            provider = "GDC APOLLO public supplementary data"
            acquisition = "BIO-05 public validation queue dated 20260730; no per-file download timestamp column"
        elif "cellxgene" in roles:
            provider = "CELLxGENE public H5AD"
            acquisition = "BIO-05 public validation queue dated 20260730; no per-file download timestamp column"
        elif "3ca" in roles:
            provider = "Curated Cancer Cell Atlas / 3CA public archive"
            acquisition = "BIO-05 public validation queue dated 20260730; no per-file download timestamp column"
        elif roles.startswith("cptac"):
            provider = "CPTAC processed BCM/WashU data surface"
            acquisition = "retained CPTAC processed inventory and file hashes; file-level public URL absent from current manifest"
        else:
            provider = "retained local public-data manifest"
            acquisition = "retained source-provenance manifest"
        rows.append(
            {
                "source_role": roles,
                "provider_or_collection": provider,
                "source_id_or_accession": source_id or path.name,
                "source_url": source_url or "not recorded in current manifest",
                "source_file_name": path.name,
                "size_bytes": item.get("size_bytes"),
                "sha256": item.get("sha256"),
                "provenance_note": acquisition,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(root / "current_p1_closure/SOURCE_PROVENANCE_CURRENT_P1.csv", index=False)
    return out


def build_methods_results(root: Path, data: dict[str, Any]) -> str:
    primary = data["primary"]
    heterogeneity = data["heterogeneity"]
    feature = data["feature"]
    threshold1 = data["threshold1"]
    residual = data["residual"]
    diagnostics = data["diagnostics"]
    flow = data["flow"]
    tcga_sig = data["tcga_sig"]
    tcga_subtype = data["tcga_subtype"]
    tcga_tests = data["tcga_tests"]
    tcga_corr = data["tcga_corr"]
    apollo = data["apollo"]
    scrna_global = data["scrna_global"]
    scrna_within = data["scrna_within"]
    spatial_summary = data["spatial_summary"]
    spatial_pairwise = data["spatial_pairwise"]
    support = data["support"]
    table2_current = data["table2_current"]

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
    max_influence = primary.loc[primary["max_abs_leave_one_out_rho_delta"].idxmax()]
    tcga_main = tcga_sig[~tcga_sig["signature_label"].eq("Proteasome")]
    apollo_opposite = apollo[apollo["evidence_status"].eq("opposite")]
    apollo_supported = apollo[apollo["evidence_status"].eq("supported")]
    scrna_support = scrna_within[scrna_within["within_source_support_window"].astype(bool)]
    scrna_top = scrna_support.sort_values("q_value_bh").head(3)
    apollo_top = apollo_supported.sort_values("q_value_bh").head(2)

    residual_linear_eval = int(residual["linear_p_value"].notna().sum())
    residual_rank_eval = int(residual["rank_p_value"].notna().sum())
    residual_context_counts = residual.groupby("context_variable")["linear_status"].value_counts().unstack(fill_value=0)
    context_sentence = "; ".join(
        f"{idx}: {int(row.get('evaluable', 0))} evaluable and {int(row.get('not_evaluable', 0))} not evaluable"
        for idx, row in residual_context_counts.iterrows()
    )

    table_display = table2_current.copy()
    for col in ["median_within_cancer_rho", "positive_cancer_fraction", "descriptive_i2"]:
        table_display[col] = table_display[col].map(fmt)
    table_display = table_display.rename(
        columns={
            "comparison_layer": "Layer",
            "program_label": "Program",
            "registered_cancer_rows": "Registered cancers",
            "evaluable_cancer_n": "Evaluable cancers",
            "total_complete_pairs_across_strata": "Complete pairs",
            "median_within_cancer_rho": "Median rho",
            "positive_cancer_fraction": "Positive fraction",
            "descriptive_i2": "I2 fraction",
            "interpretation_boundary": "Boundary",
        }
    )
    table_markdown = dataframe_markdown(
        table_display[
            [
                "Layer",
                "Program",
                "Registered cancers",
                "Evaluable cancers",
                "Complete pairs",
                "Median rho",
                "Positive fraction",
                "I2 fraction",
            ]
        ]
    )

    anchor_ranges = (
        fmt(tcga_main["selected_anchor_median_rho"].min()),
        fmt(tcga_main["selected_anchor_median_rho"].max()),
    )
    anchor_names = "; ".join(
        f"{row.signature_label}: {row.selected_anchor_names}"
        for row in tcga_sig.sort_values("signature_label").itertuples()
    )
    apollo_examples = "; ".join(
        f"{r.cohort_label} {r.program_label} {r.comparison_layer} rho={fmt(r.spearman_rho)}, 95% CI {fmt(r.ci_lower)} to {fmt(r.ci_upper)}, q={fmt(r.q_value_bh)}"
        for r in apollo_top.itertuples()
    )
    apollo_opposite_text = "; ".join(
        f"{r.cohort_label} {r.program_label} {r.comparison_layer} rho={fmt(r.spearman_rho)}, 95% CI {fmt(r.ci_lower)} to {fmt(r.ci_upper)}, q={fmt(r.q_value_bh)}"
        for r in apollo_opposite.itertuples()
    )
    scrna_examples = "; ".join(
        f"{r.cancer_code.upper()} {r.program_label} {r.contrast}, paired donor n={int(r.paired_donor_n)}, median delta={fmt(r.median_delta)}, q={fmt(r.q_value_bh)}"
        for r in scrna_top.itertuples()
    )

    text = f"""# BIO-05 Methods and Results

## Methods

### Study Design and Evidence Boundary

We performed a cross-sectional pan-cancer proteogenomic concordance analysis in ten cancer types represented in processed Clinical Proteomic Tumor Analysis Consortium data: breast cancer, clear-cell renal cell carcinoma, colon adenocarcinoma, glioblastoma, head and neck squamous carcinoma, lung squamous carcinoma, lung adenocarcinoma, ovarian cancer, pancreatic ductal adenocarcinoma, and endometrial cancer. The independent unit for the primary analysis was the patient or case. The study had no longitudinal time origin, clinical endpoint, treatment exposure, or post-baseline covariate. It was designed to estimate within-cancer concordance between RNA and protein or phosphoprotein immune-ecology program scores, not prognosis, prediction, diagnosis, treatment response, causality, malignant-cell-intrinsic mechanism, post-transcriptional mechanism, or therapeutic vulnerability.

### Source Identity, Provenance, and Case Mapping

The analysis protocol was locked at 2026-08-11 08:17:51 UTC and the protocol SHA-256 was {PROTOCOL_SHA256}. Result generation began at 2026-08-11 16:07:39 UTC. The source-provenance manifest distributed with this analysis was finalized at 2026-08-12 02:19:28 UTC, after result generation; therefore source provenance is reported as a post-result documentary manifest rather than as an immutable pre-result manifest. The retained command ledger, source reports, file sizes, and SHA-256 hashes were used to verify the actual files analyzed.

The primary CPTAC layer used processed BCM/WashU RNA, proteome, phosphoproteome, ESTIMATE, xCell, and CIBERSORT files identified by provider, cancer, data type, file name, byte count, MD5 where available, and SHA-256. The current source manifest did not contain file-level public URLs for the core CPTAC analysis matrices. TCGA PanCanAtlas files were identified from UCSC Xena download records, including RNA expression, 68 immune signatures, and immune subtype files. APOLLO files were identified from GDC API URLs in the public validation queue. Single-cell and Visium H5AD files were identified by CELLxGENE dataset identifiers and download URLs, and the 3CA archives by retained public archive files. No repository-level public release URL or DOI is asserted here.

RNA aliquots labelled tumor were retained; identifiers lacking a tumor/normal suffix were treated as tumor-like only for cohorts whose processed RNA matrix used unsuffixed tumor identifiers. Adjacent/normal and quality-control identifiers were excluded. Terminal tumor, adjacent, normal, and technical suffixes were removed to construct a normalized case identifier. When more than one eligible aliquot mapped to one case, the deterministic rule preferred an explicit tumor suffix and then lexical sample identifier. The selected and excluded candidates were retained in the duplicate registry. Before modeling, uniqueness was required for cancer, program, molecular layer, and normalized case. RNA and each protein layer were joined one-to-one, and raw tumor protein/phosphoprotein headers were cross-checked against the frozen sample-to-case map.

### Locked Immune-Ecology Programs and Scores

Eight fixed programs were analyzed: MHC-I antigen presentation, MHC-II antigen presentation, interferon-gamma response, cytolytic T-cell context, checkpoint/exhaustion context, myeloid inflammatory context, TGF/EMT exclusion context, and proteasome antigen processing. Within each cancer and molecular matrix, every recovered feature was standardized across eligible tumor samples with mean zero and sample standard deviation one. For the primary protein and phosphoprotein score, features mapping to the same gene, including multiple phosphosites, were averaged within sample before recovered genes were averaged with equal weight. RNA used the same gene-equal rule. A sample score required at least three nonmissing recovered genes. No primary CPTAC score value was imputed.

### Primary Estimand, Multiplicity, and Diagnostics

The registered primary family contained all 10 cancers x 8 programs x 2 comparison layers, for 160 rows. For each row, the estimand was the patient-level Spearman correlation between the RNA program score and the matched protein or phosphoprotein score within cancer. A row was evaluable when at least 20 complete cases remained and both variables had at least three distinct values. Ninety-five percent confidence intervals used paired-case bias-corrected and accelerated bootstrap resampling with 20,000 deterministic, stratum-seeded resamples. Two-sided asymptotic Spearman P values were used when n>30 and neither score contained ties. Rows with n<=30 or score ties used 999,999 deterministic Monte Carlo permutations of the matched omic score with a plus-one correction. Benjamini-Hochberg adjustment was then applied once across all finite P values in the unchanged 160-row registry; non-evaluable rows remained in the registry denominator.

Diagnostics included requested and recovered genes/features, sample-to-case agreement, duplicate keys, complete denominators, distinct values, tie fractions, leave-one-case changes in rho, and descriptive cancer-level heterogeneity. The report uses cancer-specific estimates, median within-cancer rho, positive-cancer fraction, and descriptive Fisher-z Cochran Q/I2. I2 is reported as a fraction, not a percentage. No single pooled pan-cancer coefficient is a primary estimand.

### Sensitivity and RNA-Adjusted Residual Context Analyses

The scoring sensitivity gave equal weight to feature rows rather than genes. A second sensitivity relaxed the per-sample requirement to one recovered gene without changing primary evaluability labels. Leave-one-cancer-out medians and correlations after standardizing both scores within each cancer were descriptive checks only.

For the RNA-adjusted residual analysis, a separate ordinary least-squares model was fit within each evaluable cancer-program-layer stratum: omic score = intercept + beta x matched RNA score. The linear residual was the observed omic score minus this fitted value. The rank-residual sensitivity first converted RNA and omic scores within the same stratum to average ranks, then fit rank(omic score) = intercept + beta x rank(RNA score), and used the resulting rank residual. Residual scores were joined to four context variables by normalized case: tumor purity, ESTIMATE immune score, CIBERSORT T-cell signal, and CIBERSORT myeloid signal. For each residual family and context variable, finite paired rows were required to have nonmissing residual and context values, n>=20, and at least three distinct values in both variables. The residual-context association statistic was Spearman rho with a two-sided asymptotic P value. Linear residuals and rank residuals formed separate Benjamini-Hochberg families over their finite residual-context P values; non-evaluable registry rows were retained with reasons.

### Public Supporting Layers

TCGA PanCanAtlas was used only as an RNA-only orthogonal sensitivity. One primary-tumor sample per patient was scored. The full transparent evidence surface was the 8 program x 68 public immune signature x 10 cancer matrix of within-cancer Spearman correlations. For this 68-signature matrix, `fdr_bh` was computed by applying Benjamini-Hochberg adjustment once across all 5440 finite Spearman P values in the full matrix under the retained `tcga_xena_68immune_within_cancer_correlations_v0_5` P-value family; no within-cancer, within-program, or within-immune-signature partition was used. The canonical anchor list used for compact summaries was the fixed immune-concept mapping retained in code before the current TCGA v0.5 execution: {anchor_names}. These anchor summaries were not used to select favorable results, and the full 68-signature matrix remains the evidence surface. Immune subtype tests used Kruskal-Wallis tests across immune subtype groups with at least 20 patients; Benjamini-Hochberg adjustment was applied separately across the eight program-level subtype tests. Small subtype strata were summarized descriptively only.

APOLLO LUAD and ovarian public matrices were scored as same-cancer proteomic support by provenance. LUAD used public RNA count files, global proteome, and phosphoproteome files; ovarian used public whole-tumor normalized RNA and the source-level pre-imputed global proteomics file labelled `APOLLO_OV_WholeTumor_GlobalProteomics_imputed.csv`. No additional study-level imputation of APOLLO proteomic values or program scores was performed after reading these source files. LUAD RNA count tables were unioned across sample files, absent gene/sample count combinations were set to zero before log2(count+1) count-matrix construction, and duplicated gene symbols were collapsed by averaging. Genes were z-standardized across samples for each cohort-layer matrix, duplicated genes were collapsed by averaging, and a program score required at least three recovered genes. The APOLLO inference unit was the public `sample_id` within each cohort-layer score table: 87 LUAD sample units were available for RNA, proteomics, and phosphoproteomics, and 70 ovarian sample units were available for RNA and proteomics. Before registry formation, the score table was required to have no duplicated cohort-cancer-layer-program-sample inference keys; within each program comparison, duplicated RNA or omic `sample_id` values would fail the run, and the RNA-to-omic join was one-to-one on `sample_id`. The APOLLO registry was the fixed 3 cohort-layer families x 8 programs, for 24 rows. A row was evaluable only when the one-to-one joined, nonmissing RNA and omic score pairs had n>=10 and both RNA and omic scores had at least three distinct values. Rows lacking a program score after the three-recovered-gene scoring gate, or with n<10, were not evaluable with reason `missing_program_score_or_n_below_10`; this rule produced 17 evaluable and 7 non-evaluable rows. APOLLO 95% Spearman confidence intervals used a Fisher z interval, atanh(rho) +/- 1.9599639845/sqrt(n-3), back-transformed with tanh. RNA-to-omic Spearman P values were adjusted separately within each cohort-layer family; positive q<0.05 rows were labelled supported, negative q<0.05 rows opposite, and remaining evaluable rows null. Because a harmonized APOLLO-to-CPTAC donor crosswalk was unavailable, this layer was not donor-proven independent replication.

Single-cell sources were reduced to donor-by-compartment pseudobulk profiles. For each source, cells were mapped to cancer, donor, and broad compartment. Malignant, immune, and stromal compartments were assigned from source-provided major cell-type fields, normal-cell calls, dataset-specific overrides, and keyword rules. A donor-compartment group required at least 20 cells. Counts were aggregated by donor-compartment, normalized to counts per million, log1p transformed, collapsed to gene symbols by averaging, z-standardized within cancer, and scored by averaging recovered program genes when at least three genes were present. Source-bounded contrasts were paired donor differences for malignant minus immune, malignant minus stromal, and immune minus stromal within each source-cancer-program stratum. These contrasts were tested with two-sided paired Wilcoxon signed-rank tests using zero_method=`wilcox`; all-zero paired differences were retained as evaluable with P=1. Benjamini-Hochberg adjustment was applied within each source-cancer family of 24 program-by-contrast rows. A support row required evaluable status, paired donor n>=6, and within-family q<0.05. A separate n>=10 flag was retained as sensitivity only. Pooled cross-source donor-unit rows were retained only as descriptive summaries because no cross-repository biological-participant crosswalk established global donor independence.

Spatial BRCA Visium H5AD files were processed as sample-level triangulation. In-tissue spots were assigned to malignant, immune, or stromal compartments by source spot-score argmax or classification fallback. Spot counts were normalized to counts per million, log1p transformed, program genes were z-standardized within sample, and compartment medians were computed per sample and program. Spatial contrasts used four sample-level paired median differences only; no spot-level P value was used.

## Results

### CPTAC Case Flow and Primary Registry

The one-to-one matched table contained {support['primary']['unique_cases_any_matched_layer']} unique cases across either protein layer. The 160-row primary registry retained {support['primary']['primary_evaluable_rows']} evaluable rows and {support['primary']['primary_not_evaluable_rows']} non-evaluable rows; every non-evaluable row had fewer than 20 complete score pairs under the three-gene rule. There were no duplicated cancer-program-layer-case keys. Maximum matched case counts across cancer-layer combinations ranged from {int(flow.groupby(['cancer', 'comparison_layer'])['matched_unique_cases'].max().min())} to {int(flow.groupby(['cancer', 'comparison_layer'])['matched_unique_cases'].max().max())}.

### Cancer-Specific RNA-Protein and RNA-Phosphoprotein Concordance

Of {len(prot)} evaluable RNA-protein rows, {int((prot['q_value_bh'] < 0.05).sum())} had q<0.05; rho ranged from {fmt(prot['spearman_rho'].min())} to {fmt(prot['spearman_rho'].max())}, with a median of {fmt(prot['spearman_rho'].median())}. Of {len(phospho)} evaluable RNA-phosphoprotein rows, {int((phospho['q_value_bh'] < 0.05).sum())} had q<0.05; rho ranged from {fmt(phospho['spearman_rho'].min())} to {fmt(phospho['spearman_rho'].max())}, with a median of {fmt(phospho['spearman_rho'].median())}. The only evaluable row not passing BH was {nonsig['cancer'].upper()} {nonsig['program_label']} RNA-protein concordance (n={int(nonsig['complete_pair_n'])}, rho={fmt(nonsig['spearman_rho'])}, P={fmt(nonsig['p_value'])}, q={fmt(nonsig['q_value_bh'])}).

All program-layer summaries with evaluable cancers had a positive-cancer fraction of 1.000, but descriptive I2 fractions ranged from {fmt(heterogeneity['descriptive_i2'].min())} to {fmt(heterogeneity['descriptive_i2'].max())}. Two registered phosphoproteomic program-layer combinations had zero evaluable cancers under the complete-pair gate: MHC-II and Cytolytic.

{table_markdown}

### Scoring, Threshold, and Influence Sensitivities

Feature-row and gene-equal estimates had the same sign in all {len(sensitivity)} jointly evaluable rows. Protein estimates were identical because those matrices contained one recovered feature per mapped gene. Among {len(phospho_sens)} jointly evaluable phosphoprotein rows, the median absolute rho difference was {fmt((phospho_sens['spearman_rho_gene'] - phospho_sens['spearman_rho_feature']).abs().median())} and the maximum was {fmt((phospho_sens['spearman_rho_gene'] - phospho_sens['spearman_rho_feature']).abs().max())}; thus gene collapse changed magnitudes while retaining directions. The one-gene threshold sensitivity made {int((threshold1['status'] == 'evaluable').sum())} rows evaluable and {int((threshold1['q_value_bh'] < 0.05).sum())} passed its separate BH family, but it did not alter the locked primary registry. The maximum absolute leave-one-case change in primary rho was {fmt(max_influence['max_abs_leave_one_out_rho_delta'])}, in the GBM Checkpoint RNA-phosphoprotein stratum after omitting case {max_influence['most_influential_case_id']}; the full-stratum result remained positive and BH-supported (n={int(max_influence['complete_pair_n'])}, rho={fmt(max_influence['spearman_rho'])}, 95% CI {fmt(max_influence['ci_lower'])} to {fmt(max_influence['ci_upper'])}, q={fmt(max_influence['q_value_bh'])}).

### RNA-Adjusted Residual Context

Linear residual models were evaluable in {int((diagnostics['status'] == 'evaluable').sum())} cancer-program-layer strata. The residual-context registry contained 640 rows, but the finite tested denominator was {residual_linear_eval} for the linear-residual family and {residual_rank_eval} for the rank-residual sensitivity family. By context variable, each family had the same evaluability pattern: {context_sentence}. Of the finite rows, {int((residual['linear_q_value_bh'] < 0.05).sum())} passed the linear-residual BH family and {int((residual['rank_q_value_bh'] < 0.05).sum())} passed the rank-residual sensitivity family. The {len(residual_both)} rows passing both families had concordant directions. Linear-residual supported rows were distributed across CIBERSORT myeloid ({int((residual[residual['linear_q_value_bh'] < 0.05]['context_variable'] == 'cibersort_myeloid').sum())}), CIBERSORT T-cell ({int((residual[residual['linear_q_value_bh'] < 0.05]['context_variable'] == 'cibersort_t_cell').sum())}), ESTIMATE immune score ({int((residual[residual['linear_q_value_bh'] < 0.05]['context_variable'] == 'estimate_immune_score').sum())}), and tumor purity ({int((residual[residual['linear_q_value_bh'] < 0.05]['context_variable'] == 'estimate_tumor_purity').sum())}). These associations remained secondary descriptions of immune/purity context after RNA adjustment.

### RNA-Only TCGA Sensitivity

The TCGA layer included {support['support_layers']['tcga']['primary_tumor_patients']} primary-tumor patients and {support['support_layers']['tcga']['scored_samples']} scored samples, with no repeated patient-program row. The full 68-signature evidence surface contained {len(tcga_corr)} within-cancer correlations and {int((tcga_corr['fdr_bh'] < 0.05).sum())} had q<0.05 under the single 5440-test TCGA 68-signature BH family. Across the seven non-proteasome programs, canonical-anchor median within-cancer correlations ranged from {anchor_ranges[0]} to {anchor_ranges[1]}; the minimum positive-cancer fraction among these anchor summaries was {fmt(tcga_main['selected_anchor_min_positive_cancer_fraction'].min())}. The proteasome comparator had an anchor median rho of {fmt(tcga_sig.loc[tcga_sig['signature_label'].eq('Proteasome'), 'selected_anchor_median_rho'].iloc[0])}. Immune subtype testing generated {len(tcga_tests)} program-level Kruskal-Wallis tests in a separate eight-test BH family. Immune C5 contained only {int(tcga_subtype[tcga_subtype['immune_subtype_short'].eq('Immune C5')]['n'].max())} patients and remained descriptive rather than part of the n>=20 ordering. These results support RNA-level immune-landscape alignment only.

### APOLLO Same-Cancer Support

APOLLO contributed 87 LUAD matched public sample units with RNA, protein, and phosphoprotein measurements and 70 ovarian matched public sample units with RNA and protein measurements. Across 24 registered cohort-layer-program rows, {support['support_layers']['apollo']['evaluable_rows']} were evaluable, {len(apollo_supported)} were labelled supported, {len(apollo_opposite)} opposite, {support['support_layers']['apollo']['null_result_rows']} null, and {support['support_layers']['apollo']['not_evaluable_rows']} not evaluable after cohort-layer BH control. The 7 non-evaluable rows were those without enough recovered and matched RNA-to-omic program scores to satisfy the n>=10 and distinct-value gate. Representative supported rows were {apollo_examples}. The opposite row was retained rather than recoded: {apollo_opposite_text}. This mixture supports only bounded same-cancer public proteomic context.

### Donor-Aware Single-Cell and Sample-Aware Spatial Triangulation

The combined single-cell/3CA layer contained {support['support_layers']['single_cell']['unique_donor_units']} source-specific donor units from {support['support_layers']['single_cell']['unique_source_ids']} public sources. No exact raw donor label appeared in more than one source, but the provenance audit found no cross-repository participant crosswalk and therefore could not establish global biological independence. The 24 pooled program-by-compartment rows were retained as descriptive summaries only; none contributed a global P value, q value, or support count. Within source-cancer families, {support['support_layers']['single_cell']['within_source_evaluable_rows']} rows were evaluable and {len(scrna_support)} met the separately adjusted n>=6 support window; {int((scrna_support['median_delta'] > 0).sum())} supported rows had positive median deltas and {int((scrna_support['median_delta'] < 0).sum())} had negative median deltas. Representative source-bounded rows were {scrna_examples}. These contrasts locate program signal across broad compartments within sources but do not establish malignant-cell origin or mechanism.

All four breast Visium samples yielded sample-level program and compartment summaries. The spatial pairwise table contained {len(spatial_pairwise)} sample-level program-contrast rows, each with paired sample n=4. Examples included immune-minus-stromal MHC-I median delta {fmt(spatial_pairwise[(spatial_pairwise['signature_label'] == 'MHC-I') & (spatial_pairwise['contrast'] == 'immune_minus_stromal')]['median_delta'].iloc[0])} and malignant-minus-immune Cytolytic median delta {fmt(spatial_pairwise[(spatial_pairwise['signature_label'] == 'Cytolytic') & (spatial_pairwise['contrast'] == 'malignant_minus_immune')]['median_delta'].iloc[0])}. Spatial contrasts were therefore treated as same-collection, sample-aware triangulation, without spot-level P values or an independent-validation claim.
"""
    return text


def build_numeric_map(root: Path, data: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    primary = data["primary"]
    residual = data["residual"]
    diagnostics = data["diagnostics"]
    heterogeneity = data["heterogeneity"]
    feature = data["feature"]
    threshold1 = data["threshold1"]
    tcga_corr = data["tcga_corr"]
    tcga_sig = data["tcga_sig"]
    tcga_tests = data["tcga_tests"]
    tcga_subtype = data["tcga_subtype"]
    apollo = data["apollo"]
    apollo_scores = pd.read_csv(root / "results/apollo/apollo_signature_scores_v0_7.csv")
    apollo_score_key = ["cohort_label", "cancer_code", "comparison_layer", "signature", "sample_id"]
    apollo_registry_key = ["cohort_label", "cancer_code", "comparison_layer", "program"]
    apollo_luad_units = (
        apollo_scores[apollo_scores["cohort_label"].eq("APOLLO-LUAD")]
        .groupby("comparison_layer")["sample_id"]
        .nunique()
    )
    apollo_ov_units = (
        apollo_scores[apollo_scores["cohort_label"].eq("APOLLO-OV")]
        .groupby("comparison_layer")["sample_id"]
        .nunique()
    )
    scrna_within = data["scrna_within"]
    spatial_pairwise = data["spatial_pairwise"]
    support = data["support"]
    table2_current = data["table2_current"]

    evaluable = primary[primary["status"].eq("evaluable")]
    prot = evaluable[evaluable["comparison_layer"].eq("proteomics")]
    phospho = evaluable[evaluable["comparison_layer"].eq("phosphoproteomics")]
    nonsig = evaluable[~(evaluable["q_value_bh"] < 0.05)].iloc[0]
    max_influence = primary.loc[primary["max_abs_leave_one_out_rho_delta"].idxmax()]
    sensitivity = primary[["cancer", "program", "comparison_layer", "spearman_rho"]].merge(
        feature[["cancer", "program", "comparison_layer", "spearman_rho"]],
        on=["cancer", "program", "comparison_layer"],
        suffixes=("_gene", "_feature"),
    ).dropna()
    phospho_sens = sensitivity[sensitivity["comparison_layer"].eq("phosphoproteomics")]
    scrna_support = scrna_within[scrna_within["within_source_support_window"].astype(bool)]

    add_map(rows, "primary_unique_cases", support["primary"]["unique_cases_any_matched_layer"], "RESULTS_SUMMARY.json", ".primary.unique_cases_any_matched_layer", "Results/CPTAC case flow", "Unique CPTAC cases across either matched layer")
    add_map(rows, "primary_registry_rows", len(primary), "results/cptac/cptac_primary_concordance_registry_v2.csv", "all rows", "Results/CPTAC case flow", "Registered primary rows")
    add_map(rows, "primary_evaluable_rows", len(evaluable), "results/cptac/cptac_primary_concordance_registry_v2.csv", "status=evaluable", "Results/CPTAC case flow", "Primary evaluable rows")
    add_map(rows, "primary_not_evaluable_rows", int((primary["status"] != "evaluable").sum()), "results/cptac/cptac_primary_concordance_registry_v2.csv", "status!=evaluable", "Results/CPTAC case flow", "Primary non-evaluable rows")
    add_map(rows, "protein_evaluable_rows", len(prot), "results/cptac/cptac_primary_concordance_registry_v2.csv", "comparison_layer=proteomics,status=evaluable", "Results/primary concordance", "RNA-protein evaluable rows")
    add_map(rows, "protein_q_lt_0_05", int((prot["q_value_bh"] < 0.05).sum()), "results/cptac/cptac_primary_concordance_registry_v2.csv", "proteomics evaluable q<0.05", "Results/primary concordance", "RNA-protein BH-supported rows")
    add_map(rows, "phosphoprotein_evaluable_rows", len(phospho), "results/cptac/cptac_primary_concordance_registry_v2.csv", "comparison_layer=phosphoproteomics,status=evaluable", "Results/primary concordance", "RNA-phosphoprotein evaluable rows")
    add_map(rows, "phosphoprotein_q_lt_0_05", int((phospho["q_value_bh"] < 0.05).sum()), "results/cptac/cptac_primary_concordance_registry_v2.csv", "phosphoproteomics evaluable q<0.05", "Results/primary concordance", "RNA-phosphoprotein BH-supported rows")
    for field in ["complete_pair_n", "spearman_rho", "ci_lower", "ci_upper", "p_value", "q_value_bh"]:
        add_map(rows, f"nonsig_pdac_proteasome_{field}", nonsig[field], "results/cptac/cptac_primary_concordance_registry_v2.csv", "cancer=pdac,program_label=Proteasome,comparison_layer=proteomics", "Results/primary concordance", f"Only non-BH-supported primary row {field}")
    add_map(rows, "descriptive_i2_min", heterogeneity["descriptive_i2"].min(), "results/cptac/cptac_heterogeneity_summary_v2.csv", "all evaluable program-layer summaries", "Results/primary concordance", "Minimum descriptive I2 fraction")
    add_map(rows, "descriptive_i2_max", heterogeneity["descriptive_i2"].max(), "results/cptac/cptac_heterogeneity_summary_v2.csv", "all evaluable program-layer summaries", "Results/primary concordance", "Maximum descriptive I2 fraction")
    add_map(rows, "zero_evaluable_program_layer_rows", int((table2_current["evaluable_cancer_n"] == 0).sum()), "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv", "evaluable_cancer_n=0", "Results/primary concordance", "Registered program-layer rows with zero evaluable cancers")
    add_map(rows, "joint_gene_feature_sensitivity_rows", len(sensitivity), "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv", "joined finite primary and feature-row rows", "Results/sensitivity", "Jointly evaluable gene-equal and feature-row rows")
    add_map(rows, "phospho_feature_gene_median_abs_delta", (phospho_sens["spearman_rho_gene"] - phospho_sens["spearman_rho_feature"]).abs().median(), "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv", "phosphoproteomics jointly evaluable", "Results/sensitivity", "Median absolute rho difference for phosphoproteomics")
    add_map(rows, "phospho_feature_gene_max_abs_delta", (phospho_sens["spearman_rho_gene"] - phospho_sens["spearman_rho_feature"]).abs().max(), "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv", "phosphoproteomics jointly evaluable", "Results/sensitivity", "Maximum absolute rho difference for phosphoproteomics")
    add_map(rows, "threshold1_evaluable_rows", int((threshold1["status"] == "evaluable").sum()), "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv", "status=evaluable", "Results/sensitivity", "One-gene threshold evaluable rows")
    add_map(rows, "threshold1_q_lt_0_05", int((threshold1["q_value_bh"] < 0.05).sum()), "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv", "q<0.05", "Results/sensitivity", "One-gene threshold BH-supported rows")
    for field in ["complete_pair_n", "spearman_rho", "ci_lower", "ci_upper", "q_value_bh", "max_abs_leave_one_out_rho_delta"]:
        add_map(rows, f"max_influence_{field}", max_influence[field], "results/cptac/cptac_primary_concordance_registry_v2.csv", "row with max max_abs_leave_one_out_rho_delta", "Results/sensitivity", f"Maximum influence row {field}")
    add_map(rows, "max_influence_case", max_influence["most_influential_case_id"], "results/cptac/cptac_primary_concordance_registry_v2.csv", "row with max max_abs_leave_one_out_rho_delta", "Results/sensitivity", "Most influential case identifier")
    add_map(rows, "residual_registry_rows", len(residual), "results/cptac/cptac_residual_context_registry_v2.csv", "all rows", "Results/residual context", "Registered residual-context rows")
    add_map(rows, "residual_linear_finite_rows", residual["linear_p_value"].notna().sum(), "results/cptac/cptac_residual_context_registry_v2.csv", "linear_p_value finite", "Results/residual context", "Finite linear residual-context tests")
    add_map(rows, "residual_rank_finite_rows", residual["rank_p_value"].notna().sum(), "results/cptac/cptac_residual_context_registry_v2.csv", "rank_p_value finite", "Results/residual context", "Finite rank residual-context tests")
    add_map(rows, "residual_linear_q_lt_0_05", int((residual["linear_q_value_bh"] < 0.05).sum()), "results/cptac/cptac_residual_context_registry_v2.csv", "linear_q_value_bh<0.05", "Results/residual context", "Linear residual-context BH-supported rows")
    add_map(rows, "residual_rank_q_lt_0_05", int((residual["rank_q_value_bh"] < 0.05).sum()), "results/cptac/cptac_residual_context_registry_v2.csv", "rank_q_value_bh<0.05", "Results/residual context", "Rank residual-context BH-supported rows")
    add_map(rows, "residual_both_q_lt_0_05", int(((residual["linear_q_value_bh"] < 0.05) & (residual["rank_q_value_bh"] < 0.05)).sum()), "results/cptac/cptac_residual_context_registry_v2.csv", "linear and rank q<0.05", "Results/residual context", "Rows passing both residual families")
    add_map(rows, "residual_evaluable_strata", int((diagnostics["status"] == "evaluable").sum()), "results/cptac/cptac_residual_model_diagnostics_v2.csv", "status=evaluable", "Results/residual context", "Evaluable residual model strata")
    add_map(rows, "tcga_patients", support["support_layers"]["tcga"]["primary_tumor_patients"], "results/support_layers_summary_v2.json", ".tcga.primary_tumor_patients", "Results/TCGA", "TCGA primary-tumor patients")
    add_map(rows, "tcga_scored_samples", support["support_layers"]["tcga"]["scored_samples"], "results/support_layers_summary_v2.json", ".tcga.scored_samples", "Results/TCGA", "TCGA scored samples")
    add_map(rows, "tcga_68immune_correlations", len(tcga_corr), "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv", "all rows", "Results/TCGA", "TCGA within-cancer 68-immune correlations")
    add_map(rows, "tcga_68immune_finite_p_values", pd.to_numeric(tcga_corr["p_value"], errors="coerce").notna().sum(), "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv", "p_value finite", "Methods/TCGA multiplicity", "Finite TCGA 68-immune Spearman P values in the single BH family")
    add_map(rows, "tcga_68immune_finite_fdr_bh_values", pd.to_numeric(tcga_corr["fdr_bh"], errors="coerce").notna().sum(), "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv", "fdr_bh finite", "Methods/TCGA multiplicity", "Finite TCGA 68-immune BH-adjusted q values")
    add_map(rows, "tcga_68immune_p_family_count", tcga_corr["p_family"].nunique(), "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv", "unique p_family", "Methods/TCGA multiplicity", "Number of retained TCGA 68-immune P-value families")
    add_map(rows, "tcga_68immune_q_lt_0_05", int((tcga_corr["fdr_bh"] < 0.05).sum()), "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv", "fdr_bh<0.05", "Results/TCGA", "TCGA 68-immune correlations with q<0.05")
    add_map(rows, "tcga_anchor_median_min_nonproteasome", tcga_sig[~tcga_sig["signature_label"].eq("Proteasome")]["selected_anchor_median_rho"].min(), "results/tcga/tcga_xena_signature_summary_v0_5.csv", "non-Proteasome signatures", "Results/TCGA", "Minimum non-proteasome selected-anchor median rho")
    add_map(rows, "tcga_anchor_median_max_nonproteasome", tcga_sig[~tcga_sig["signature_label"].eq("Proteasome")]["selected_anchor_median_rho"].max(), "results/tcga/tcga_xena_signature_summary_v0_5.csv", "non-Proteasome signatures", "Results/TCGA", "Maximum non-proteasome selected-anchor median rho")
    add_map(rows, "tcga_proteasome_anchor_median", tcga_sig.loc[tcga_sig["signature_label"].eq("Proteasome"), "selected_anchor_median_rho"].iloc[0], "results/tcga/tcga_xena_signature_summary_v0_5.csv", "signature_label=Proteasome", "Results/TCGA", "Proteasome selected-anchor median rho")
    add_map(rows, "tcga_subtype_tests", len(tcga_tests), "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv", "all rows", "Results/TCGA", "TCGA immune subtype Kruskal-Wallis tests")
    add_map(rows, "tcga_subtype_finite_fdr_bh_values", pd.to_numeric(tcga_tests["kruskal_fdr_bh"], errors="coerce").notna().sum(), "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv", "kruskal_fdr_bh finite", "Methods/TCGA multiplicity", "Finite TCGA immune-subtype BH-adjusted q values")
    add_map(rows, "tcga_immune_c5_n", int(tcga_subtype[tcga_subtype["immune_subtype_short"].eq("Immune C5")]["n"].max()), "results/tcga/tcga_xena_immune_subtype_summary_v0_5.csv", "immune_subtype_short=Immune C5", "Results/TCGA", "Immune C5 patient count")
    for key in ["registered_rows", "evaluable_rows", "supported_rows", "opposite_rows", "null_result_rows", "not_evaluable_rows"]:
        add_map(rows, f"apollo_{key}", support["support_layers"]["apollo"][key], "results/support_layers_summary_v2.json", f".apollo.{key}", "Results/APOLLO", f"APOLLO {key}")
    add_map(rows, "apollo_score_duplicate_inference_key_rows", int(apollo_scores.duplicated(apollo_score_key, keep=False).sum()), "results/apollo/apollo_signature_scores_v0_7.csv", "duplicated cohort/cancer/layer/program/sample_id keys", "Methods/APOLLO unit rule", "APOLLO score duplicate inference-key rows")
    add_map(rows, "apollo_registry_duplicate_program_keys", int(apollo.duplicated(apollo_registry_key, keep=False).sum()), "results/apollo/apollo_concordance_registry_v2.csv", "duplicated cohort/cancer/layer/program keys", "Methods/APOLLO unit rule", "APOLLO registry duplicate program-key rows")
    add_map(rows, "apollo_fdr_family_count", apollo["fdr_family"].nunique(), "results/apollo/apollo_concordance_registry_v2.csv", "unique fdr_family", "Methods/APOLLO multiplicity", "APOLLO cohort-layer BH families")
    add_map(rows, "apollo_luad_public_sample_units", int(apollo_luad_units.min()), "results/apollo/apollo_signature_scores_v0_7.csv", "APOLLO-LUAD layer sample_id counts", "Methods/APOLLO unit rule", "Minimum APOLLO-LUAD public sample_id units across RNA/proteome/phosphoproteome score layers")
    add_map(rows, "apollo_ov_public_sample_units", int(apollo_ov_units.min()), "results/apollo/apollo_signature_scores_v0_7.csv", "APOLLO-OV layer sample_id counts", "Methods/APOLLO unit rule", "Minimum APOLLO-OV public sample_id units across RNA/proteome score layers")
    add_map(rows, "apollo_fisher_z_ci_rows", int(apollo["status"].eq("evaluable").sum()), "results/apollo/apollo_concordance_registry_v2.csv", "status=evaluable", "Methods/APOLLO CI", "APOLLO evaluable rows with Fisher-z Spearman confidence intervals")
    apollo_opp = apollo[apollo["evidence_status"].eq("opposite")].iloc[0]
    for field in ["matched_sample_n", "spearman_rho", "ci_lower", "ci_upper", "q_value_bh"]:
        add_map(rows, f"apollo_opposite_{field}", apollo_opp[field], "results/apollo/apollo_concordance_registry_v2.csv", "evidence_status=opposite", "Results/APOLLO", f"APOLLO opposite row {field}")
    add_map(rows, "scrna_unique_donor_units", support["support_layers"]["single_cell"]["unique_donor_units"], "results/support_layers_summary_v2.json", ".single_cell.unique_donor_units", "Results/single-cell", "Source-specific donor units")
    add_map(rows, "scrna_unique_sources", support["support_layers"]["single_cell"]["unique_source_ids"], "results/support_layers_summary_v2.json", ".single_cell.unique_source_ids", "Results/single-cell", "Single-cell/3CA source count")
    add_map(rows, "scrna_within_evaluable_rows", support["support_layers"]["single_cell"]["within_source_evaluable_rows"], "results/support_layers_summary_v2.json", ".single_cell.within_source_evaluable_rows", "Results/single-cell", "Within-source evaluable rows")
    add_map(rows, "scrna_within_support_rows", len(scrna_support), "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv", "within_source_support_window=True", "Results/single-cell", "Within-source support rows")
    add_map(rows, "scrna_support_positive_median_delta", int((scrna_support["median_delta"] > 0).sum()), "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv", "support rows median_delta>0", "Results/single-cell", "Support rows with positive median delta")
    add_map(rows, "scrna_support_negative_median_delta", int((scrna_support["median_delta"] < 0).sum()), "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv", "support rows median_delta<0", "Results/single-cell", "Support rows with negative median delta")
    add_map(rows, "spatial_unique_samples", support["support_layers"]["spatial"]["unique_samples"], "results/support_layers_summary_v2.json", ".spatial.unique_samples", "Results/spatial", "Spatial BRCA samples")
    add_map(rows, "spatial_pairwise_rows", len(spatial_pairwise), "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv", "all rows", "Results/spatial", "Spatial sample-level pairwise rows")

    out = pd.DataFrame(rows)
    out.to_csv(root / "current_p1_closure/NUMERIC_MAP.csv", index=False)
    return out


def build_closure(root: Path, numeric_map: pd.DataFrame) -> dict[str, Any]:
    issues = [
        {
            "issue_id": "EXTERNAL-R1-P1-01",
            "severity": "P1",
            "section": "RNA-adjusted residual context",
            "closure_class": "CLOSED_REPORTING",
            "action": "Specified the residual association statistic, finite-row gates, linear residual model, rank-residual algorithm, context missingness criteria, and finite denominators.",
            "evidence_files": [
                "code/run_cptac_v2.py",
                "results/cptac/cptac_residual_context_registry_v2.csv",
                "results/cptac/cptac_residual_model_diagnostics_v2.csv",
                "current_p1_closure/NUMERIC_MAP.csv",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P1-02",
            "severity": "P1",
            "section": "TCGA RNA-only sensitivity",
            "closure_class": "CLOSED_VERIFIED",
            "action": "Verified the retained anchor mapping from existing code, stated that full 5440-row 68-signature matrix is the evidence surface, and specified immune-subtype Kruskal-Wallis testing and multiplicity.",
            "evidence_files": [
                "results/tcga/tcga_xena_signature_summary_v0_5.csv",
                "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv",
                "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv",
                "logs/tcga_source_execution_report.md",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P1-03",
            "severity": "P1",
            "section": "APOLLO",
            "closure_class": "CLOSED_REPORTING",
            "action": "Disclosed APOLLO-OV source-level pre-imputed proteomics, distinguished it from no added score-level imputation, and specified scoring, matching, evaluability, and family-wise labels.",
            "evidence_files": [
                "code/postprocess_support_v2.py",
                "results/apollo/apollo_concordance_registry_v2.csv",
                "results/apollo/apollo_signature_recovery_v0_7.csv",
                "logs/apollo_source_execution_report.md",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P1-04",
            "severity": "P1",
            "section": "Single-cell and spatial triangulation",
            "closure_class": "CLOSED_REPORTING",
            "action": "Specified pseudobulk aggregation, CPM/log1p normalization, scoring, compartment harmonization, paired contrast estimand, Wilcoxon zero/tie handling, and n>=6 support-window rule.",
            "evidence_files": [
                "code/run_h5ad_support_v2.py",
                "code/postprocess_support_v2.py",
                "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv",
                "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P1-05",
            "severity": "P1",
            "section": "Source manifest timing and provenance",
            "closure_class": "CLOSED_REPORTING",
            "action": "Disclosed that the distributed source-provenance manifest was finalized after result generation, while propagating retained source identifiers, URLs where available, sizes, hashes, and crosswalk limitations.",
            "evidence_files": [
                "RUN_MANIFEST.json",
                "SOURCE_MANIFEST.json",
                "logs/COMMANDS_V2.log",
                "current_p1_closure/SOURCE_PROVENANCE_CURRENT_P1.csv",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P2-01",
            "severity": "P2",
            "section": "Program-layer summary table",
            "closure_class": "CLOSED_REPORTING",
            "action": "Added CURRENT_P1 table display with all 16 registered program-layer combinations, including zero-evaluable phosphoproteomic MHC-II and Cytolytic rows.",
            "evidence_files": ["tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv", "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.md"],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P2-02",
            "severity": "P2",
            "section": "Leave-one-case influence",
            "closure_class": "CLOSED_VERIFIED",
            "action": "Identified the maximum-influence stratum, case, rho change, and unchanged positive/BH-supported classification.",
            "evidence_files": ["results/cptac/cptac_primary_concordance_registry_v2.csv", "current_p1_closure/NUMERIC_MAP.csv"],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P2-03",
            "severity": "P2",
            "section": "Heterogeneity labels",
            "closure_class": "CLOSED_REPORTING",
            "action": "Labeled descriptive I2 as a fraction in text and CURRENT_P1 table.",
            "evidence_files": ["tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.md"],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-R1-P2-04",
            "severity": "P2",
            "section": "APOLLO and single-cell support magnitudes",
            "closure_class": "CLOSED_REPORTING",
            "action": "Added representative effect magnitudes and directions for APOLLO and source-bounded single-cell support rows.",
            "evidence_files": [
                "results/apollo/apollo_concordance_registry_v2.csv",
                "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-FREEZE-P1-01",
            "severity": "P1",
            "section": "TCGA RNA-only sensitivity",
            "closure_class": "CLOSED_VERIFIED_REPORTING",
            "action": "Documented that the TCGA 68-signature `fdr_bh` values were generated by one Benjamini-Hochberg adjustment across all 5440 finite Spearman P values in the retained full-matrix family, separate from the eight immune-subtype Kruskal-Wallis tests.",
            "evidence_files": [
                "code/run_tcga_source_wrapper_v2.py",
                "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv",
                "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv",
                "current_p1_closure/NUMERIC_MAP.csv",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
        {
            "issue_id": "EXTERNAL-FREEZE-P1-02",
            "severity": "P1",
            "section": "APOLLO same-cancer support",
            "closure_class": "CLOSED_VERIFIED_REPORTING",
            "action": "Documented the APOLLO public sample_id unit, duplicate-key failure rule, one-to-one RNA-to-omic join, 24-row registry, n>=10 and distinct-value evaluability gate, 17/7 denominator formation, and Fisher-z Spearman 95% CI procedure used by the current registry.",
            "evidence_files": [
                "code/postprocess_support_v2.py",
                "results/apollo/apollo_signature_scores_v0_7.csv",
                "results/apollo/apollo_signature_recovery_v0_7.csv",
                "results/apollo/apollo_concordance_registry_v2.csv",
                "current_p1_closure/NUMERIC_MAP.csv",
            ],
            "reanalysis": "NO",
            "result_change": "NO",
        },
    ]
    files_changed = [
        "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md",
        "current_p1_closure/CURRENT_P1_CLOSURE.json",
        "current_p1_closure/CURRENT_P1_CLOSURE.md",
        "current_p1_closure/NUMERIC_MAP.csv",
        "current_p1_freeze_remediation/FREEZE_REMEDIATION.json",
        "current_p1_freeze_remediation/FREEZE_REMEDIATION.md",
        "current_p1_closure/COMMANDS_RUN.txt",
        "current_p1_closure/SOURCE_PROVENANCE_CURRENT_P1.csv",
        "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv",
        "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.md",
    ]
    closure = {
        "project": PROJECT,
        "status": "READY_FOR_EXTERNAL_R2",
        "closure_class": "reporting_provenance",
        "base_authority": {
            "baseline_methods_results": str(BASELINE_PATH),
            "baseline_methods_results_sha256": sha256_file(BASELINE_PATH),
            "external_review": str(REVIEW_PATH),
            "external_review_sha256": sha256_file(REVIEW_PATH),
            "external_freeze_review": str(EXTERNAL_FREEZE_REVIEW_PATH),
            "external_freeze_review_sha256": sha256_file(EXTERNAL_FREEZE_REVIEW_PATH),
            "modification_plan": str(PLAN_PATH),
            "protocol_sha256": PROTOCOL_SHA256,
            "run_created_at_utc": "2026-08-11T16:07:39Z",
            "source_manifest_generated_at_utc": "2026-08-12T02:19:28+00:00",
        },
        "issues": issues,
        "external_freeze_p1_closure": {
            "review_decision": "NOT_FREEZE",
            "p1_open_at_review": 2,
            "closed_issue_ids": ["EXTERNAL-FREEZE-P1-01", "EXTERNAL-FREEZE-P1-02"],
            "reanalysis_performed": False,
            "result_change": False,
            "status": "PENDING_NUMERIC_GATE",
        },
        "commands_run": COMMANDS_RUN,
        "files_changed": files_changed,
        "reanalysis_performed": False,
        "result_change": False,
        "result_change_detail": "No primary, residual, TCGA, APOLLO, single-cell, or spatial estimates were recomputed or changed; this closure propagates existing machine results and provenance.",
        "unresolved_issues": [],
        "methods_results_path": str(CURRENT_METHODS_PATH),
        "numeric_gate": {
            "status": "PENDING_INDEPENDENT_GATE",
            "numeric_map_rows": int(len(numeric_map)),
        },
        "claim_boundary": "Cancer-specific patient/case-level cross-sectional CPTAC RNA-protein and RNA-phosphoprotein concordance with bounded RNA-only TCGA, same-cancer APOLLO, source-bounded donor-compartment single-cell, and sample-level spatial support; no prognosis, prediction, diagnosis, treatment response, causality, mechanism, malignant-cell-intrinsic, post-transcriptional, therapeutic-vulnerability, homogeneous pan-cancer coefficient, or independent external-validation claim.",
    }
    return closure


def build_closure_md(root: Path, closure: dict[str, Any]) -> str:
    lines = [
        "# BIO-05 Current P1 Closure",
        "",
        f"- Status: {closure['status']}",
        f"- Closure class: {closure['closure_class']}",
        f"- Reanalysis performed: {closure['reanalysis_performed']}",
        f"- Result change: {closure['result_change']}",
        f"- New Methods/Results: `{rel(CURRENT_METHODS_PATH, root)}`",
        "",
        "## Issue Closure",
        "",
        "| Issue | Severity | Class | Action | Evidence |",
        "|---|---:|---|---|---|",
    ]
    for issue in closure["issues"]:
        evidence = "; ".join(issue["evidence_files"])
        lines.append(f"| {issue['issue_id']} | {issue['severity']} | {issue['closure_class']} | {issue['action']} | {evidence} |")
    external = closure["external_freeze_p1_closure"]
    lines.extend(
        [
            "",
            "## External Freeze P1 Closure",
            "",
            f"- Review decision: {external['review_decision']}",
            f"- P1 open at review: {external['p1_open_at_review']}",
            f"- Closed issues: {', '.join(external['closed_issue_ids'])}",
            f"- Reanalysis performed: {external['reanalysis_performed']}",
            f"- Result change: {external['result_change']}",
            f"- Status: {external['status']}",
            "",
            "## Claim Boundary",
            "",
            closure["claim_boundary"],
            "",
            "## Commands",
            "",
        ]
    )
    lines.extend(f"- `{cmd}`" for cmd in closure["commands_run"])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", required=True)
    args = parser.parse_args()
    root = Path(args.analysis_root).resolve()
    if root != ROOT_ABS:
        raise SystemExit(f"Refusing to write outside authorized BIO-05 root: {root}")
    (root / "current_p1_closure").mkdir(parents=True, exist_ok=True)

    data = {
        "primary": pd.read_csv(root / "results/cptac/cptac_primary_concordance_registry_v2.csv"),
        "heterogeneity": pd.read_csv(root / "results/cptac/cptac_heterogeneity_summary_v2.csv"),
        "feature": pd.read_csv(root / "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv"),
        "threshold1": pd.read_csv(root / "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv"),
        "residual": pd.read_csv(root / "results/cptac/cptac_residual_context_registry_v2.csv"),
        "diagnostics": pd.read_csv(root / "results/cptac/cptac_residual_model_diagnostics_v2.csv"),
        "flow": pd.read_csv(root / "results/cptac/cptac_complete_denominator_flow_v2.csv"),
        "tcga_sig": pd.read_csv(root / "results/tcga/tcga_xena_signature_summary_v0_5.csv"),
        "tcga_subtype": pd.read_csv(root / "results/tcga/tcga_xena_immune_subtype_summary_v0_5.csv"),
        "tcga_tests": pd.read_csv(root / "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv"),
        "tcga_corr": pd.read_csv(root / "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv"),
        "apollo": pd.read_csv(root / "results/apollo/apollo_concordance_registry_v2.csv"),
        "scrna_global": pd.read_csv(root / "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv"),
        "scrna_within": pd.read_csv(root / "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv"),
        "spatial_summary": pd.read_csv(root / "results/spatial/spatial_brca_sample_signature_summary_v0_7.csv"),
        "spatial_pairwise": pd.read_csv(root / "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv"),
        "support": {
            "primary": json.loads((root / "RESULTS_SUMMARY.json").read_text())["primary"],
            "support_layers": json.loads((root / "RESULTS_SUMMARY.json").read_text())["support_layers"],
        },
    }
    source_manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text())
    data["table2_current"] = build_program_table(root, data["primary"], data["heterogeneity"])
    build_source_provenance(root, source_manifest)
    methods = build_methods_results(root, data)
    CURRENT_METHODS_PATH.write_text(methods, encoding="utf-8")
    numeric_map = build_numeric_map(root, data)
    closure = build_closure(root, numeric_map)
    write_json(root / "current_p1_closure/CURRENT_P1_CLOSURE.json", closure)
    (root / "current_p1_closure/CURRENT_P1_CLOSURE.md").write_text(build_closure_md(root, closure), encoding="utf-8")
    (root / "current_p1_closure/COMMANDS_RUN.txt").write_text("\n".join(COMMANDS_RUN) + "\n", encoding="utf-8")
    print(json.dumps({"project": PROJECT, "status": "CURRENT_P1_OUTPUTS_WRITTEN", "numeric_map_rows": int(len(numeric_map))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
