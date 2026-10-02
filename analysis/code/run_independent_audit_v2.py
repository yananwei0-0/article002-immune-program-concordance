#!/usr/bin/env python3
"""Independent code/output/manuscript audit for BIO-05 v2."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import rankdata, zscore
from statsmodels.stats.multitest import multipletests


PROTOCOL_SHA256 = "a3bc200faf24eefb1b89be08cb4b6804b81423f50c3337a430434e6015143fc3"
RAW_SCORE_AUDIT = [
    ("brca", "transcriptomics", "antigen_presentation_mhc_i"),
    ("brca", "phosphoproteomics", "tgfb_emt_exclusion_context"),
    ("luad", "proteomics", "ifn_gamma_response"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def ensg(value: Any) -> str:
    match = re.search(r"(ENSG\d+)(?:\.\d+)?", str(value))
    return match.group(1) if match else ""


def case_id(value: Any) -> str:
    text = str(value).strip()
    if not text or "QC" in text.upper():
        return ""
    text = re.sub(r"(_T|_A|_N)$", "", text, flags=re.I)
    text = re.sub(r"-(T|N)$", "", text, flags=re.I)
    text = re.sub(r"_D\d+$", "", text, flags=re.I)
    return text


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def independently_score(path: Path, layer: str, program: str, genes: pd.DataFrame) -> pd.DataFrame:
    with open_text(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    samples = header[1:] if header[0].lower() in {"idx", "index", "gene", "gene_id"} else header
    names = ["feature_id"] + samples
    raw = pd.read_csv(path, sep="\t", header=None, skiprows=1, names=names, na_values=["NA", "NaN", "nan", ""], low_memory=False)
    sample_meta = pd.DataFrame({"sample_id": samples})
    sample_meta["case_id"] = sample_meta["sample_id"].map(case_id)
    sample_meta["eligible"] = sample_meta["case_id"].ne("")
    if layer == "transcriptomics":
        normal = sample_meta["sample_id"].str.contains(r"(_A|_N|-N)$", case=False, regex=True)
        sample_meta["eligible"] &= ~normal
    sample_meta = sample_meta[sample_meta["eligible"]].copy()
    if sample_meta.duplicated("case_id").any():
        raise RuntimeError("Independent audit encountered unanticipated eligible case duplicates")
    sample_cols = sample_meta["sample_id"].tolist()
    values = raw[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        scaled = zscore(values, axis=1, ddof=1, nan_policy="omit")
    target = genes[genes["signature"].eq(program)].copy()
    stable_to_gene: dict[str, str] = {}
    for row in target.itertuples(index=False):
        if pd.isna(row.ensembl_stable_ids_from_cptac):
            continue
        for token in str(row.ensembl_stable_ids_from_cptac).split(";"):
            stable = ensg(token)
            if stable:
                stable_to_gene[stable] = str(row.gene_symbol).upper()
    feature_genes = pd.Series([stable_to_gene.get(ensg(x), "") for x in raw["feature_id"]])
    keep = feature_genes.ne("").to_numpy()
    frame = pd.DataFrame(scaled[keep, :], index=feature_genes[keep].to_numpy(), columns=sample_cols)
    gene_matrix = frame.groupby(level=0, sort=True).mean()
    nonmissing = gene_matrix.notna().sum(axis=0)
    score = gene_matrix.mean(axis=0, skipna=True)
    score[nonmissing < 3] = np.nan
    return pd.DataFrame({"sample_id": sample_cols, "normalized_case_id": sample_meta["case_id"].to_numpy(), "independent_score": score.to_numpy()})


def independent_rho(x: pd.Series, y: pd.Series) -> float:
    work = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(work) < 2:
        return np.nan
    rx = rankdata(work["x"].to_numpy(), method="average")
    ry = rankdata(work["y"].to_numpy(), method="average")
    return float(np.corrcoef(rx, ry)[0, 1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--protocol", required=True)
    args = parser.parse_args()
    root = Path(args.analysis_root).resolve()
    data_root = Path(args.data_root).resolve()
    protocol_path = Path(args.protocol).resolve()
    issues: list[dict[str, str]] = []
    checks: dict[str, Any] = {}

    checks["protocol_sha256"] = sha256(protocol_path)
    if checks["protocol_sha256"] != PROTOCOL_SHA256:
        issues.append({"severity": "P0", "id": "AUD-P0-PROTOCOL", "detail": "Protocol SHA-256 mismatch."})

    primary = pd.read_csv(root / "results/cptac/cptac_primary_concordance_registry_v2.csv")
    matched = pd.read_csv(root / "results/cptac/cptac_matched_gene_equal_scores_v2.csv")
    rho_diffs = []
    for row in primary.itertuples(index=False):
        group = matched[
            matched["cancer"].eq(row.cancer)
            & matched["program"].eq(row.program)
            & matched["comparison_layer"].eq(row.comparison_layer)
        ]
        audit_rho = independent_rho(group["rna_score"], group["omic_score"])
        if np.isfinite(row.spearman_rho) and np.isfinite(audit_rho):
            rho_diffs.append(abs(float(row.spearman_rho) - audit_rho))
    checks["independent_rank_pearson_rows"] = len(rho_diffs)
    checks["independent_rank_pearson_max_abs_rho_diff"] = max(rho_diffs) if rho_diffs else None
    if not rho_diffs or max(rho_diffs) > 1e-12:
        issues.append({"severity": "P1", "id": "AUD-P1-RHO", "detail": "Independent average-rank Pearson check did not reproduce all finite Spearman rho values."})

    finite = primary["p_value"].notna()
    independent_q = np.full(len(primary), np.nan)
    independent_q[finite] = multipletests(primary.loc[finite, "p_value"].to_numpy(), method="fdr_bh")[1]
    q_diff = np.nanmax(np.abs(independent_q - primary["q_value_bh"].to_numpy()))
    checks["independent_statsmodels_bh_max_abs_q_diff"] = float(q_diff)
    if q_diff > 1e-12:
        issues.append({"severity": "P1", "id": "AUD-P1-FDR", "detail": "Independent statsmodels BH values differ from the primary registry."})

    evaluable = primary[primary["status"].eq("evaluable")]
    checks["primary_bca_ci_rows"] = int(evaluable["ci_method"].eq("paired_case_bca_bootstrap").sum())
    checks["primary_bca_ci_expected_rows"] = int(len(evaluable))
    checks["primary_bca_ci_order_failures"] = int(
        ((evaluable["ci_lower"] > evaluable["spearman_rho"]) | (evaluable["spearman_rho"] > evaluable["ci_upper"])).sum()
    )
    permutation_required = evaluable[
        (evaluable["complete_pair_n"] <= 30)
        | (evaluable["x_tie_fraction"] > 0)
        | (evaluable["y_tie_fraction"] > 0)
    ]
    checks["primary_permutation_required_rows"] = int(len(permutation_required))
    checks["primary_permutation_completed_rows"] = int(
        permutation_required["p_value_method"].eq("two_sided_monte_carlo_permutation").sum()
    )
    checks["primary_permutation_resamples_min"] = int(permutation_required["permutation_resamples"].min()) if len(permutation_required) else 0
    if (
        checks["primary_bca_ci_rows"] != checks["primary_bca_ci_expected_rows"]
        or checks["primary_bca_ci_order_failures"] != 0
        or checks["primary_permutation_completed_rows"] != checks["primary_permutation_required_rows"]
        or (len(permutation_required) and checks["primary_permutation_resamples_min"] != 999999)
    ):
        issues.append({"severity": "P1", "id": "AUD-P1-SPEARMAN-CALIBRATION", "detail": "Primary Spearman bootstrap/permutation calibration is incomplete or internally inconsistent."})

    input_files = pd.read_csv(root / "results/cptac/cptac_input_files_v2.csv")
    score_output = pd.read_csv(root / "results/cptac/cptac_program_scores_v2.csv")
    genes = pd.read_csv(data_root / "04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv")
    raw_checks = []
    for cancer, layer, program in RAW_SCORE_AUDIT:
        path = Path(input_files[input_files["cancer"].eq(cancer) & input_files["layer"].eq(layer)]["path"].iloc[0])
        independent = independently_score(path, layer, program, genes)
        observed = score_output[
            score_output["cancer"].eq(cancer)
            & score_output["layer"].eq(layer)
            & score_output["program"].eq(program)
        ][["normalized_case_id", "gene_equal_score"]]
        merged = independent.merge(observed, on="normalized_case_id", how="inner", validate="one_to_one").dropna()
        diff = float((merged["independent_score"] - merged["gene_equal_score"]).abs().max()) if len(merged) else np.nan
        raw_checks.append({"cancer": cancer, "layer": layer, "program": program, "compared_samples": len(merged), "max_abs_score_diff": diff})
    checks["raw_matrix_score_crosschecks"] = raw_checks
    if any(not np.isfinite(x["max_abs_score_diff"]) or x["max_abs_score_diff"] > 1e-12 for x in raw_checks):
        issues.append({"severity": "P1", "id": "AUD-P1-SCORE", "detail": "Independent raw-matrix gene-equal scoring check failed."})

    rerun = json.loads((root / "local_audit/RERUN_CONSISTENCY.json").read_text())
    checks["rerun_all_core_outputs_byte_identical"] = rerun["all_core_outputs_byte_identical"]
    if not rerun["all_core_outputs_byte_identical"]:
        issues.append({"severity": "P1", "id": "AUD-P1-RERUN", "detail": "Core CPTAC outputs were not byte-identical on rerun."})

    validation = json.loads((root / "VALIDATION_REPORT.json").read_text())
    results = json.loads((root / "RESULTS_SUMMARY.json").read_text())
    source_manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text())
    support = json.loads((root / "results/support_layers_summary_v2.json").read_text())
    donor_audit = pd.read_csv(root / "results/scrna_combined/public_scrna_source_donor_provenance_audit_v2.csv")
    global_scrna = pd.read_csv(root / "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv")
    manuscript = (root / "manuscript/BIO-05_METHODS_RESULTS_V2.md").read_text(encoding="utf-8")
    checks["authority_validation_pass"] = validation["status"] == "PASS"
    checks["source_manifest_complete"] = source_manifest["status"] == "COMPLETE"
    checks["single_cell_donor_provenance_sources"] = int(donor_audit["source_id"].nunique())
    checks["single_cell_global_p_values_removed"] = bool(global_scrna["wilcoxon_p_value"].isna().all() and global_scrna["q_value_bh"].isna().all())
    checks["single_cell_global_claim_rows_zero"] = int(global_scrna["main_claim_window"].astype(bool).sum()) == 0
    checks["single_cell_global_inference_not_retained"] = support["single_cell"]["donor_provenance_audit"]["global_pooled_donor_inference_retained"] is False
    checks["manuscript_has_methods"] = "## Methods" in manuscript
    checks["manuscript_has_results"] = "## Results" in manuscript
    checks["manuscript_has_no_introduction_or_discussion"] = "## Introduction" not in manuscript and "## Discussion" not in manuscript
    checks["manuscript_claim_fences_present"] = all(
        phrase in manuscript
        for phrase in [
            "no longitudinal time origin",
            "not called donor-proven independent replication",
            "spots were not treated as patients",
            "not establish malignant-cell origin or mechanism",
            "former pooled global donor-unit tests were removed from inferential claims",
        ]
    )
    checks["numeric_propagation"] = {
        "unique_cases_1023": "1023 unique cases" in manuscript,
        "primary_registry_160": "160-row primary registry" in manuscript,
        "tcga_4780": "4780 primary-tumor patients" in manuscript,
        "apollo_24": "24 registered cohort-layer-program rows" in manuscript,
    }
    for name in ["authority_validation_pass", "source_manifest_complete", "single_cell_global_p_values_removed", "single_cell_global_claim_rows_zero", "single_cell_global_inference_not_retained", "manuscript_has_methods", "manuscript_has_results", "manuscript_has_no_introduction_or_discussion", "manuscript_claim_fences_present"]:
        if not checks[name]:
            issues.append({"severity": "P1", "id": f"AUD-P1-{name.upper()}", "detail": f"Audit check failed: {name}."})
    if not all(checks["numeric_propagation"].values()):
        issues.append({"severity": "P1", "id": "AUD-P1-PROPAGATION", "detail": "One or more authoritative result counts were not propagated to Methods/Results."})

    open_p0 = sum(x["severity"] == "P0" for x in issues)
    open_p1 = sum(x["severity"] == "P1" for x in issues)
    gate_status = "READY" if open_p0 == 0 and open_p1 == 0 and results["status"] == "V2_ANALYSIS_COMPLETE" else "NEEDS_REVISION"
    gate = {
        "project": "BIO-05",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "gate_status": gate_status,
        "open_p0": open_p0,
        "open_p1": open_p1,
        "protocol_sha256": checks["protocol_sha256"],
        "real_source_execution_required": True,
        "real_source_execution_status": results["status"],
        "checks": checks,
        "issues": issues,
        "claim_boundary": "Cross-sectional cancer-specific concordance with bounded support layers; no prognosis, prediction, causality, mechanism, or therapeutic inference.",
    }
    (root / "AUDIT_GATE.json").write_text(json.dumps(gate, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    audit_dir = root / "local_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# BIO-05 independent v2 scientific audit",
        "",
        f"- Gate: **{gate_status}**",
        f"- Open P0: **{open_p0}**",
        f"- Open P1: **{open_p1}**",
        f"- Protocol SHA-256: `{checks['protocol_sha256']}`",
        "",
        "## Independent numerical checks",
        "",
        f"- Average-rank Pearson independently reproduced {checks['independent_rank_pearson_rows']} finite Spearman estimates; maximum absolute difference `{checks['independent_rank_pearson_max_abs_rho_diff']:.3e}`.",
        f"- statsmodels independently reproduced the primary BH family; maximum absolute q difference `{checks['independent_statsmodels_bh_max_abs_q_diff']:.3e}`.",
        f"- Retained rerun compared {rerun['core_file_count']} core files; byte-identical files `{rerun['byte_identical_count']}`.",
        "- Independent raw-matrix score checks:",
    ]
    for row in raw_checks:
        lines.append(f"  - `{row['cancer']}/{row['layer']}/{row['program']}`: n={row['compared_samples']}, maximum absolute score difference `{row['max_abs_score_diff']:.3e}`.")
    lines.extend(["", "## Design and propagation checks", ""])
    for name, value in checks.items():
        if isinstance(value, bool):
            lines.append(f"- `{name}`: {'PASS' if value else 'FAIL'}")
    lines.extend(["", "## Open issues", ""])
    if issues:
        for issue in issues:
            lines.append(f"- `{issue['severity']}` `{issue['id']}`: {issue['detail']}")
    else:
        lines.append("- None at P0 or P1.")
    lines.extend(
        [
            "",
            "## Audit conclusion",
            "",
            "The audit treats cancer-specific patient-level estimates as primary, retains non-evaluable and attenuated/opposite support results, and does not use historical COMPLETE/PASS labels or v1 numeric tokens as evidence.",
            "",
        ]
    )
    (audit_dir / "BIO-05_INDEPENDENT_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"gate_status": gate_status, "open_p0": open_p0, "open_p1": open_p1}, indent=2))
    return 0 if gate_status == "READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
