#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT = "BIO-05"
ROOT_ABS = Path(expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/manuscript27_current_p1_closure_20260815/analyses/BIO-05"))
METHODS_REL = "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(clean(payload), indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")


def update_closure_md(root: Path, status: str, failures: list[str]) -> None:
    md_path = root / "current_p1_closure/CURRENT_P1_CLOSURE.md"
    text = md_path.read_text(encoding="utf-8")
    external_status = "READY_FOR_EXTERNAL_R2" if status == "PASS" and not failures else "BLOCKED_REQUIRED_SOURCE"
    text = text.replace("- Status: PENDING_NUMERIC_GATE", f"- Status: {external_status}")
    base = text.split("\n## Numeric Gate\n", 1)[0].rstrip()
    failure_text = "None" if not failures else "; ".join(failures)
    section = (
        "\n\n## Numeric Gate\n\n"
        f"- Status: {status}\n"
        f"- Failures: {failure_text}\n"
        "- Gate output: `current_p1_closure/NUMERIC_GATE_CURRENT_P1.json`\n"
    )
    md_path.write_text(base + section, encoding="utf-8")


def fail_if(condition: bool, failures: list[str], label: str) -> None:
    if condition:
        failures.append(label)


def contains_all(text: str, tokens: list[str]) -> bool:
    return all(token in text for token in tokens)


def bh_fdr(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
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


def write_freeze_remediation(root: Path, gate: dict[str, Any], closure: dict[str, Any], failures: list[str]) -> None:
    out_dir = root / "current_p1_freeze_remediation"
    out_dir.mkdir(parents=True, exist_ok=True)
    freeze_issues = [issue for issue in closure.get("issues", []) if str(issue.get("issue_id", "")).startswith("EXTERNAL-FREEZE-P1")]
    evidence_files = sorted(
        {
            "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md",
            "current_p1_closure/CURRENT_P1_CLOSURE.json",
            "current_p1_closure/CURRENT_P1_CLOSURE.md",
            "current_p1_closure/NUMERIC_MAP.csv",
            "current_p1_closure/NUMERIC_GATE_CURRENT_P1.json",
            *[path for issue in freeze_issues for path in issue.get("evidence_files", [])],
        }
    )
    status = "READY_FOR_RESUBMISSION" if gate["status"] == "PASS" and not failures else "BLOCKED_REQUIRED_SOURCE"
    payload = {
        "project": PROJECT,
        "generated_at_utc": now_utc(),
        "review_decision": "NOT_FREEZE",
        "issues": freeze_issues,
        "evidence_files": evidence_files,
        "actions": [issue["action"] for issue in freeze_issues],
        "recalculation_performed": {
            "source_data_or_model_rerun": False,
            "endpoint_result_recalculation": False,
            "deterministic_gate_verification": [
                "TCGA standard BH recomputed from existing 5440 Spearman P values for verification only.",
                "APOLLO Fisher-z confidence intervals recomputed from existing registry rho/n for verification only.",
            ],
        },
        "result_change": False,
        "result_change_detail": "No TCGA or APOLLO estimates, P values, q values, labels, or source-model outputs were changed.",
        "files_changed": sorted(
            set(closure.get("files_changed", []))
            | {
                "current_p1_freeze_remediation/FREEZE_REMEDIATION.json",
                "current_p1_freeze_remediation/FREEZE_REMEDIATION.md",
            }
        ),
        "numeric_gate": {
            "status": gate["status"],
            "failures": failures,
            "path": "current_p1_closure/NUMERIC_GATE_CURRENT_P1.json",
            "external_freeze_checks": {
                key: value
                for key, value in gate["checks"].items()
                if key.startswith("tcga_") or key.startswith("apollo_") or key.startswith("text_contains_freeze")
            },
        },
        "unresolved_issues": failures,
        "status": status,
    }
    write_json(out_dir / "FREEZE_REMEDIATION.json", payload)

    lines = [
        "# BIO-05 Current-P1 Freeze Remediation",
        "",
        f"- Project: {PROJECT}",
        f"- Review decision: {payload['review_decision']}",
        f"- Status: {payload['status']}",
        f"- Result change: {payload['result_change']}",
        f"- Numeric gate: {gate['status']}",
        f"- Unresolved issues: {'None' if not failures else '; '.join(failures)}",
        "",
        "## Issues",
        "",
        "| Issue | Action | Evidence |",
        "|---|---|---|",
    ]
    for issue in freeze_issues:
        lines.append(
            f"| {issue['issue_id']} | {issue['action']} | {'; '.join(issue.get('evidence_files', []))} |"
        )
    lines.extend(
        [
            "",
            "## Recalculation",
            "",
            "- Source data/model rerun: False",
            "- Endpoint result recalculation: False",
            "- Gate verification: TCGA BH from existing P values; APOLLO Fisher-z CI from existing rho/n.",
            "",
            "## Files Changed",
            "",
        ]
    )
    lines.extend(f"- `{path}`" for path in payload["files_changed"])
    lines.append("")
    (out_dir / "FREEZE_REMEDIATION.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", required=True)
    args = parser.parse_args()
    root = Path(args.analysis_root).resolve()
    if root != ROOT_ABS:
        raise SystemExit(f"Refusing to run outside authorized BIO-05 root: {root}")

    failures: list[str] = []
    checks: dict[str, Any] = {}

    methods_path = root / METHODS_REL
    text = methods_path.read_text(encoding="utf-8")
    primary = pd.read_csv(root / "results/cptac/cptac_primary_concordance_registry_v2.csv")
    residual = pd.read_csv(root / "results/cptac/cptac_residual_context_registry_v2.csv")
    diagnostics = pd.read_csv(root / "results/cptac/cptac_residual_model_diagnostics_v2.csv")
    tcga_corr = pd.read_csv(root / "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv")
    tcga_tests = pd.read_csv(root / "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv")
    apollo = pd.read_csv(root / "results/apollo/apollo_concordance_registry_v2.csv")
    apollo_scores = pd.read_csv(root / "results/apollo/apollo_signature_scores_v0_7.csv")
    scrna_within = pd.read_csv(root / "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv")
    spatial_pairwise = pd.read_csv(root / "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv")
    table2 = pd.read_csv(root / "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv")
    numeric_map = pd.read_csv(root / "current_p1_closure/NUMERIC_MAP.csv")
    closure_path = root / "current_p1_closure/CURRENT_P1_CLOSURE.json"
    closure = json.loads(closure_path.read_text())

    evaluable = primary[primary["status"].eq("evaluable")]
    prot = evaluable[evaluable["comparison_layer"].eq("proteomics")]
    phospho = evaluable[evaluable["comparison_layer"].eq("phosphoproteomics")]
    max_influence = primary.loc[primary["max_abs_leave_one_out_rho_delta"].idxmax()]
    scrna_support = scrna_within[scrna_within["within_source_support_window"].astype(bool)]
    tcga_p = pd.to_numeric(tcga_corr["p_value"], errors="coerce").to_numpy(dtype=float)
    tcga_q = pd.to_numeric(tcga_corr["fdr_bh"], errors="coerce").to_numpy(dtype=float)
    tcga_q_expected = bh_fdr(tcga_p)
    tcga_q_mask = np.isfinite(tcga_q) & np.isfinite(tcga_q_expected)
    tcga_q_max_abs_diff = (
        float(np.max(np.abs(tcga_q[tcga_q_mask] - tcga_q_expected[tcga_q_mask]))) if tcga_q_mask.any() else np.nan
    )
    apollo_evaluable = apollo[apollo["status"].eq("evaluable")].copy()
    apollo_rho = pd.to_numeric(apollo_evaluable["spearman_rho"], errors="coerce").to_numpy(dtype=float)
    apollo_n = pd.to_numeric(apollo_evaluable["matched_sample_n"], errors="coerce").to_numpy(dtype=float)
    apollo_z = np.arctanh(np.clip(apollo_rho, -0.999999, 0.999999))
    apollo_se = 1 / np.sqrt(apollo_n - 3)
    apollo_ci_low_expected = np.tanh(apollo_z - 1.9599639845 * apollo_se)
    apollo_ci_high_expected = np.tanh(apollo_z + 1.9599639845 * apollo_se)
    apollo_ci_low = pd.to_numeric(apollo_evaluable["ci_lower"], errors="coerce").to_numpy(dtype=float)
    apollo_ci_high = pd.to_numeric(apollo_evaluable["ci_upper"], errors="coerce").to_numpy(dtype=float)
    apollo_ci_max_abs_diff = float(
        max(
            np.max(np.abs(apollo_ci_low - apollo_ci_low_expected)),
            np.max(np.abs(apollo_ci_high - apollo_ci_high_expected)),
        )
    )
    apollo_score_key = ["cohort_label", "cancer_code", "comparison_layer", "signature", "sample_id"]
    apollo_registry_key = ["cohort_label", "cancer_code", "comparison_layer", "program"]

    checks["methods_results_exists"] = methods_path.is_file() and methods_path.stat().st_size > 0
    checks["no_absolute_local_paths_in_methods"] = all(token not in text for token in ["/Users/", "/Volumes/", "Desktop/research"])
    checks["no_internal_review_terms_in_methods"] = all(token not in text for token in ["EXTERNAL-R1", "CURRENT_P1", "CLI", "reviewer packet", "current_p1_closure"])
    checks["no_todo_tbd_author_placeholders"] = not any(
        token in text.lower() for token in ["todo", "tbd", "[author", "author contribution", "xxx"]
    )
    checks["contains_methods_and_results_only"] = contains_all(text, ["## Methods", "## Results"]) and "## Discussion" not in text and "## Introduction" not in text
    checks["claim_boundary_negative_present"] = contains_all(text, ["not prognosis", "not donor-proven independent replication", "without spot-level P values"])
    checks["primary_registry_160"] = len(primary) == 160
    checks["primary_evaluable_129"] = len(evaluable) == 129
    checks["protein_77_76"] = len(prot) == 77 and int((prot["q_value_bh"] < 0.05).sum()) == 76
    checks["phospho_52_52"] = len(phospho) == 52 and int((phospho["q_value_bh"] < 0.05).sum()) == 52
    checks["residual_denominators"] = len(residual) == 640 and residual["linear_p_value"].notna().sum() == 516 and residual["rank_p_value"].notna().sum() == 516
    checks["residual_significance_counts"] = int((residual["linear_q_value_bh"] < 0.05).sum()) == 35 and int((residual["rank_q_value_bh"] < 0.05).sum()) == 52
    checks["residual_context_by_variable"] = residual.groupby("context_variable")["linear_p_value"].apply(lambda s: int(s.notna().sum())).to_dict()
    checks["residual_evaluable_strata_129"] = int((diagnostics["status"] == "evaluable").sum()) == 129
    checks["tcga_5440_and_8_subtype_tests"] = len(tcga_corr) == 5440 and len(tcga_tests) == 8
    checks["tcga_68immune_p_family_single"] = sorted(tcga_corr["p_family"].dropna().unique().tolist()) == [
        "tcga_xena_68immune_within_cancer_correlations_v0_5"
    ]
    checks["tcga_68immune_finite_p_and_q_5440"] = int(np.isfinite(tcga_p).sum()) == 5440 and int(np.isfinite(tcga_q).sum()) == 5440
    checks["tcga_bh_max_abs_diff"] = tcga_q_max_abs_diff
    checks["tcga_bh_matches_full_5440_family"] = tcga_q_mask.sum() == 5440 and tcga_q_max_abs_diff <= 1e-12
    checks["tcga_subtype_separate_8_test_family"] = len(tcga_tests) == 8 and pd.to_numeric(tcga_tests["kruskal_fdr_bh"], errors="coerce").notna().sum() == 8
    checks["apollo_counts"] = apollo["evidence_status"].value_counts().to_dict()
    checks["apollo_expected_counts"] = checks["apollo_counts"] == {"supported": 11, "not_evaluable": 7, "null_result": 5, "opposite": 1}
    checks["apollo_score_duplicate_inference_keys_zero"] = int(apollo_scores.duplicated(apollo_score_key, keep=False).sum()) == 0
    checks["apollo_registry_duplicate_program_keys_zero"] = int(apollo.duplicated(apollo_registry_key, keep=False).sum()) == 0
    checks["apollo_fdr_family_rows"] = apollo.groupby("fdr_family").size().to_dict()
    checks["apollo_fdr_families_3x8"] = checks["apollo_fdr_family_rows"] == {
        "apollo_APOLLO-LUAD_phosphoproteomics": 8,
        "apollo_APOLLO-LUAD_proteomics": 8,
        "apollo_APOLLO-OV_proteomics": 8,
    }
    checks["apollo_not_evaluable_rule_counts"] = (
        int(apollo["status"].eq("not_evaluable").sum()) == 7
        and set(apollo.loc[apollo["status"].eq("not_evaluable"), "reason"].dropna().unique().tolist()) == {"missing_program_score_or_n_below_10"}
    )
    checks["apollo_ci_max_abs_diff"] = apollo_ci_max_abs_diff
    checks["apollo_ci_matches_fisher_z"] = len(apollo_evaluable) == 17 and apollo_ci_max_abs_diff <= 1e-12
    checks["scrna_support_count"] = len(scrna_support) == 116
    checks["spatial_sample_level_rows"] = len(spatial_pairwise) == 24 and set(spatial_pairwise["paired_sample_n"].astype(int)) == {4}
    checks["table2_has_16_rows"] = len(table2) == 16
    checks["table2_zero_evaluable_rows"] = int((table2["evaluable_cancer_n"] == 0).sum()) == 2
    checks["table2_zero_rows_are_expected"] = set(table2.loc[table2["evaluable_cancer_n"] == 0, "program_label"].tolist()) == {"MHC-II", "Cytolytic"}
    checks["numeric_map_has_required_rows"] = len(numeric_map) >= 70
    checks["numeric_map_sources_exist"] = all((root / str(path)).exists() for path in numeric_map["source_file"].dropna().unique() if str(path).endswith((".csv", ".json")))
    checks["text_contains_key_numeric_tokens"] = contains_all(
        text,
        [
            "160-row primary registry",
            "77 evaluable RNA-protein rows",
            "52 evaluable RNA-phosphoprotein rows",
            "finite tested denominator was 516",
            "5440 within-cancer correlations",
            "Across 24 registered cohort-layer-program rows",
            "345 source-specific donor units",
            "paired sample n=4",
        ],
    )
    checks["text_contains_freeze_p1_closure_tokens"] = contains_all(
        text,
        [
            "all 5440 finite Spearman P values",
            "no within-cancer, within-program, or within-immune-signature partition",
            "single 5440-test TCGA 68-signature BH family",
            "public `sample_id` within each cohort-layer score table",
            "17 evaluable and 7 non-evaluable rows",
            "APOLLO 95% Spearman confidence intervals used a Fisher z interval",
        ],
    )
    checks["max_influence_text_matches_source"] = (
        "GBM Checkpoint RNA-phosphoprotein" in text
        and str(max_influence["most_influential_case_id"]) in text
        and "0.118" in text
    )
    checks["source_manifest_timing_disclosed"] = "finalized at 2026-08-12 02:19:28 UTC, after result generation" in text

    for key, value in checks.items():
        if isinstance(value, bool) and not value:
            failures.append(key)
    if checks.get("residual_context_by_variable") != {
        "cibersort_myeloid": 129,
        "cibersort_t_cell": 129,
        "estimate_immune_score": 129,
        "estimate_tumor_purity": 129,
    }:
        failures.append("residual_context_by_variable")

    gate = {
        "project": PROJECT,
        "generated_at_utc": now_utc(),
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "checks": checks,
        "audited_files": [
            METHODS_REL,
            "current_p1_closure/NUMERIC_MAP.csv",
            "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_CURRENT_P1.csv",
            "results/cptac/cptac_primary_concordance_registry_v2.csv",
            "results/cptac/cptac_residual_context_registry_v2.csv",
            "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv",
            "results/tcga/tcga_xena_immune_subtype_tests_v0_5.csv",
            "results/apollo/apollo_signature_scores_v0_7.csv",
            "results/apollo/apollo_concordance_registry_v2.csv",
            "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv",
            "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv",
        ],
    }
    write_json(root / "current_p1_closure/NUMERIC_GATE_CURRENT_P1.json", gate)

    closure["numeric_gate"] = gate
    if failures:
        closure["status"] = "BLOCKED_REQUIRED_SOURCE"
        closure["unresolved_issues"] = failures
    else:
        closure["status"] = "READY_FOR_EXTERNAL_R2"
        closure["unresolved_issues"] = []
    if "external_freeze_p1_closure" in closure:
        closure["external_freeze_p1_closure"]["status"] = (
            "READY_FOR_EXTERNAL_R2" if not failures else "BLOCKED_REQUIRED_SOURCE"
        )
    files = set(closure.get("files_changed", []))
    files.add("current_p1_closure/NUMERIC_GATE_CURRENT_P1.json")
    files.add("current_p1_closure/build_current_p1_outputs.py")
    files.add("current_p1_closure/run_current_p1_numeric_gate.py")
    files.add("current_p1_freeze_remediation/FREEZE_REMEDIATION.json")
    files.add("current_p1_freeze_remediation/FREEZE_REMEDIATION.md")
    closure["files_changed"] = sorted(files)
    write_freeze_remediation(root, gate, closure, failures)
    write_json(closure_path, closure)
    update_closure_md(root, gate["status"], failures)

    print(json.dumps({"project": PROJECT, "numeric_gate": gate["status"], "failures": failures}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
