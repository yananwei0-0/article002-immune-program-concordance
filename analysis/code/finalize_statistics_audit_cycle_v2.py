#!/usr/bin/env python3
"""Finalize the BIO-05 statistics-audit cycle without touching other projects."""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT = "BIO-05"
PROTOCOL_SHA256 = "a3bc200faf24eefb1b89be08cb4b6804b81423f50c3337a430434e6015143fc3"
ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
CYCLE_ROOT = Path(expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/sci27_dual_pro_statistics_cycle_20260812"))
REVISION_ROOT = CYCLE_ROOT / "cli_revisions" / PROJECT
FINAL_PATH = Path(
    expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/manuscript27_v2_scientific_rebuild_20260811/"
    "scientific_reaudit_waves/wave2/BIO-05/final.json")
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def output_manifest(root: Path) -> list[dict[str, Any]]:
    records = []
    for folder in ["code", "results", "tables", "figures", "manuscript", "local_audit"]:
        base = root / folder
        if not base.exists():
            continue
        for path in sorted(x for x in base.rglob("*") if x.is_file() and "__pycache__" not in x.parts):
            records.append(
                {
                    "relative_path": str(path.relative_to(root)),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    return records


def main() -> int:
    gate = json.loads((ANALYSIS_ROOT / "AUDIT_GATE.json").read_text())
    validation = json.loads((ANALYSIS_ROOT / "VALIDATION_REPORT.json").read_text())
    results = json.loads((ANALYSIS_ROOT / "RESULTS_SUMMARY.json").read_text())
    rerun = json.loads((ANALYSIS_ROOT / "local_audit/RERUN_CONSISTENCY.json").read_text())
    if gate["gate_status"] != "READY" or gate["open_p0"] or gate["open_p1"]:
        raise RuntimeError("BIO-05 independent gate is not READY with P0/P1 zero")
    if validation["status"] != "PASS" or results["status"] != "V2_ANALYSIS_COMPLETE":
        raise RuntimeError("BIO-05 authority files are not complete")
    if not rerun["all_core_outputs_byte_identical"] or rerun["byte_identical_count"] != 7:
        raise RuntimeError("BIO-05 revised pipeline rerun is not 7/7 byte-identical")

    run_manifest_path = ANALYSIS_ROOT / "RUN_MANIFEST.json"
    run_manifest = json.loads(run_manifest_path.read_text())
    run_manifest["completed_at_utc"] = now_utc()
    run_manifest["output_manifest"] = output_manifest(ANALYSIS_ROOT)
    run_manifest["validation_status"] = "PASS"
    run_manifest["statistics_audit_cycle"] = {
        "cycle": "STAT-R1",
        "status": "P0_P1_CLOSED",
        "audit_source": str(CYCLE_ROOT / "initial_audit_responses/BIO-05_INITIAL_STAT_AUDIT.md"),
        "revision_record": str(REVISION_ROOT / "finding_adjudication.json"),
        "rerun_core_outputs_byte_identical": "7/7",
    }
    write_json(run_manifest_path, run_manifest)

    authority_paths = {
        "SOURCE_MANIFEST.json": ANALYSIS_ROOT / "SOURCE_MANIFEST.json",
        "RUN_MANIFEST.json": run_manifest_path,
        "RESULTS_SUMMARY.json": ANALYSIS_ROOT / "RESULTS_SUMMARY.json",
        "VALIDATION_REPORT.json": ANALYSIS_ROOT / "VALIDATION_REPORT.json",
        "AUDIT_GATE.json": ANALYSIS_ROOT / "AUDIT_GATE.json",
    }
    authority_hashes = {name: sha256(path) for name, path in authority_paths.items()}

    adjudication = {
        "project": PROJECT,
        "round": "STAT-R1",
        "generated_at_utc": now_utc(),
        "protocol_sha256": PROTOCOL_SHA256,
        "analysis_plan_deviation": False,
        "status": "P0_P1_CLOSED",
        "open_p0": 0,
        "open_p1": 0,
        "findings": [
            {
                "finding_id": "BIO05-STAT-R1-P1-01",
                "severity": "P1",
                "decision": "ACCEPT_AND_REANALYZE",
                "status": "CLOSED",
                "issue": "Fisher-z interval was not a calibrated primary confidence interval for Spearman rho; small-n/tie P-value calibration required audit.",
                "actions": [
                    "Preserved all locked scores, estimands, strata, thresholds, and the 160-row multiplicity family.",
                    "Recomputed all 129 evaluable primary 95% intervals using 20,000 deterministic paired-case BCa bootstrap resamples.",
                    "Verified that no evaluable primary row contained score ties and that only two rows had n<=30.",
                    "Recomputed those two P values with 999,999 deterministic two-sided Monte Carlo permutations and plus-one correction; retained Monte Carlo standard errors.",
                    "Recomputed BH q values once across the unchanged finite primary family and propagated results to CSV, workbook, Methods/Results, authority files, and gate.",
                ],
                "key_results": {
                    "primary_evaluable_rows": 129,
                    "bca_interval_rows": 129,
                    "tied_evaluable_rows": 0,
                    "permutation_calibrated_rows": 2,
                    "coad_checkpoint_phosphoproteomics": {
                        "n": 25,
                        "rho": 0.446154,
                        "bca_95_ci": [-0.021654, 0.772278],
                        "asymptotic_p": 0.025385,
                        "permutation_p": 0.026419,
                        "bh_q": 0.027048,
                    },
                    "gbm_checkpoint_phosphoproteomics": {
                        "n": 26,
                        "rho": 0.506325,
                        "bca_95_ci": [0.061938, 0.784390],
                        "asymptotic_p": 0.008307,
                        "permutation_p": 0.009007,
                        "bh_q": 0.009603,
                    },
                    "primary_bh_q_lt_0_05_rows": 128,
                    "rerun_byte_identical_core_files": "7/7",
                },
                "evidence": [
                    "results/cptac/cptac_primary_concordance_registry_v2.csv",
                    "results/cptac/cptac_internal_validation_v2.json",
                    "local_audit/RERUN_CONSISTENCY.json",
                    "AUDIT_GATE.json",
                ],
            },
            {
                "finding_id": "BIO05-STAT-R1-P1-02",
                "severity": "P1",
                "decision": "ACCEPT_AND_RESTRICT_CLAIM",
                "status": "CLOSED",
                "issue": "Source-prefixed donor labels did not establish cross-source biological independence for pooled global single-cell inference.",
                "actions": [
                    "Audited all 10 sources by accession, publication/collection identity, cancer, donor namespace, and exact raw donor-label overlap.",
                    "Confirmed zero exact raw donor-label overlaps, while explicitly recording that a cross-repository biological-participant crosswalk is unavailable.",
                    "Verified that the unusual six-label ccRCC donor field maps one-to-one to six HCA donor UUIDs in the retained H5AD.",
                    "Removed all pooled global donor-unit P values, q values, and support claims rather than assuming independence.",
                    "Retained 24 pooled rows only as descriptive effect summaries and retained source-cancer donor-paired Wilcoxon families as bounded supporting evidence.",
                ],
                "key_results": {
                    "audited_sources": 10,
                    "source_specific_donor_units": 345,
                    "exact_raw_donor_labels_reused_across_sources": 0,
                    "globally_unique_biological_donor_crosswalk_available": False,
                    "global_descriptive_rows": 24,
                    "global_finite_p_rows": 0,
                    "global_finite_q_rows": 0,
                    "global_support_rows": 0,
                    "within_source_evaluable_rows": 256,
                    "within_source_support_rows": 116,
                },
                "evidence": [
                    "results/scrna_combined/public_scrna_source_donor_provenance_audit_v2.csv",
                    "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv",
                    "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv",
                    "results/support_layers_summary_v2.json",
                    "AUDIT_GATE.json",
                ],
            },
        ],
        "p2_adjudication": [
            {
                "finding": "Score-composition/missingness robustness",
                "decision": "OPTIONAL_NOT_FREEZE_BLOCKING",
                "reason": "The registered three-gene primary rule, one-gene sensitivity, feature-row sensitivity, gene recovery, nonmissing counts, and primary denominator registry remain available; no P0/P1 contradiction was identified.",
            },
            {
                "finding": "Cross-cancer heterogeneity presentation",
                "decision": "ALREADY_BOUNDED_NOT_FREEZE_BLOCKING",
                "reason": "Cancer-specific rho remains primary; median rho, positive-cancer fraction, Q, and I2 are explicitly descriptive and no homogeneous pooled effect is claimed.",
            },
        ],
        "authority_hashes": authority_hashes,
        "gate": {"status": "READY", "open_p0": 0, "open_p1": 0},
    }
    write_json(REVISION_ROOT / "finding_adjudication.json", adjudication)

    summary = f"""# BIO-05 STAT-R1 local revision summary

- Status: **P0/P1 CLOSED**
- Machine/independent gate: **READY**
- Open P0/P1: **0/0**
- Locked protocol SHA-256: `{PROTOCOL_SHA256}`; unchanged
- Analysis-plan deviation: **false**

## P1-1: Spearman uncertainty calibration

All 129 evaluable primary rows now use deterministic paired-case BCa bootstrap 95% intervals with 20,000 resamples. The two rows with n<=30 use 999,999-permutation two-sided P values; no primary row had score ties. The primary BH family remained the same 160 registered rows, and 128 rows remained q<0.05. The revised pipeline reproduced 7/7 core files byte-for-byte.

## P1-2: cross-source donor independence

Ten sources and 345 source-specific donor units were audited. No exact donor label was reused, but no cross-repository participant crosswalk could establish global biological independence. Pooled global P/q values and the former 21/24 support claim were therefore removed. The 24 pooled rows are descriptive only; 256 source-cancer rows remain evaluable and 116 meet the separately controlled source-bounded support window.

## Final boundary

The evidence supports cancer-specific CPTAC RNA-protein/phosphoprotein concordance and bounded source-specific public triangulation. It does not establish a homogeneous pan-cancer effect, independent validation of the full framework, malignant-cell mechanism, prognosis, prediction, treatment response, or clinical utility.
"""
    (REVISION_ROOT / "summary.md").write_text(summary, encoding="utf-8")

    commands = """# BIO-05 STAT-R1 command ledger

1. Read initial_audit_responses/BIO-05_INITIAL_STAT_AUDIT.md and all BIO-05 protocol/authority/code/results/gate assets.
2. python3 -m py_compile code/run_cptac_v2.py code/postprocess_support_v2.py code/build_final_assets_v2.py code/run_independent_audit_v2.py
3. python3 code/run_cptac_v2.py --analysis-root <BIO-05> --data-root <BIO-05 data root> --protocol <locked protocol>
4. python3 code/postprocess_support_v2.py --analysis-root <BIO-05>
5. python3 code/build_final_assets_v2.py --analysis-root <BIO-05> --protocol <locked protocol>
6. python3 code/run_independent_audit_v2.py --analysis-root <BIO-05> --data-root <BIO-05 data root> --protocol <locked protocol>
7. python3 code/check_rerun_consistency_v2.py --analysis-root <BIO-05> --data-root <BIO-05 data root> --protocol <locked protocol>
8. Re-ran independent audit after revised rerun evidence; gate READY, P0=0, P1=0.
9. Inspected the supplementary workbook and authority propagation; no stale global single-cell inferential claim remained.
10. python3 code/finalize_statistics_audit_cycle_v2.py

No model, score, threshold, hypothesis, multiplicity family, or sample was changed to obtain significance.
"""
    (REVISION_ROOT / "commands.log").write_text(commands, encoding="utf-8")

    modified = """code/run_cptac_v2.py
code/postprocess_support_v2.py
code/build_final_assets_v2.py
code/run_independent_audit_v2.py
code/finalize_statistics_audit_cycle_v2.py
results/cptac/cptac_primary_concordance_registry_v2.csv
results/cptac/cptac_feature_row_sensitivity_registry_v2.csv
results/cptac/cptac_threshold1_sensitivity_registry_v2.csv
results/cptac/cptac_internal_validation_v2.json
results/cptac/cptac_results_summary_v2.json
results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv
results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv
results/scrna_combined/public_scrna_source_donor_provenance_audit_v2.csv
results/support_layers_summary_v2.json
tables/SUPPLEMENT_TABLE_S1_PRIMARY_REGISTRY_V2.csv
tables/SUPPLEMENT_TABLE_S2_SCORING_SENSITIVITY_V2.csv
tables/BIO-05_SUPPLEMENTARY_DATA_V2.xlsx
manuscript/BIO-05_METHODS_RESULTS_V2.md
manuscript/BIO-05_FIGURE_LEGENDS_V2.md
SOURCE_MANIFEST.json
RUN_MANIFEST.json
RESULTS_SUMMARY.json
VALIDATION_REPORT.json
AUDIT_GATE.json
local_audit/RERUN_CONSISTENCY.json
local_audit/BIO-05_INDEPENDENT_AUDIT.md
logs/COMMANDS_V2.log
scientific_reaudit_waves/wave2/BIO-05/final.json
outputs/sci27_dual_pro_statistics_cycle_20260812/cli_revisions/BIO-05/finding_adjudication.json
outputs/sci27_dual_pro_statistics_cycle_20260812/cli_revisions/BIO-05/commands.log
outputs/sci27_dual_pro_statistics_cycle_20260812/cli_revisions/BIO-05/modified_files.txt
outputs/sci27_dual_pro_statistics_cycle_20260812/cli_revisions/BIO-05/summary.md
"""
    (REVISION_ROOT / "modified_files.txt").write_text(modified, encoding="utf-8")

    final = {
        "project": PROJECT,
        "verdict": "READY",
        "open_p0": 0,
        "open_p1": 0,
        "protocol_sha256": PROTOCOL_SHA256,
        "analysis_plan_deviation": False,
        "statistics_audit_round": "STAT-R1",
        "key_results": [
            "CPTAC included 1,023 independent cases; 129/160 primary rows were evaluable and 128 had BH q<0.05.",
            "All 129 primary 95% intervals use 20,000-resample paired-case BCa bootstrap; 2/2 n<=30 rows use 999,999-permutation calibrated P values.",
            "The revised pipeline reproduced 7/7 core files byte-for-byte; independent rho and BH maximum absolute differences were 1.11e-16 and 6.05e-17.",
            "Across 10 single-cell/3CA sources, 345 source-specific donor units were observed, but global donor independence was not established; pooled global P/q values and support claims were removed.",
            "Within-source single-cell evidence retained 256 evaluable rows and 116 source-bounded support rows; all opposite, null, and non-evaluable evidence remained visible.",
        ],
        "authority_hashes": authority_hashes,
        "audit_gate": "AUDIT_GATE.json: READY; real-source execution complete; open P0=0, open P1=0.",
        "revision_record": str(REVISION_ROOT / "finding_adjudication.json"),
        "summary": "BIO-05 STAT-R1 closed both P1 findings by inferential recalibration and claim restriction without changing the locked scientific estimand or targeting significance.",
    }
    write_json(FINAL_PATH, final)
    print(json.dumps({"project": PROJECT, "status": "READY", "open_p0": 0, "open_p1": 0}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
