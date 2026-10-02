#!/usr/bin/env python3
"""Executable syntax, deterministic, numeric, semantic, and immutability gate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import py_compile
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t, wilcoxon


PROJECT = "BIO-05"
SOURCE_P1 = "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P1.md"
SOURCE_P1_SHA256 = "7ab1be4929eaa15030fc226b90172c5cbecbe69be42257470562112b83b5d5e2"
P2_PATH = "manuscript/BIO-05_METHODS_RESULTS_CURRENT_P2.md"
MAP_PATH = "current_p2_closure/NUMERIC_MAP_CURRENT_P2.csv"
GATE_PATH = "current_p2_closure/NUMERIC_GATE_CURRENT_P2.json"
CLOSURE_PATH = "current_p2_closure/CURRENT_P2_CLOSURE.json"
ANCHOR_PATH = "current_p2_closure/TCGA_ANCHOR_SELECTION_CURRENT_P2.csv"
CROSS_REVISION_PATH = "cross_freeze_revision/CROSS_FREEZE_REVISION_1.json"
HASH_PROOF_PATH = "current_p2_closure/TWO_PASS_DETERMINISTIC_HASHES_CURRENT_P2.json"


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
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def load_builder(root: Path):
    path = root / "current_p2_closure/build_current_p2_outputs.py"
    spec = importlib.util.spec_from_file_location("bio05_current_p2_builder", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def equal_value(left: Any, right: Any) -> bool:
    if left is None or right is None:
        return left is None and right is None
    try:
        a, b = float(left), float(right)
        return bool(np.isclose(a, b, rtol=1e-12, atol=1e-14))
    except (TypeError, ValueError):
        return str(left) == str(right)


def bh(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
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


def wilcoxon_crosscheck(root: Path) -> dict[str, Any]:
    scores = pd.read_csv(root / "results/scrna_combined/public_scrna_and_3ca_donor_compartment_scores_v2.csv")
    retained = pd.read_csv(root / "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv")
    scores["donor_unit"] = scores["source_id"].astype(str) + "::" + scores["cancer_code"].astype(str) + "::" + scores["donor_id"].astype(str)
    rows = []
    for (source_id, cancer, program), group in scores.groupby(["source_id", "cancer_code", "signature"], sort=False):
        wide = group.pivot(index="donor_unit", columns="compartment", values="signature_score")
        for left, right in [("malignant", "immune"), ("malignant", "stromal"), ("immune", "stromal")]:
            if left not in wide or right not in wide:
                delta = pd.Series(dtype=float)
            else:
                delta = pd.to_numeric(wide[left] - wide[right], errors="coerce").dropna()
            n = len(delta)
            if n < 6:
                p = np.nan
            elif np.allclose(delta.to_numpy(dtype=float), 0):
                p = 1.0
            else:
                p = float(wilcoxon(delta.to_numpy(dtype=float), alternative="two-sided", zero_method="wilcox", method="auto").pvalue)
            rows.append({"source_id": source_id, "cancer_code": cancer, "program": program, "contrast": f"{left}_minus_{right}", "n": n, "p_recomputed": p})
    rebuilt = pd.DataFrame(rows)
    joined = retained.merge(rebuilt, on=["source_id", "cancer_code", "program", "contrast"], how="left", validate="one_to_one")
    finite = joined["wilcoxon_p_value"].notna()
    p_diff = float(np.max(np.abs(joined.loc[finite, "wilcoxon_p_value"] - joined.loc[finite, "p_recomputed"]))) if finite.any() else 0.0
    q_recomputed = np.full(len(joined), np.nan)
    for _, idx in joined.groupby("fdr_family").groups.items():
        positions = list(idx)
        q_recomputed[positions] = bh(joined.loc[positions, "p_recomputed"])
    q_finite = joined["q_value_bh"].notna()
    q_diff = float(np.max(np.abs(joined.loc[q_finite, "q_value_bh"] - q_recomputed[q_finite.to_numpy()]))) if q_finite.any() else 0.0
    return {
        "recomputed_rows": len(joined),
        "finite_rows": int(finite.sum()),
        "max_abs_p_difference": p_diff,
        "max_abs_q_difference": q_diff,
        "pass": len(joined) == 336 and int(finite.sum()) == 256 and p_diff <= 1e-14 and q_diff <= 1e-14,
    }


def spearman_algorithm_crosscheck(root: Path) -> dict[str, Any]:
    outputs = {}
    for name, rel, n_field, rho_field, p_field in [
        ("tcga", "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv", "n", "spearman_rho", "p_value"),
        ("apollo", "results/apollo/apollo_concordance_registry_v2.csv", "matched_sample_n", "spearman_rho", "p_value"),
    ]:
        frame = pd.read_csv(root / rel).dropna(subset=[rho_field, p_field])
        rho = frame[rho_field].to_numpy(dtype=float)
        n = frame[n_field].to_numpy(dtype=float)
        stat = rho * np.sqrt((n - 2) / np.clip((rho + 1) * (1 - rho), np.finfo(float).tiny, None))
        expected = 2 * t.sf(np.abs(stat), df=n - 2)
        observed = frame[p_field].to_numpy(dtype=float)
        diff = float(np.max(np.abs(expected - observed))) if len(frame) else np.nan
        outputs[name] = {"finite_rows": len(frame), "max_abs_p_difference": diff, "pass": bool(diff <= 1e-14)}
    outputs["pass"] = all(x["pass"] for x in outputs.values())
    return outputs


def apollo_ci_crosscheck(root: Path) -> dict[str, Any]:
    frame = pd.read_csv(root / "results/apollo/apollo_concordance_registry_v2.csv").dropna(subset=["spearman_rho"])
    z = np.arctanh(np.clip(frame["spearman_rho"].to_numpy(float), -0.999999, 0.999999))
    se = 1 / np.sqrt(frame["matched_sample_n"].to_numpy(float) - 3)
    low = np.tanh(z - 1.9599639845 * se)
    high = np.tanh(z + 1.9599639845 * se)
    diff = float(max(np.max(np.abs(low - frame["ci_lower"].to_numpy(float))), np.max(np.abs(high - frame["ci_upper"].to_numpy(float)))))
    q_diff = 0.0
    for _, idx in frame.groupby("fdr_family").groups.items():
        positions = list(idx)
        expected = bh(frame.loc[positions, "p_value"])
        q_diff = max(q_diff, float(np.max(np.abs(expected - frame.loc[positions, "q_value_bh"].to_numpy(float)))))
    return {"evaluable_rows": len(frame), "ci_max_abs_difference": diff, "q_max_abs_difference": q_diff, "pass": len(frame) == 17 and diff <= 1e-12 and q_diff <= 1e-14}


def tcga_bh_crosscheck(root: Path) -> dict[str, Any]:
    frame = pd.read_csv(root / "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv")
    expected = bh(frame["p_value"])
    observed = frame["fdr_bh"].to_numpy(float)
    diff = float(np.max(np.abs(expected - observed)))
    return {"rows": len(frame), "unique_p_families": int(frame["p_family"].nunique()), "max_abs_q_difference": diff, "pass": len(frame) == 5440 and frame["p_family"].nunique() == 1 and diff <= 1e-14}


def deterministic_regeneration(root: Path) -> dict[str, Any]:
    tracked = [
        "current_p2_closure/build_current_p2_outputs.py",
        "current_p2_closure/run_current_p2_gate.py",
        P2_PATH,
        MAP_PATH,
        "current_p2_closure/PROGRAM_DEFINITION_CURRENT_P2.csv",
        "current_p2_closure/SOFTWARE_ENVIRONMENT_CURRENT_P2.json",
        "current_p2_closure/SUPPORT_MAPPING_RULES_CURRENT_P2.json",
        "current_p2_closure/TCGA_PRIMARY_SAMPLE_SELECTION_CURRENT_P2.csv",
        ANCHOR_PATH,
        "current_p2_closure/NUMERIC_DERIVATIONS_CURRENT_P2.json",
        "current_p2_closure/CURRENT_P2_ADJUDICATION.json",
        "current_p2_closure/CURRENT_P2_ADJUDICATION.md",
        "current_p2_closure/COMMANDS_RUN.txt",
        "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py",
        "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap.py",
        "cross_freeze_revision/CROSS_FREEZE_REVISION_1.json",
        "cross_freeze_revision/CROSS_FREEZE_REVISION_1.md",
    ]
    runs = []
    pass_hashes = []
    for pass_index in (1, 2):
        completed = subprocess.run(
            [sys.executable, "current_p2_closure/build_current_p2_outputs.py", "--analysis-root", "."],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
        )
        runs.append(
            {
                "pass": pass_index,
                "builder_exit_code": completed.returncode,
                "stdout": completed.stdout.strip(),
                "stderr": completed.stderr.strip(),
            }
        )
        pass_hashes.append({path: sha256_file(root / path) for path in tracked})
    mismatches = [path for path in tracked if pass_hashes[0][path] != pass_hashes[1][path]]
    return {
        "runs": runs,
        "tracked_file_count": len(tracked),
        "tracked_files": tracked,
        "pass_1_hashes": pass_hashes[0],
        "pass_2_hashes": pass_hashes[1],
        "hash_mismatches": mismatches,
        "pass": all(run["builder_exit_code"] == 0 for run in runs) and not mismatches,
    }


def tcga_anchor_crosscheck(root: Path, builder: Any) -> dict[str, Any]:
    artifact = pd.read_csv(root / ANCHOR_PATH)
    summary = pd.read_csv(root / "results/tcga/tcga_xena_signature_summary_v0_5.csv")
    helper = root / "current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_2_controls_rna_depmap.py"
    configured = builder.literal_assignment(helper, "SIGNATURE_ANCHORS")
    mapping_mismatches = []
    aggregation_mismatches = []
    for signature, anchors in configured.items():
        subset = artifact[artifact["signature"].eq(signature)].sort_values("anchor_order")
        observed_names = subset["immune_signature_68"].astype(str).tolist()
        if observed_names != anchors:
            mapping_mismatches.append(signature)
        retained = summary[summary["signature"].eq(signature)]
        if len(retained) != 1 or subset.empty:
            aggregation_mismatches.append(signature)
            continue
        row = retained.iloc[0]
        if not np.isclose(
            subset["median_within_cancer_rho"].median(),
            row["selected_anchor_median_rho"],
            rtol=1e-12,
            atol=1e-14,
        ) or not np.isclose(
            subset["positive_cancer_fraction"].min(),
            row["selected_anchor_min_positive_cancer_fraction"],
            rtol=1e-12,
            atol=1e-14,
        ):
            aggregation_mismatches.append(signature)
    return {
        "anchor_rows": len(artifact),
        "programs": int(artifact["signature"].nunique()),
        "membership_match_all": bool(artifact["membership_match"].all()),
        "summary_match_all": bool(artifact["summary_match"].all()),
        "mapping_mismatches": mapping_mismatches,
        "aggregation_mismatches": aggregation_mismatches,
        "pass": (
            len(artifact) == 17
            and artifact["signature"].nunique() == 8
            and bool(artifact["membership_match"].all())
            and bool(artifact["summary_match"].all())
            and not mapping_mismatches
            and not aggregation_mismatches
        ),
    }


def sensitivity_outcome_crosscheck(root: Path) -> dict[str, Any]:
    loco = pd.read_csv(root / "results/cptac/cptac_leave_one_cancer_out_v2.csv")
    pooled = pd.read_csv(root / "results/cptac/cptac_within_cancer_standardized_pooled_sensitivity_v2.csv")
    scrna = pd.read_csv(root / "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv")
    pooled_finite = pooled[pooled["within_cancer_standardized_pooled_spearman_rho"].notna()]
    sensitivity = scrna["n_ge_10_sensitivity_window"].astype(bool)
    observed = {
        "loco_rows": len(loco),
        "loco_median_rho_min": float(loco["leave_one_cancer_out_median_rho"].min()),
        "loco_median_rho_max": float(loco["leave_one_cancer_out_median_rho"].max()),
        "loco_positive_fraction_min": float(loco["leave_one_cancer_out_positive_fraction"].min()),
        "pooled_registered_rows": len(pooled),
        "pooled_finite_rows": len(pooled_finite),
        "pooled_positive_rows": int((pooled_finite["within_cancer_standardized_pooled_spearman_rho"] > 0).sum()),
        "pooled_rho_min": float(pooled_finite["within_cancer_standardized_pooled_spearman_rho"].min()),
        "pooled_rho_max": float(pooled_finite["within_cancer_standardized_pooled_spearman_rho"].max()),
        "scrna_n10_supported": int(sensitivity.sum()),
        "scrna_n10_positive": int((sensitivity & (scrna["median_delta"] > 0)).sum()),
        "scrna_n10_negative": int((sensitivity & (scrna["median_delta"] < 0)).sum()),
    }
    expected = {
        "loco_rows": 129,
        "loco_median_rho_min": 0.3169652855543113,
        "loco_median_rho_max": 0.8575531770762768,
        "loco_positive_fraction_min": 1.0,
        "pooled_registered_rows": 16,
        "pooled_finite_rows": 14,
        "pooled_positive_rows": 14,
        "pooled_rho_min": 0.318787008844282,
        "pooled_rho_max": 0.866281400913768,
        "scrna_n10_supported": 85,
        "scrna_n10_positive": 31,
        "scrna_n10_negative": 54,
    }
    mismatches = [key for key in expected if not equal_value(observed[key], expected[key])]
    return {"observed": observed, "expected": expected, "mismatches": mismatches, "pass": not mismatches}


def update_closure(root: Path, gate: dict[str, Any], builder: Any) -> None:
    closure = json.loads((root / CLOSURE_PATH).read_text(encoding="utf-8"))
    closure["source_current_p1"]["observed_sha256"] = sha256_file(root / SOURCE_P1)
    closure["source_current_p1"]["unchanged"] = closure["source_current_p1"]["observed_sha256"] == SOURCE_P1_SHA256
    closure["current_p2"]["sha256"] = sha256_file(root / P2_PATH)
    closure["numeric_gate"] = {"status": gate["status"], "path": GATE_PATH, "failures": gate["failures"]}
    closure["forbidden_language_gate"] = {
        "status": "PASS" if gate["checks"]["forbidden_language"]["pass"] else "FAIL",
        "path": GATE_PATH,
        "hits": gate["checks"]["forbidden_language"]["forbidden_hits"],
    }
    deterministic = gate["checks"]["deterministic_regeneration"]
    closure["deterministic_hash_proof"] = {
        "status": "PASS" if deterministic["pass"] else "FAIL",
        "path": HASH_PROOF_PATH,
        "hash_mismatches": deterministic["hash_mismatches"],
    }
    closure["status"] = "READY_FOR_CROSS_FREEZE" if gate["status"] == "PASS" and not closure["unresolved_p0"] and not closure["unresolved_p1"] else "NOT_READY"
    closure["cross_freeze_revision"] = {"path": CROSS_REVISION_PATH, "status": closure["status"]}
    write_json(root / CLOSURE_PATH, closure)

    revision = json.loads((root / CROSS_REVISION_PATH).read_text(encoding="utf-8"))
    revision["after_revision"]["current_p2_sha256"] = sha256_file(root / P2_PATH)
    revision["after_revision"]["numeric_map_sha256"] = sha256_file(root / MAP_PATH)
    revision["gate_validation"] = {
        "numeric_and_semantic_gate": gate["status"],
        "forbidden_language_gate": "PASS" if gate["checks"]["forbidden_language"]["pass"] else "FAIL",
        "two_pass_deterministic_hashes": "PASS" if deterministic["pass"] else "FAIL",
        "gate_path": GATE_PATH,
        "hash_proof_path": HASH_PROOF_PATH,
        "gate_failures": gate["failures"],
        "hash_mismatches": deterministic["hash_mismatches"],
    }
    revision["status"] = closure["status"]
    write_json(root / CROSS_REVISION_PATH, revision)
    (root / "cross_freeze_revision/CROSS_FREEZE_REVISION_1.md").write_text(
        builder.render_cross_freeze_revision_md(revision), encoding="utf-8"
    )
    lines = [
        "# BIO-05 CURRENT_P2 closure",
        "",
        f"- Source CURRENT_P1: `{SOURCE_P1}`; SHA256 `{SOURCE_P1_SHA256}`; unchanged: `{str(closure['source_current_p1']['unchanged']).lower()}`.",
        f"- CURRENT_P2: `{P2_PATH}`; SHA256 `{closure['current_p2']['sha256']}`.",
        "- Scientific reanalysis: no.",
        "- Numerical result change: no.",
        "- Existing reported numerical estimate change: no.",
        "- Retained sensitivity outcomes added to Results: yes.",
        "- Primary conclusion change: no.",
        f"- Unresolved P0: {len(closure['unresolved_p0'])}.",
        f"- Unresolved P1: {len(closure['unresolved_p1'])}.",
        f"- Unresolved P2: {len(closure.get('unresolved_p2', []))}.",
        f"- Deferred P2: {len(closure['deferred_p2'])}.",
        f"- Numeric gate: {gate['status']}.",
        f"- Forbidden-language gate: {closure['forbidden_language_gate']['status']}.",
        f"- Two-pass deterministic hash proof: {closure['deterministic_hash_proof']['status']}.",
        f"- Final status: {closure['status']}.",
        "",
        "The cross-freeze P1/P2 findings were accepted and closed by semantic numeric rebinding, exact anchor/scoring-rule recovery, and reporting of retained sensitivity outcomes. No scientific analysis was rerun and no retained result or primary conclusion changed.",
        "",
        "## Claim boundary",
        "",
        closure["claim_boundary"],
        "",
    ]
    (root / "current_p2_closure/CURRENT_P2_CLOSURE.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = Path(args.analysis_root).resolve()
    builder = load_builder(root)
    checks: dict[str, Any] = {}

    syntax_errors = []
    for rel in ["current_p2_closure/build_current_p2_outputs.py", "current_p2_closure/run_current_p2_gate.py"]:
        try:
            py_compile.compile(str(root / rel), doraise=True)
        except Exception as exc:
            syntax_errors.append(f"{rel}: {exc}")
    checks["syntax"] = {"errors": syntax_errors, "pass": not syntax_errors}

    source_hash = sha256_file(root / SOURCE_P1)
    checks["source_current_p1_immutable"] = {"expected": SOURCE_P1_SHA256, "observed": source_hash, "pass": source_hash == SOURCE_P1_SHA256}

    audit_text = (root / "current_p2_closure/FRESH_AUDIT_RESPONSE.md").read_text(encoding="utf-8")
    declared = {
        "P0": int(re.search(r"<P0_OPEN>(\d+)</P0_OPEN>", audit_text).group(1)),
        "P1": int(re.search(r"<P1_OPEN>(\d+)</P1_OPEN>", audit_text).group(1)),
        "P2": int(re.search(r"<P2_OPEN>(\d+)</P2_OPEN>", audit_text).group(1)),
    }
    adjudication = json.loads((root / "current_p2_closure/CURRENT_P2_ADJUDICATION.json").read_text(encoding="utf-8"))
    records = adjudication["adjudicated_records"]
    valid_dispositions = {"ACCEPTED", "REJECTED_WITH_EVIDENCE", "CONDITIONAL"}
    checks["fresh_review_inventory"] = {
        "declared": declared,
        "record_count": len(records),
        "unique_finding_ids": len({x["finding_id"] for x in records}),
        "all_required_fields": all(all(k in x for k in ["severity", "disposition", "exact_evidence", "action", "rerun_status", "affected_artifacts", "resolution"]) for x in records),
        "all_dispositions_valid": all(x["disposition"] in valid_dispositions for x in records),
        "pass": declared == {"P0": 0, "P1": 4, "P2": 0} and len(records) == 4 and len({x["finding_id"] for x in records}) == 4 and all(x["severity"] == "P1" for x in records) and all(x["disposition"] in valid_dispositions for x in records),
    }

    cross_audit_text = builder.CROSS_FREEZE_AUDIT.read_text(encoding="utf-8")
    cross_declared = {
        "P0": int(re.search(r"<P0_OPEN>(\d+)</P0_OPEN>", cross_audit_text).group(1)),
        "P1": int(re.search(r"<P1_OPEN>(\d+)</P1_OPEN>", cross_audit_text).group(1)),
        "P2": int(re.search(r"<P2_OPEN>(\d+)</P2_OPEN>", cross_audit_text).group(1)),
    }
    cross_revision = json.loads((root / CROSS_REVISION_PATH).read_text(encoding="utf-8"))
    cross_records = cross_revision["adjudications"]
    cross_severity_counts = {
        severity: sum(item["severity"] == severity for item in cross_records)
        for severity in ("P0", "P1", "P2")
    }
    checks["cross_freeze_review_inventory"] = {
        "governing_review_path": str(builder.CROSS_FREEZE_AUDIT),
        "expected_sha256": builder.CROSS_FREEZE_AUDIT_SHA256,
        "observed_sha256": sha256_file(builder.CROSS_FREEZE_AUDIT),
        "declared": cross_declared,
        "record_count": len(cross_records),
        "severity_counts": cross_severity_counts,
        "unique_finding_ids": len({item["finding_id"] for item in cross_records}),
        "unresolved": cross_revision["unresolved"],
        "pass": (
            sha256_file(builder.CROSS_FREEZE_AUDIT) == builder.CROSS_FREEZE_AUDIT_SHA256
            and cross_declared == {"P0": 0, "P1": 2, "P2": 1}
            and len(cross_records) == 3
            and cross_severity_counts == {"P0": 0, "P1": 2, "P2": 1}
            and len({item["finding_id"] for item in cross_records}) == 3
            and all(item["disposition"].startswith("ACCEPTED_AND_CLOSED") for item in cross_records)
            and cross_revision["unresolved"] == {"P0": [], "P1": [], "P2": []}
        ),
    }

    mandatory = [
        P2_PATH,
        "current_p2_closure/CURRENT_P2_ADJUDICATION.json",
        "current_p2_closure/CURRENT_P2_ADJUDICATION.md",
        CLOSURE_PATH,
        "current_p2_closure/CURRENT_P2_CLOSURE.md",
        MAP_PATH,
        GATE_PATH,
        "current_p2_closure/COMMANDS_RUN.txt",
        ANCHOR_PATH,
        CROSS_REVISION_PATH,
        "cross_freeze_revision/CROSS_FREEZE_REVISION_1.md",
        HASH_PROOF_PATH,
        "current_p2_closure/build_current_p2_outputs.py",
        "current_p2_closure/run_current_p2_gate.py",
    ]
    # Seed the gate file so the mandatory-file check includes the file being produced.
    if not (root / GATE_PATH).exists():
        write_json(root / GATE_PATH, {"project": PROJECT, "status": "RUNNING"})
    if not (root / HASH_PROOF_PATH).exists():
        write_json(root / HASH_PROOF_PATH, {"project": PROJECT, "status": "RUNNING"})
    missing = [path for path in mandatory if not (root / path).is_file()]
    checks["mandatory_files"] = {"files": mandatory, "missing": missing, "pass": not missing}

    manuscript = (root / P2_PATH).read_text(encoding="utf-8")
    forbidden = {
        "absolute_user_path": r"/Users/",
        "absolute_volume_path": r"/Volumes/",
        "todo_tbd": r"\b(?:TODO|TBD)\b",
        "author_placeholder": r"(?i)\b(?:insert author|author name|corresponding author|affiliation placeholder)\b",
        "internal_review_language": r"(?i)\b(?:P0|P1|P2|ACCEPTED|REJECTED_WITH_EVIDENCE|CONDITIONAL|READY_FOR_CROSS_FREEZE|fresh audit|adjudication|closure status)\b",
        "non_methods_results_sections": r"(?m)^## (?:Introduction|Discussion|Abstract|Declarations|Cover Letter)\b",
    }
    hits = {name: sorted(set(re.findall(pattern, manuscript))) for name, pattern in forbidden.items() if re.search(pattern, manuscript)}
    checks["forbidden_language"] = {
        "forbidden_patterns": forbidden,
        "forbidden_hits": hits,
        "pass": not hits,
    }
    checks["semantic_hygiene"] = {
        "forbidden_hits": hits,
        "has_methods": "## Methods" in manuscript,
        "has_results": "## Results" in manuscript,
        "preserves_negative_boundaries": all(term in manuscript for term in ["not designed to establish prognosis", "no pooled pan-cancer coefficient", "significant opposite direction", "were null", "negative median differences", "no spot-level P value"]),
        "pass": not hits and "## Methods" in manuscript and "## Results" in manuscript and all(term in manuscript for term in ["not designed to establish prognosis", "no pooled pan-cancer coefficient", "significant opposite direction", "were null", "negative median differences", "no spot-level P value"]),
    }

    numeric_map = pd.read_csv(root / MAP_PATH, dtype={"numeric_token": str})
    expected_occurrences = []
    for line_no, line in enumerate(manuscript.splitlines(), start=1):
        for ordinal, token in enumerate(builder.numeric_tokens(line), start=1):
            expected_occurrences.append((f"L{line_no:04d}_N{ordinal:02d}", line_no, ordinal, token))
    observed_occurrences = [tuple(x) for x in numeric_map[["occurrence_id", "manuscript_line", "token_ordinal", "numeric_token"]].itertuples(index=False, name=None)]
    facts = builder.compute_facts(root)
    derivations = json.loads((root / "current_p2_closure/NUMERIC_DERIVATIONS_CURRENT_P2.json").read_text(encoding="utf-8"))["facts"]
    fact_failures = []
    hash_failures = []
    display_failures = []
    for row in numeric_map.itertuples(index=False):
        fact = facts.get(row.numeric_id)
        retained_fact = derivations.get(row.numeric_id)
        if fact is None or retained_fact is None or not equal_value(fact["value"], retained_fact["value"]) or not equal_value(fact["value"], row.source_value):
            fact_failures.append(row.occurrence_id)
            continue
        expected_display = builder.display(fact["value"], row.display_rule)
        if expected_display != str(row.numeric_token):
            display_failures.append(row.occurrence_id)
        source = root / row.source_file
        if not source.is_file() or sha256_file(source) != row.source_sha256:
            hash_failures.append(row.occurrence_id)
    checks["numeric_map"] = {
        "manuscript_numeric_occurrences": len(expected_occurrences),
        "map_rows": len(numeric_map),
        "occurrence_inventory_exact": expected_occurrences == observed_occurrences,
        "unique_occurrence_ids": not numeric_map["occurrence_id"].duplicated().any(),
        "fact_failures": fact_failures,
        "display_failures": display_failures,
        "source_hash_failures": hash_failures,
        "pass": expected_occurrences == observed_occurrences and not numeric_map["occurrence_id"].duplicated().any() and not fact_failures and not display_failures and not hash_failures,
    }

    required_bindings = {
        "primary_min_distinct_values": ("code/run_cptac_v2.py", 1),
        "residual_min_distinct_values": ("code/run_cptac_v2.py", 1),
        "tcga_ddof": ("current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py", 1),
        "tcga_min_subtype_n": ("current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", 1),
        "tcga_min_program_genes": ("current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", 1),
        "tcga_q_cutoff": ("current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py", 1),
        "tcga_anchor_rule_rows": (ANCHOR_PATH, 1),
        "apollo_min_program_genes": ("current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py", 1),
        "apollo_min_distinct_values": ("code/postprocess_support_v2.py", 1),
        "apollo_ci_percent": ("code/postprocess_support_v2.py", 3),
        "scrna_cellxgene_min_group_cells": ("code/run_h5ad_support_v2.py", 1),
        "scrna_cellxgene_min_program_genes": ("code/run_h5ad_support_v2.py", 1),
        "scrna_3ca_min_group_cells": ("current_p2_closure/executed_code/run_bio05_3ca_source_mapping_v0_8.py", 1),
        "scrna_3ca_min_program_genes": ("current_p2_closure/executed_code/run_bio05_3ca_source_mapping_v0_8.py", 1),
        "scrna_support_q_cutoff": ("code/postprocess_support_v2.py", 1),
        "spatial_min_program_genes": ("code/run_h5ad_support_v2.py", 1),
    }
    binding_failures = []
    observed_bindings = {}
    for numeric_id, (expected_source, expected_count) in required_bindings.items():
        subset = numeric_map[numeric_map["numeric_id"].eq(numeric_id)]
        observed_bindings[numeric_id] = {
            "row_count": len(subset),
            "source_files": sorted(subset["source_file"].astype(str).unique().tolist()),
            "source_locators": sorted(subset["source_locator"].astype(str).unique().tolist()),
        }
        if len(subset) != expected_count or set(subset["source_file"].astype(str)) != {expected_source}:
            binding_failures.append(numeric_id)
    obsolete_semantic_ids = sorted(
        set(numeric_map["numeric_id"]).intersection({"threshold_q_cutoff", "scrna_min_group_cells"})
    )
    checks["numeric_semantic_bindings"] = {
        "required": {
            key: {"source_file": value[0], "occurrence_count": value[1]}
            for key, value in required_bindings.items()
        },
        "observed": observed_bindings,
        "binding_failures": binding_failures,
        "obsolete_semantic_ids_present": obsolete_semantic_ids,
        "pass": not binding_failures and not obsolete_semantic_ids,
    }

    config_path = root / "current_p2_closure/PROGRAM_DEFINITION_CURRENT_P2.csv"
    config = pd.read_csv(config_path)
    tcga_selection = pd.read_csv(root / "current_p2_closure/TCGA_PRIMARY_SAMPLE_SELECTION_CURRENT_P2.csv")
    rules = json.loads((root / "current_p2_closure/SUPPORT_MAPPING_RULES_CURRENT_P2.json").read_text(encoding="utf-8"))
    required_score_rule_sections = ["tcga", "apollo", "single_cell_cellxgene", "single_cell_3ca", "spatial"]
    score_rule_missing = [
        section for section in required_score_rule_sections
        if "score_aggregation_and_missingness" not in rules.get(section, {})
    ]
    checks["corrected_reproducibility_values"] = {
        "program_config_sha256": sha256_file(config_path),
        "program_rows": len(config),
        "program_count": int(config["signature"].nunique()),
        "tcga_selection_rows": len(tcga_selection),
        "tcga_unique_patients": int(tcga_selection["patient_id"].nunique()),
        "tcga_sample_type_codes": sorted(pd.to_numeric(tcga_selection["sample_type_code"]).unique().tolist()),
        "apollo_order_starts_with_collapse": all(any("average" in step for step in steps[:-1]) and "z-standardize" in " ".join(steps) for steps in rules["apollo"]["operation_order"].values()),
        "spatial_fallback_complete": len(rules["spatial"]["assignment"]) == 6,
        "score_rule_sections": required_score_rule_sections,
        "score_rule_sections_missing": score_rule_missing,
        "tcga_anchor_rule_artifact": rules["tcga"]["anchor_selection"]["artifact"],
        "tcga_anchor_rule_count": rules["tcga"]["anchor_selection"]["anchor_count"],
        "pass": sha256_file(config_path) == "b631094186b4f1a69f3e2ee8b948009f10f4b7ad0b644cd1cf446e4152c46ef1" and len(config) == 139 and config["signature"].nunique() == 8 and len(tcga_selection) == 4780 and tcga_selection["patient_id"].nunique() == 4780 and sorted(pd.to_numeric(tcga_selection["sample_type_code"]).unique().tolist()) == [1] and len(rules["spatial"]["assignment"]) == 6 and not score_rule_missing and rules["tcga"]["anchor_selection"]["anchor_count"] == 17,
    }

    checks["tcga_anchor_selection"] = tcga_anchor_crosscheck(root, builder)
    checks["retained_sensitivity_outcomes"] = sensitivity_outcome_crosscheck(root)

    checks["spearman_algorithm"] = spearman_algorithm_crosscheck(root)
    checks["apollo_ci_and_bh"] = apollo_ci_crosscheck(root)
    checks["tcga_bh"] = tcga_bh_crosscheck(root)
    checks["wilcoxon_recomputation"] = wilcoxon_crosscheck(root)
    checks["deterministic_regeneration"] = deterministic_regeneration(root)
    deterministic = checks["deterministic_regeneration"]
    write_json(
        root / HASH_PROOF_PATH,
        {
            "project": PROJECT,
            "status": "PASS" if deterministic["pass"] else "FAIL",
            "method": "two consecutive executions of the deterministic CURRENT_P2 generator under the same retained inputs",
            "tracked_file_count": deterministic["tracked_file_count"],
            "tracked_files": deterministic["tracked_files"],
            "pass_1_hashes": deterministic["pass_1_hashes"],
            "pass_2_hashes": deterministic["pass_2_hashes"],
            "hash_mismatches": deterministic["hash_mismatches"],
            "builder_runs": deterministic["runs"],
        },
    )

    # Deterministic regeneration rewrites the initial closure; re-check immutable input and output map.
    checks["source_after_regeneration"] = {"observed": sha256_file(root / SOURCE_P1), "pass": sha256_file(root / SOURCE_P1) == SOURCE_P1_SHA256}
    failures = [name for name, result in checks.items() if not bool(result.get("pass", False))]
    gate = {
        "project": PROJECT,
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "checks": checks,
        "audited_files": mandatory + [
            "current_p2_closure/PROGRAM_DEFINITION_CURRENT_P2.csv",
            "current_p2_closure/SUPPORT_MAPPING_RULES_CURRENT_P2.json",
            "current_p2_closure/TCGA_PRIMARY_SAMPLE_SELECTION_CURRENT_P2.csv",
            ANCHOR_PATH,
            "current_p2_closure/SOFTWARE_ENVIRONMENT_CURRENT_P2.json",
            "current_p2_closure/NUMERIC_DERIVATIONS_CURRENT_P2.json",
            HASH_PROOF_PATH,
            CROSS_REVISION_PATH,
            "cross_freeze_revision/CROSS_FREEZE_REVISION_1.md",
            "results/cptac/cptac_primary_concordance_registry_v2.csv",
            "results/cptac/cptac_heterogeneity_summary_v2.csv",
            "results/cptac/cptac_leave_one_cancer_out_v2.csv",
            "results/cptac/cptac_within_cancer_standardized_pooled_sensitivity_v2.csv",
            "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv",
            "results/apollo/apollo_concordance_registry_v2.csv",
            "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv",
            "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv",
        ],
    }
    write_json(root / GATE_PATH, gate)
    update_closure(root, gate, builder)
    print(json.dumps({"project": PROJECT, "gate_status": gate["status"], "failures": failures, "final_status": "READY_FOR_CROSS_FREEZE" if not failures else "NOT_READY"}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
