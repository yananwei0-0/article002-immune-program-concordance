#!/usr/bin/env python3
"""Run the CPTAC pipeline in a retained audit sandbox and compare core outputs."""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


CORE_FILES = [
    "results/cptac/cptac_program_scores_v2.csv",
    "results/cptac/cptac_matched_gene_equal_scores_v2.csv",
    "results/cptac/cptac_primary_concordance_registry_v2.csv",
    "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv",
    "results/cptac/cptac_residual_context_registry_v2.csv",
    "tables/TABLE1_CPTAC_CASE_FLOW_V2.csv",
    "tables/TABLE2_CPTAC_PROGRAM_SUMMARY_V2.csv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--protocol", required=True)
    args = parser.parse_args()
    root = Path(args.analysis_root).resolve()
    rerun = root / "local_audit" / "rerun_execution"
    rerun.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(root / "code/run_cptac_v2.py"),
        "--analysis-root",
        str(rerun),
        "--data-root",
        str(Path(args.data_root).resolve()),
        "--protocol",
        str(Path(args.protocol).resolve()),
    ]
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", expand_legacy_paths("${SCI27_TMP_ROOT}/bio05_mpl_cache"))
    completed = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True)
    (root / "logs/rerun_consistency.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (root / "logs/rerun_consistency.stderr.log").write_text(completed.stderr, encoding="utf-8")
    comparisons = []
    for rel in CORE_FILES:
        original = root / rel
        repeated = rerun / rel
        comparisons.append(
            {
                "relative_path": rel,
                "original_exists": original.is_file(),
                "rerun_exists": repeated.is_file(),
                "original_sha256": sha256(original) if original.is_file() else None,
                "rerun_sha256": sha256(repeated) if repeated.is_file() else None,
                "byte_identical": original.is_file() and repeated.is_file() and sha256(original) == sha256(repeated),
            }
        )
    payload = {
        "command": command,
        "exit_code": completed.returncode,
        "core_file_count": len(CORE_FILES),
        "byte_identical_count": sum(x["byte_identical"] for x in comparisons),
        "all_core_outputs_byte_identical": completed.returncode == 0 and all(x["byte_identical"] for x in comparisons),
        "comparisons": comparisons,
    }
    out = root / "local_audit/RERUN_CONSISTENCY.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if payload["all_core_outputs_byte_identical"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
