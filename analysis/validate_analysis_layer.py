#!/usr/bin/env python3
"""Static and provenance validation for the public scientific-analysis layer."""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
EXPECTED_PROTOCOL_SHA256 = "f419f5ac589bf6d6ec61abb6663a47a5114b86a7d82980327b298adccb46aefc"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    required = [
        ANALYSIS / "README.md",
        ANALYSIS / "INPUT_LAYOUT.md",
        ANALYSIS / "METHOD_TO_CODE_MAP.csv",
        ANALYSIS / "run_pipeline.py",
        ANALYSIS / "compare_recomputed_to_locked.py",
        ANALYSIS / "code/run_cptac_v2.py",
        ANALYSIS / "code/run_h5ad_support_v2.py",
        ANALYSIS / "code/postprocess_support_v2.py",
        ANALYSIS / "definitions/BIO-05_V2_LOCKED_PROTOCOL.json",
        ANALYSIS / "definitions/signature_gene_sets_v0_1.csv",
        ANALYSIS / "provenance/HISTORICAL_LOCKED_TABLE_PARITY.json",
        ANALYSIS / "provenance/PORTABLE_REAL_SOURCE_SMOKE_20261002.json",
        ANALYSIS / "provenance/RUN_MANIFEST_PUBLIC.json",
        ANALYSIS / "input_manifests/public_validation_queue.tsv",
        ANALYSIS / "input_manifests/brca_spatial_queue.tsv",
    ]
    require(all(path.is_file() and path.stat().st_size > 0 for path in required), "Required analysis-layer file missing")

    protocol = ANALYSIS / "definitions/BIO-05_V2_LOCKED_PROTOCOL.json"
    require(sha256(protocol) == EXPECTED_PROTOCOL_SHA256, "Public locked-protocol hash mismatch")

    parity = json.loads((ANALYSIS / "provenance/HISTORICAL_LOCKED_TABLE_PARITY.json").read_text(encoding="utf-8"))
    require(parity.get("status") == "PASS", "Historical locked-table parity is not PASS")
    require(parity.get("comparisons_passed") == parity.get("comparisons_total") == 14, "Expected 14/14 parity comparisons")
    smoke = json.loads((ANALYSIS / "provenance/PORTABLE_REAL_SOURCE_SMOKE_20261002.json").read_text(encoding="utf-8"))
    require(smoke.get("status") == "PASS", "Portable real-source smoke is not PASS")
    require(smoke.get("cptac", {}).get("locked_table_comparisons") == "6/6 PASS", "CPTAC portable smoke parity failed")
    require(smoke.get("tcga", {}).get("exact_historical_hash_matches") == 7, "TCGA portable smoke hash count mismatch")
    require(smoke.get("apollo", {}).get("core_concordance_columns_exact") is True, "APOLLO portable smoke core comparison failed")
    require(len(smoke.get("apollo", {}).get("exact_historical_hash_matches", {})) == 4, "APOLLO portable smoke hash count mismatch")

    python_files = sorted((ANALYSIS / "code").rglob("*.py")) + [
        ANALYSIS / "run_pipeline.py",
        ANALYSIS / "compare_recomputed_to_locked.py",
    ]
    for path in python_files:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    queue_rows = 0
    for queue_path in [
        ANALYSIS / "input_manifests/public_validation_queue.tsv",
        ANALYSIS / "input_manifests/brca_spatial_queue.tsv",
    ]:
        with queue_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        require(rows, f"Empty input queue: {queue_path.name}")
        for row in rows:
            remote_path = Path(row["remote_path"])
            require(not remote_path.is_absolute(), f"Queue path must be relative: {remote_path}")
            require(str(row["url"]).startswith("https://"), f"Queue URL must use HTTPS: {row['url']}")
            expected_size = str(row["expected_size"]).strip()
            require(not expected_size or expected_size.isdigit(), f"Queue size is not numeric: {row['expected_size']}")
        queue_rows += len(rows)

    local_path_pattern = re.compile(r"/(?:Users|Volumes)/[A-Za-z0-9._ -]+/[A-Za-z0-9._ -]+")
    secret_pattern = re.compile(
        r"AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|sk-[A-Za-z0-9]{20,}|BEGIN (?:RSA |OPENSSH )?PRIVATE KEY"
    )
    local_hits: list[str] = []
    secret_hits: list[str] = []
    suffixes = {".py", ".md", ".txt", ".csv", ".tsv", ".json", ".toml", ".lock", ".in", ".example"}
    for path in ANALYSIS.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in suffixes or "__pycache__" in path.parts:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if local_path_pattern.search(text):
            local_hits.append(path.relative_to(ROOT).as_posix())
        if secret_pattern.search(text):
            secret_hits.append(path.relative_to(ROOT).as_posix())
    require(not local_hits, f"Concrete machine-local paths found: {local_hits}")
    require(not secret_hits, f"Possible credentials found: {secret_hits}")

    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "required_files": len(required),
        "python_files_parsed": len(python_files),
        "input_queue_rows": queue_rows,
        "locked_table_parity": "14/14 PASS",
        "portable_real_source_smoke": "PASS",
        "protocol_sha256": EXPECTED_PROTOCOL_SHA256,
        "concrete_local_path_hits": 0,
        "credential_like_hits": 0,
        "boundary": "Static/public-package validation plus retained historical real-source parity; not a fresh download-and-clean rerun.",
    }
    output = ANALYSIS / "validation/ANALYSIS_LAYER_VALIDATION.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
