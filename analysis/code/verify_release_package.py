#!/usr/bin/env python3
"""Verify package checksums and the current Methods/Results authority hash."""
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

failures = []
rows = list(csv.DictReader((ROOT / "CHECKSUMS_SHA256.csv").open(encoding="utf-8")))
for row in rows:
    path = ROOT / row["relative_path"]
    if not path.exists(): failures.append({"check": "checksum", "file": row["relative_path"], "reason": "missing"})
    elif digest(path) != row["sha256"]: failures.append({"check": "checksum", "file": row["relative_path"], "reason": "mismatch"})
authority = json.loads((ROOT / "provenance/CURRENT_AUTHORITY_RECORD.json").read_text(encoding="utf-8"))
methods = list((ROOT / "manuscript").glob("*_METHODS_RESULTS_CURRENT_P2.md"))
if len(methods) != 1 or digest(methods[0]) != authority["current_methods_results_sha256"]:
    failures.append({"check": "methods_results_authority", "reason": "missing_or_hash_mismatch"})
gate = json.loads((ROOT / "provenance/NUMERIC_GATE_CURRENT_P2.json").read_text(encoding="utf-8"))
gate_status = gate.get("status", gate.get("overall"))
if gate_status != "PASS": failures.append({"check": "current_p2_numeric_gate", "reason": str(gate_status)})
print(json.dumps({"status": "PASS" if not failures else "FAIL", "checksum_entries": len(rows), "failures": failures}, indent=2))
raise SystemExit(1 if failures else 0)
