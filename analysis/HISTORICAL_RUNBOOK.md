# Historical staging runbook

> This file is retained to document the earlier packaging state. It is not the
> current execution guide. Use `analysis/README.md`, `analysis/INPUT_LAYOUT.md`,
> and `analysis/run_pipeline.py` for the public `v1.0.0-rc3` workflow.

Run all commands from the repository root. Python dependency/path smoke has passed; a new full scientific rerun is not implied by that smoke test.

## 1. Python environment

The lock was resolved with Python 3.14. Use Python 3.14 for the closest reproduction target.

```bash
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install --index-url https://mirrors.aliyun.com/pypi/simple -r environment/requirements.lock
python -m pip install --no-deps -e .
python -m pip check
```

## 3. Paths and source data

Create local input directories, obtain raw files from `DATA_ACCESS.md`, and match them to `provenance/SOURCE_MANIFEST_PUBLIC.csv` by file identity/checksum.

```bash
cp config/paths.env.example .env
export SCI27_REPO_ROOT="$PWD"
export SCI27_RESEARCH_ROOT="/absolute/path/to/research_inputs"
export SCI27_DATA_ROOT="/absolute/path/to/public_data_lake"
export SCI27_BIOBANK_ROOT="/absolute/path/to/large_public_cache"
export SCI27_TMP_ROOT="${TMPDIR:-/tmp}/sci27-bio-05"
```

## 4. Package integrity

```bash
python code/verify_release_package.py
```

## 5. Analysis entrypoint inventory

The commands below are executable candidates recovered from the current analysis package. Check each script's arguments and source bindings before a full rerun.

| Role | Command template |
|---|---|
| Analysis candidate | `python code/run_cptac_v2.py` |
| Analysis candidate | `python code/run_h5ad_support_v2.py` |
| Analysis candidate | `python code/run_tcga_source_wrapper_v2.py` |
| Audit/closure | `python code/current_p2_closure/executed_code/run_bio05_3ca_source_mapping_v0_8.py` |
| Audit/closure | `python code/current_p2_closure/executed_code/run_bio05_apollo_external_cohort_v0_7.py` |
| Audit/closure | `python code/current_p2_closure/executed_code/run_bio05_scrna_source_mapping_v0_7.py` |
| Audit/closure | `python code/run_independent_audit_v2.py` |
| Audit/closure | `python code/check_rerun_consistency_v2.py` |
| Audit/closure | `python code/finalize_statistics_audit_cycle_v2.py` |
| Audit/closure | `python code/current_p1_closure/run_current_p1_numeric_gate.py` |
| Audit/closure | `python code/current_p2_closure/executed_code/run_pan_cancer_proteomic_immune_evasion_discovery_v0_1.py` |

The complete executable inventory is `provenance/ENTRYPOINT_INVENTORY.csv`.

## Evidence boundary

Historical numerical reproduction tier: `DETERMINISTIC_RERUN_PARITY_PASS`.
The original staging gate passed. Current release evidence is recorded in
`provenance/HISTORICAL_LOCKED_TABLE_PARITY.json`,
`provenance/PORTABLE_REAL_SOURCE_SMOKE_20261002.json`, and
`validation/ANALYSIS_LAYER_VALIDATION.json`.
