# Analysis provenance

Current `v1.0.0-rc3` evidence:

- `HISTORICAL_LOCKED_TABLE_PARITY.json` — 14/14 retained real-source registries matched to locked publication tables.
- `PORTABLE_REAL_SOURCE_SMOKE_20261002.json` — fresh path-portable CPTAC, TCGA, and APOLLO runs.
- `PATH_PORTABILITY_TRANSFORMATION.json` — original and public locked-protocol hashes.
- `RUN_MANIFEST_PUBLIC.json` — path-redacted manifest of the original final scientific run.

The other JSON/CSV files in this directory are retained staging audit records.
Their historical fields such as `public_release_status` describe the state at
the time they were created and do not override the repository-root
`RELEASE_STATUS.json` or the current validation reports. The authoritative
redistribution and source-terms audit is in the repository-root `provenance/`
directory.
