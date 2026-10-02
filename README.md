# Article 002 — scientific-analysis and publication reproducibility release

Submission release candidate: `v1.0.0-rc3` (2026-10-02)

Public repository: <https://github.com/yananwei0-0/article002-immune-program-concordance>

This repository reproduces the principal scientific result registries, four
main figures, six supplementary figures, and 22 machine-readable analytical
tables for:

> RNA–protein and RNA–phosphoprotein concordance of tumor immune programs across ten CPTAC cancers

## Reproducibility boundary

The release contains two connected layers:

1. `analysis/` recomputes the principal CPTAC, TCGA, APOLLO, single-cell/3CA,
   and spatial statistics from prepared public-source matrices. It includes
   program scoring, correlations, bootstrap intervals, conditional permutation
   tests, BH correction, robustness analyses, and support-layer assembly.
2. The repository-root scripts read the locked analytical workbook and rebuild
   all publication figures and CSV tables.

The upstream public molecular files are not redistributed, and source download
and source-specific preparation are not fully automated. Official locators,
expected file sizes, checksums where available, and the required input layout
are documented. This is therefore a scientific-analysis reproduction package
from prepared public inputs, not a redistribution of raw participant-, donor-,
cell-, spot-, or facility-level data.

## Contents

- `analysis/README.md` — scientific-analysis scope, execution, and parity gate.
- `analysis/run_pipeline.py` — unified stage runner for the principal analyses.
- `analysis/compare_recomputed_to_locked.py` — 14-registry numerical parity gate.
- `analysis/definitions/` — fixed programs, public path-sanitized protocol, anchors, and mapping rules.
- `analysis/input_manifests/` — public single-cell/spatial URLs, expected sizes, and relative landing paths.
- `plot_all_figures.py` — validates the locked registries and generates all 10 figures.
- `export_tables.py` — exports all 22 analytical worksheets as deterministic CSV files.
- `data/Article2_Consolidated_Tables.xlsx` — 23-sheet consolidated numerical source.
- `tables/` — three main-table CSVs, 16 supplementary-table CSVs, three supporting CSVs, and an export manifest.
- `generated/` — PNG (300 dpi), editable SVG, TIFF (600 dpi), and a portable run manifest.
- `manuscript/FULL_MANUSCRIPT_RELEASE_CANDIDATE.md` — corrected submission text with populated main tables and supplementary-table index.
- `provenance/` and `analysis/provenance/` — source manifests, execution records, and numerical-parity evidence.
- `validation/` and `analysis/validation/` — machine-readable release and analysis-layer audits.

## Recompute the main analyses

See `analysis/INPUT_LAYOUT.md`, then from the repository root:

```bash
python3 -m venv .venv-analysis
source .venv-analysis/bin/activate
python -m pip install -r analysis/environment/requirements.lock
python -m pip install -e analysis

python analysis/run_pipeline.py \
  --stage all \
  --data-root /path/to/prepared_cptac_tcga_project \
  --external-root /path/to/public_external_inputs \
  --analysis-root /path/to/recomputed_article002
```

Inspect the exact commands first by adding `--dry-run`. Individual stages may
be selected with `--stage cptac`, `--stage tcga`, and so on.

The final pipeline stage compares 14 recomputed registries directly with the
locked release tables. Retained historical real-source outputs pass all 14
comparisons, including all 5,440 TCGA correlations; the machine-readable report
is `analysis/provenance/HISTORICAL_LOCKED_TABLE_PARITY.json`.

## Rebuild figures and tables

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python plot_all_figures.py
python export_tables.py
python build_previews.py
python analysis/validate_analysis_layer.py
python validate_release.py
```

Generate a subset of figures with:

```bash
python plot_all_figures.py \
  --only Fig2_heterogeneity_and_scoring_robustness FigS3_TCGA_full_correlation_atlas
```

## Validated invariants

- 14/14 historical scientific registries match the locked publication tables column-for-column within declared numerical tolerances;
- 160 CPTAC registered rows, including 129 evaluable rows;
- 77 evaluable RNA–protein rows: all positive and 76 with `q < 0.05`;
- 52 evaluable RNA–phosphoprotein rows: all positive and all with `q < 0.05`;
- 5,440 TCGA correlations, including 4,850 with `q < 0.05`;
- 17 fixed TCGA program–metric anchors;
- APOLLO status counts of 11 supported, five null, one opposite, and seven not evaluable;
- 336 source-bounded single-cell rows, including 256 evaluable and 116 with within-family support;
- 24 spatial contrasts based on four paired breast-cancer samples;
- 10 figures in each of the PNG, SVG, and TIFF output directories;
- 22 CSV exports matching their workbook worksheets.

## Data, license, and citation

Public source identities and file-level provenance are summarized in
`DATA_ACCESS.md`, `DATA_SOURCES.md`, and the two provenance directories. None of
the upstream assets is redistributed.

Repository software is licensed under MIT. Original figure artwork, table
presentation, and documentation are licensed under CC BY 4.0, subject to the
exclusions in `LICENSE.md`. Upstream data are not relicensed. `CITATION.cff`
uses a collective author label so the code release can be cited without
publishing personal author metadata; named authors and an archival DOI may be
added later without changing the scientific outputs.
