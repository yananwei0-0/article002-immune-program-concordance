# Article 002 — publication-figure reproducibility release

Release candidate: `v1.0.0-rc1` (2026-09-29)

This repository candidate reproduces the four main figures, six supplementary figures, and 22 machine-readable analytical tables for:

> RNA–protein and RNA–phosphoprotein concordance of tumor immune programs across ten CPTAC cancers

The technical package is complete and validated. It is **not yet approved for public release** because author metadata, the final license, repository URL, and archived release DOI still require author input. See `PUBLIC_RELEASE_BLOCKERS.md`.

## Reproducibility boundary

This is the complete publication-figure and table-export layer. It reads already-computed, locked result registries and rebuilds the publication figures and CSV tables. It does not replace the upstream pipelines that derived scores, correlations, bootstrap intervals, multiple-testing corrections, or compartment contrasts from raw molecular data.

No raw participant-, donor-, cell-, or facility-level data are redistributed here.

## Contents

- `plot_all_figures.py` — validates the locked registries and generates all 10 figures.
- `export_tables.py` — exports all 22 analytical worksheets as deterministic CSV files.
- `data/Article2_Consolidated_Tables.xlsx` — consolidated numerical source with 23 worksheets: one contents sheet and 22 data sheets.
- `tables/` — three main-table CSVs, 16 supplementary-table CSVs, three supporting CSVs, and an export manifest.
- `generated/` — PNG (300 dpi), editable SVG, TIFF (600 dpi), and a portable run manifest.
- `previews/` — contact sheets for rapid visual review.
- `manuscript/FULL_MANUSCRIPT_RELEASE_CANDIDATE.md` — corrected text with populated main tables and supplementary-table index.
- `provenance/` — public source manifests and the single-cell/spatial citation maps.
- `validation/` — machine-readable validation report and human-readable release audit.
- `FIGURE_DATA_MAP.csv` — figure-to-worksheet mapping.
- `requirements.txt` — portable dependency ranges.
- `requirements-lock.txt` — exact tested plotting environment.

## Rebuild

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python plot_all_figures.py
python export_tables.py
python build_previews.py
python validate_release.py
```

For a broader compatible environment, install `requirements.txt` instead of the lock file. Generate a subset of figures with:

```bash
python plot_all_figures.py \
  --only Fig2_heterogeneity_and_scoring_robustness FigS3_TCGA_full_correlation_atlas
```

Use a different workbook or output directory with:

```bash
python plot_all_figures.py --workbook /path/to/results.xlsx --outdir /path/to/figures
```

## Validated invariants

The release validation checks the following locked quantities:

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

## Data and provenance

Public source identities and file-level provenance are summarized in `DATA_ACCESS.md`, `DATA_SOURCES.md`, and `provenance/`. The source manifest retains `NOT_ASSESSED` redistribution flags where source-specific terms still require author review; the raw assets themselves are not included.

## Release policy

Do not make this candidate public or create a final GitHub/Zenodo release until every item in `PUBLIC_RELEASE_BLOCKERS.md` is resolved and `RELEASE_STATUS.json` is changed to `public_release_ready: true` after a final validation run.
