# Scientific-analysis reproduction layer

This directory contains the code used to recompute the study's principal
analytical registries from prepared public-source matrices. It is separate from
the repository-root publication layer, which converts the locked registries
into journal figures and tables.

## What this layer recomputes

| Layer | Recomputed quantities | Principal locked output |
| --- | --- | --- |
| CPTAC | fixed-program scores; matched RNA–protein and RNA–phosphoprotein Spearman correlations; 20,000-resample paired-case BCa intervals; exact Monte Carlo permutation tests where required; BH correction; leave-one-case, feature-row, coverage-threshold, pooled, residual, and context analyses | `tables/supplementary/S3_Primary.csv` through `S5_ResidualContext.csv` |
| TCGA | primary-tumor sample selection; eight RNA program scores; 5,440 within-cancer correlations against 68 immune metrics; BH correction; immune-subtype summaries | `tables/supplementary/S6a_TCGA_Full.csv` through `S6d_TCGA_Subtypes.csv` |
| APOLLO | LUAD/OV program scores; matched RNA–protein and RNA–phosphoprotein correlations; bootstrap intervals; discovery-consistency classification | `tables/supplementary/S7_APOLLO.csv` |
| Single-cell and 3CA | memory-bounded H5AD reading; donor-compartment pseudobulk scores; paired donor contrasts; source-bounded Wilcoxon tests and BH correction | `tables/supplementary/S8a_scRNA_Global.csv`, `S8b_scRNA_Source.csv` |
| Spatial BRCA | spot-level fixed-program scores; sample-compartment summaries; paired four-sample contrasts | `tables/supplementary/S9_Spatial_BRCA.csv` |

The repository does not redistribute the upstream molecular files and does not
automate every source download or source-specific harmonization step. Users
must obtain the public files under their original terms and place the prepared
inputs in the layout described in `INPUT_LAYOUT.md`. This boundary does not
affect the included statistical analysis code: the score construction,
correlations, intervals, multiple-testing corrections, sensitivity analyses,
and support-layer summaries are all present.

## Environment

From the repository root:

```bash
python3 -m venv .venv-analysis
source .venv-analysis/bin/activate
python -m pip install -r analysis/environment/requirements.lock
python -m pip install -e analysis
```

Python 3.11–3.14 is supported by the packaged project metadata. The exact
historical environment is recorded under `environment/`; a compatible public
lock is provided for reruns.

## Preflight and run

Inspect the staged commands without running the large analyses:

```bash
python analysis/run_pipeline.py \
  --stage all \
  --data-root /path/to/prepared_cptac_tcga_project \
  --external-root /path/to/public_external_inputs \
  --analysis-root /path/to/recomputed_article002 \
  --dry-run
```

Run all scientific layers and the final parity gate by removing `--dry-run`.
Individual stages can be selected, for example:

```bash
python analysis/run_pipeline.py \
  --stage cptac tcga \
  --data-root /path/to/prepared_cptac_tcga_project \
  --analysis-root /path/to/recomputed_article002
```

The H5AD stage is intentionally memory bounded. `--source-cache-dir` may point
to fast scratch storage; source bytes are size checked before analysis.

## Numerical parity gate

After recomputation, compare the scientific outputs directly with the locked
release tables:

```bash
python analysis/compare_recomputed_to_locked.py \
  --analysis-root /path/to/recomputed_article002 \
  --repo-root .
```

The gate checks 14 registries by stable keys and compares every published
column with numerical tolerances of `rtol=1e-10` and `atol=1e-12`. The retained
historical rerun evidence passed all 14 comparisons; see
`provenance/HISTORICAL_LOCKED_TABLE_PARITY.json`.

## Provenance boundary

- `definitions/` contains the fixed gene programs, public path-sanitized locked
  protocol, anchor selections, and mapping rules.
- `input_manifests/` contains public URLs, expected sizes, source identifiers,
  and repository-relative landing paths for single-cell and spatial inputs.
- `provenance/RUN_MANIFEST_PUBLIC.json` records the historical real-source run
  and hashes of its principal outputs.
- `provenance/HISTORICAL_LOCKED_TABLE_PARITY.json` independently maps those
  retained outputs to the publication tables.
- `provenance/PORTABLE_REAL_SOURCE_SMOKE_20261002.json` records the fresh
  public-path TCGA and APOLLO smoke reruns performed for this release.
- `code/current_p2_closure/executed_code/` contains the source modules used by
  the final analysis, with only path/CLI portability changes in this release.

The original locked protocol hash and the path-sanitized public hash are both
recorded in `provenance/PATH_PORTABILITY_TRANSFORMATION.json`.
