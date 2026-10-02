# Data access

Raw participant-, donor-, cell-, spot-, and facility-level files are not redistributed in this release candidate.

The locked workbook and exported CSV files contain the aggregated analytical registries required to reproduce the publication figures. Public data should be retrieved from the official resources named in the manuscript and the following provenance files:

- `provenance/SOURCE_MANIFEST_PUBLIC.csv` — file names, retained public URLs, byte sizes, checksums, and package redistribution status.
- `provenance/SOURCE_TERMS_AUDIT.md` — source-family terms review and public-package decisions.
- `provenance/SINGLE_CELL_SOURCE_CITATION_MAP.csv` — exact 10-source single-cell route and publication map.
- `provenance/BRCA_VISIUM_SOURCE_CITATION_MAP.csv` — exact four-sample spatial route, identifiers, and analyzed spot counts.

The four analyzed Visium samples are `CID4290`, `CID44971`, `CID4465`, and `CID4535`, all from CELLxGENE collection `dea97145-f712-431c-a223-6b5f565f362a` and the Wu et al. breast-cancer atlas (DOI `10.1038/s41588-021-00911-1`).

All 105 source-manifest rows are marked `NOT_REDISTRIBUTED`, and none of those
upstream asset filenames is present in the repository tree. Users should obtain
source data from the official routes and comply with the terms of the original
resource; the licenses in this repository do not relicense upstream data.
