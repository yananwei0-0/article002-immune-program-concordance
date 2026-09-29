# Data access

Raw participant-, donor-, cell-, spot-, and facility-level files are not redistributed in this release candidate.

The locked workbook and exported CSV files contain the aggregated analytical registries required to reproduce the publication figures. Public data should be retrieved from the official resources named in the manuscript and the following provenance files:

- `provenance/SOURCE_MANIFEST_PUBLIC.csv` — file names, retained public URLs, byte sizes, checksums, and redistribution-review status.
- `provenance/SINGLE_CELL_SOURCE_CITATION_MAP.csv` — exact 10-source single-cell route and publication map.
- `provenance/BRCA_VISIUM_SOURCE_CITATION_MAP.csv` — exact four-sample spatial route, identifiers, and analyzed spot counts.

The four analyzed Visium samples are `CID4290`, `CID44971`, `CID4465`, and `CID4535`, all from CELLxGENE collection `dea97145-f712-431c-a223-6b5f565f362a` and the Wu et al. breast-cancer atlas (DOI `10.1038/s41588-021-00911-1`).

Rows marked `NOT_ASSESSED` in the source manifest require source-specific terms review before any raw asset is redistributed. This candidate does not redistribute those assets.
