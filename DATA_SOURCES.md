# Data sources

The analysis uses processed public data from CPTAC, TCGA/UCSC Xena, APOLLO, CZ CELLxGENE Discover, the Curated Cancer Cell Atlas, and the Wu et al. breast-cancer single-cell/spatial atlas.

The consolidated workbook is the numerical authority for the publication
layer. Its `Contents` worksheet maps every main, supplementary, and supporting
table to a machine-readable CSV path. The `analysis/` directory contains the
scientific code that recomputes the principal registries from prepared public
inputs and checks them against those locked tables.

None of the 105 upstream assets recorded in
`provenance/SOURCE_MANIFEST_PUBLIC.csv` is included in the repository. Each is
therefore marked `NOT_REDISTRIBUTED`; this describes the package boundary and
does not replace the terms imposed by the original source. The audited
source-family decisions and official terms routes are in
`provenance/SOURCE_TERMS_AUDIT.md`.

This repository is a scientific-analysis plus figure/table reproducibility
package. It does not redistribute upstream data and does not automate every
download or source-specific preparation step.
