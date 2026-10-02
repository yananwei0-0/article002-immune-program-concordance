# Source-terms and redistribution audit

Audit date: 2026-10-02

## Release boundary

The repository contains author-generated code, figures, documentation, a locked
results workbook, and derived table exports. It does **not** contain any of the
105 upstream assets recorded in `SOURCE_MANIFEST_PUBLIC.csv`. File-name
comparison against the repository tree found zero upstream-asset matches.

Accordingly, every source-manifest row is marked `NOT_REDISTRIBUTED`. This
status describes the package contents; it is not a claim that every upstream
source uses the same license. Users must obtain upstream data from the official
source and comply with that source's current terms.

## Source-family decisions

| Source family | Terms/access finding | Public-package decision |
| --- | --- | --- |
| CPTAC/PDC | PDC states that data use is governed by CC BY 4.0 and requests CPTAC and study attribution. | Derived summaries and code may be released with attribution; upstream files are not included. |
| TCGA/UCSC Xena | Xena provides compiled files from public resources; underlying GDC access classifications and source terms still apply. | Provide source routes and citations; do not redistribute the matrices. |
| APOLLO/GDC/PDC | Proteomic files follow PDC guidance; GDC data may be open or controlled depending on the file. | Release only derived cohort-level summaries and code; no source matrices or controlled-access files. |
| CZ CELLxGENE Discover | Published datasets are described as CC BY 4.0; current terms also prohibit attempted re-identification. | Release source-level derived summaries and exact citation maps; do not include H5AD files. |
| Curated Cancer Cell Atlas | The website distributes curated study data, but no blanket repository-wide reuse license was identified on the audited pages. | Do not redistribute archives; retain original-study citations and derived summaries only. |
| Wu breast-cancer Zenodo record | The public record is identifiable by DOI `10.5281/zenodo.4739739`; its displayed rights section does not name a license. | Do not redistribute matrices, images, or metadata archives; release only derived summaries and citation mapping. |
| DepMap Public 26Q1 | Used only as a feature-label dictionary. | Do not redistribute DepMap files; retain the release citation and derived harmonized labels. |

The machine-readable version of this decision table is
`SOURCE_TERMS_AUDIT.csv`. This audit documents the release boundary and is not
legal advice.

## Official routes checked

- PDC data-use guidelines: <https://proteomic.datacommons.cancer.gov/pdc/data-use-guidelines>
- GDC documentation: <https://docs.gdc.cancer.gov/>
- UCSC Xena: <https://xena.ucsc.edu/>
- APOLLO: <https://www.cancer.gov/about-nci/organization/cbiit/projects/apollo>
- CELLxGENE terms: <https://cellxgene.cziscience.com/tos>
- CELLxGENE data-contribution/licensing documentation: <https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data>
- 3CA: <https://www.weizmann.ac.il/sites/3CA/>
- Wu breast-cancer spatial record: <https://doi.org/10.5281/zenodo.4739739>
- DepMap portal: <https://depmap.org/portal/>

