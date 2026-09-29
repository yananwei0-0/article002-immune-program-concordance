# Single-cell and BRCA Visium citation audit

Audit date: 2026-09-01

## Adjudication

The combined single-cell layer comprises 10 analysis sources: 8 CELLxGENE collections/datasets and 2 study-specific 3CA routes. Each analysis source now has a corresponding source publication, exact analysis identifier, data-version locator, cancer coverage, and source-specific donor count in `SINGLE_CELL_SOURCE_CITATION_MAP.csv`.

The four breast Visium samples are `CID4290`, `CID44971`, `CID4465`, and `CID4535`. They are four distinct CELLxGENE datasets but all belong to collection `dea97145-f712-431c-a223-6b5f565f362a` and all cite Wu et al., *A single-cell and spatially resolved atlas of human breast cancers*, Nature Genetics (2021), DOI `10.1038/s41588-021-00911-1`. Their exact dataset UUIDs and H5AD version URLs are in `BRCA_VISIUM_SOURCE_CITATION_MAP.csv`.

## Important boundaries

- The breast single-cell source and all four Visium datasets come from the same Wu et al. collection. The spatial layer is therefore same-collection triangulation, not independent validation.
- The GBmap source is linked by both the frozen H5AD citation and the current CELLxGENE collection to DOI `10.1101/2022.08.27.505439`. This is a bioRxiv preprint; it is not represented as peer reviewed.
- The frozen HNSCC H5AD did not include a publication DOI. Its exact collection UUID and title now match the final 2026 Communications Medicine article by Kroehling et al., DOI `10.1038/s43856-026-01401-3`.
- The Regner source was downloaded from the 3CA ovarian page, but the executed analysis selected six `Endometrial Cancer` samples. Its corresponding source article is the combined ovarian/endometrial study by Regner et al., DOI `10.1016/j.molcel.2021.10.013`.
- Several analyzed sources are integrated atlases. The mapping identifies the publication corresponding to the analyzed atlas object; it does not imply that every constituent primary cohort in an integrated atlas was independently analyzed as a separate source.

## Verification basis

Identity was checked against the executed queue manifests, retained H5AD `uns/title` and `uns/citation` fields, the current CELLxGENE collection API, the official 3CA study labels, DOI registry metadata, and PubMed or publisher records where available. No manuscript result, denominator, or statistical claim was changed by this citation audit.
