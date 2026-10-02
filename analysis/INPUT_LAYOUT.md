# Public input layout

The analysis runner accepts two roots. No upstream source file is committed to
this repository.

## Prepared CPTAC/TCGA root (`--data-root`)

```text
DATA_ROOT/
├── 03_downloads/
│   ├── cptac_python_processed/                 # files referenced by the CPTAC inventory
│   └── tcga_xena_pancanatlas/
│       ├── EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
│       ├── TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz
│       └── Subtype_Immune_Model_Based.txt.gz
└── 04_processed/
    ├── discovery_v0_4/
    │   ├── cptac_bcm_transcriptomics_download_manifest_v0_4.csv
    │   └── bcm_matched_correlation_summary_v0_4.csv
    ├── first_pass_qc/
    │   └── cptac_processed_inventory_v0_1.csv
    ├── harmonization_v0_1/
    │   ├── cptac_sample_map_v0_1.csv
    │   ├── tcga_sample_map_v0_1.csv
    │   └── depmap_gene_feature_map_v0_1.csv
    └── public_pan_cancer_download_manifest_v0_1.csv
```

The inventory and manifest path fields are interpreted relative to
`DATA_ROOT`. The public fixed-program definition is supplied from
`analysis/definitions/signature_gene_sets_v0_1.csv` by the runner.

## External public-source root (`--external-root`)

```text
EXTERNAL_ROOT/
├── proteogenomics/
│   ├── APOLLO_LUAD_RNA_quant_files_20181226.tar.gz
│   ├── APOLLO1_Level3_Gprot_v061418_n87.csv
│   ├── APOLLO1_Phospho_Level3_n87.csv
│   ├── APOLLO_OV_WholeTumor_RNASeq_NormalizedCounts.csv
│   └── APOLLO_OV_WholeTumor_GlobalProteomics_imputed.csv
├── 3ca/
│   ├── Data_Peng2019_Pancreas.tar.gz
│   └── Data_Regner2021_Ovarian.tar.gz
├── scRNA/                                      # H5AD names from public_validation_queue.tsv
└── spatial/                                    # H5AD names from brca_spatial_queue.tsv
```

The two queue TSV files provide official source URLs and expected byte sizes.
The broader source inventories and checksums are under `analysis/provenance/`
and the repository-root `provenance/` directory.

## Data-use boundary

Retrieve each source through its official route and comply with the governing
terms. The repository licenses apply only to repository-owned code and content;
they do not relicense CPTAC, TCGA, APOLLO, CELLxGENE, 3CA, or other upstream
data.
