#!/usr/bin/env python3
"""Run the Article 002 scientific-analysis layers from public-source inputs.

The upstream molecular files are intentionally not redistributed. This runner
keeps all generated artifacts under one analysis root and invokes the exact
study modules for CPTAC, TCGA, APOLLO, single-cell/3CA, spatial, and final
post-processing. Use ``--dry-run`` to inspect commands without executing them.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from pathlib import Path


STAGE_ORDER = ["cptac", "tcga", "apollo", "threeca", "h5ad", "postprocess", "compare"]


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        nargs="+",
        choices=["preflight", "all", *STAGE_ORDER],
        default=["preflight"],
        help="One or more stages. 'all' runs every scientific layer and the parity gate.",
    )
    parser.add_argument("--repo-root", default=str(repo_root))
    parser.add_argument("--data-root", default="", help="Root of the CPTAC/TCGA prepared public inputs.")
    parser.add_argument(
        "--external-root",
        default="",
        help="Root containing proteogenomics/, 3ca/, and inputs/{scRNA,spatial}/.",
    )
    parser.add_argument("--analysis-root", default="analysis/recomputed")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--source-cache-dir", default="")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def required_files(repo: Path, data: Path | None, external: Path | None, stages: list[str]) -> list[Path]:
    items = [
        repo / "analysis/definitions/BIO-05_V2_LOCKED_PROTOCOL.json",
        repo / "analysis/definitions/signature_gene_sets_v0_1.csv",
    ]
    if data is not None and "cptac" in stages:
        items.extend(
            [
                data / "04_processed/first_pass_qc/cptac_processed_inventory_v0_1.csv",
                data / "04_processed/discovery_v0_4/cptac_bcm_transcriptomics_download_manifest_v0_4.csv",
                data / "04_processed/harmonization_v0_1/cptac_sample_map_v0_1.csv",
            ]
        )
    if data is not None and "tcga" in stages:
        items.extend(
            [
                data / "04_processed/harmonization_v0_1/tcga_sample_map_v0_1.csv",
                data / "04_processed/harmonization_v0_1/depmap_gene_feature_map_v0_1.csv",
                data / "04_processed/public_pan_cancer_download_manifest_v0_1.csv",
                data / "03_downloads/tcga_xena_pancanatlas/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
                data / "03_downloads/tcga_xena_pancanatlas/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz",
                data / "03_downloads/tcga_xena_pancanatlas/Subtype_Immune_Model_Based.txt.gz",
            ]
        )
    if data is not None and "apollo" in stages and "cptac" not in stages:
        items.append(data / "04_processed/discovery_v0_4/bcm_matched_correlation_summary_v0_4.csv")
    if external is not None and "apollo" in stages:
        apollo = external / "proteogenomics"
        items.extend(
            apollo / name
            for name in [
                "APOLLO_LUAD_RNA_quant_files_20181226.tar.gz",
                "APOLLO1_Level3_Gprot_v061418_n87.csv",
                "APOLLO1_Phospho_Level3_n87.csv",
                "APOLLO_OV_WholeTumor_RNASeq_NormalizedCounts.csv",
                "APOLLO_OV_WholeTumor_GlobalProteomics_imputed.csv",
            ]
        )
    if external is not None and "threeca" in stages:
        threeca = external / "3ca"
        items.extend(
            threeca / name
            for name in ["Data_Peng2019_Pancreas.tar.gz", "Data_Regner2021_Ovarian.tar.gz"]
        )
    if external is not None and "h5ad" in stages:
        for queue_name, accepted_type in [
            ("public_validation_queue.tsv", "cellxgene"),
            ("brca_spatial_queue.tsv", "cellxgene_spatial"),
        ]:
            queue_path = repo / "analysis/input_manifests" / queue_name
            with queue_path.open(encoding="utf-8", newline="") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                items.extend(
                    external / row["remote_path"]
                    for row in rows
                    if row.get("source_type") == accepted_type
                )
    return items


def run(command: list[str], env: dict[str, str], dry_run: bool) -> None:
    print("+ " + shlex.join(command), flush=True)
    if not dry_run:
        subprocess.run(command, check=True, env=env)


def main() -> int:
    args = parse_args()
    repo = Path(args.repo_root).expanduser().resolve()
    analysis = Path(args.analysis_root).expanduser()
    if not analysis.is_absolute():
        analysis = (repo / analysis).resolve()
    else:
        analysis = analysis.resolve()
    data = Path(args.data_root).expanduser().resolve() if args.data_root else None
    external = Path(args.external_root).expanduser().resolve() if args.external_root else None

    preflight_only = args.stage == ["preflight"]
    stages = STAGE_ORDER if "all" in args.stage else [stage for stage in args.stage if stage != "preflight"]
    check_stages = STAGE_ORDER[:-1] if preflight_only else stages

    missing = [path for path in required_files(repo, data, external, check_stages) if not path.is_file()]
    if missing:
        print("Preflight failed; missing required file(s):", file=sys.stderr)
        for path in missing:
            print(f"  - {path}", file=sys.stderr)
        return 2
    if any(stage in check_stages for stage in ["cptac", "tcga", "apollo"]) and data is None:
        print("--data-root is required for the selected stages", file=sys.stderr)
        return 2
    if any(stage in check_stages for stage in ["apollo", "threeca", "h5ad"]) and external is None:
        print("--external-root is required for the selected stages", file=sys.stderr)
        return 2

    if preflight_only:
        print("Preflight passed for all scientific-analysis stages.")
        return 0

    analysis.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    package_path = str(repo / "analysis/src")
    env["PYTHONPATH"] = package_path + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["ARTICLE002_DATA_ROOT"] = str(data or repo)
    protocol = repo / "analysis/definitions/BIO-05_V2_LOCKED_PROTOCOL.json"
    signatures = repo / "analysis/definitions/signature_gene_sets_v0_1.csv"
    executed = repo / "analysis/code/current_p2_closure/executed_code"

    if "cptac" in stages:
        run(
            [
                args.python,
                str(repo / "analysis/code/run_cptac_v2.py"),
                "--analysis-root", str(analysis),
                "--data-root", str(data),
                "--protocol", str(protocol),
                "--signature-table", str(signatures),
            ],
            env,
            args.dry_run,
        )

    if "tcga" in stages:
        run(
            [
                args.python,
                str(executed / "run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py"),
                "--project-root", str(data),
                "--out-dir", str(analysis / "results/tcga"),
                "--figure-dir", str(analysis / "figures/tcga"),
                "--report", str(analysis / "reports/tcga_xena_sensitivity_v0_5_report.md"),
                "--summary", str(analysis / "results/tcga/tcga_source_execution_summary.json"),
            ],
            env,
            args.dry_run,
        )

    if "apollo" in stages:
        discovery_corr = (
            analysis / "results/cptac/cptac_primary_concordance_registry_v2.csv"
            if "cptac" in stages or (analysis / "results/cptac/cptac_primary_concordance_registry_v2.csv").is_file()
            else data / "04_processed/discovery_v0_4/bcm_matched_correlation_summary_v0_4.csv"
        )
        run(
            [
                args.python,
                str(executed / "run_bio05_apollo_external_cohort_v0_7.py"),
                "--project-root", str(analysis),
                "--input-dir", str(external / "proteogenomics"),
                "--signature-table", str(signatures),
                "--discovery-corr", str(discovery_corr),
                "--out-dir", str(analysis / "results/apollo"),
                "--report", str(analysis / "reports/public_apollo_external_cohort_v0_7_report.md"),
            ],
            env,
            args.dry_run,
        )

    if "threeca" in stages:
        run(
            [
                args.python,
                str(executed / "run_bio05_3ca_source_mapping_v0_8.py"),
                "--project-root", str(analysis),
                "--input-dir", str(external / "3ca"),
                "--signature-table", str(signatures),
                "--out-dir", str(analysis / "results/threeca"),
                "--report", str(analysis / "reports/public_3ca_source_mapping_v0_8_report.md"),
            ],
            env,
            args.dry_run,
        )

    if "h5ad" in stages:
        command = [
            args.python,
            str(repo / "analysis/code/run_h5ad_support_v2.py"),
            "--mode", "all",
            "--scrna-queue", str(repo / "analysis/input_manifests/public_validation_queue.tsv"),
            "--spatial-queue", str(repo / "analysis/input_manifests/brca_spatial_queue.tsv"),
            "--signature-table", str(signatures),
            "--analysis-root", str(analysis),
            "--input-root", str(external),
        ]
        if args.source_cache_dir:
            command.extend(["--source-cache-dir", str(Path(args.source_cache_dir).expanduser().resolve())])
        run(command, env, args.dry_run)

    if "postprocess" in stages:
        run(
            [args.python, str(repo / "analysis/code/postprocess_support_v2.py"), "--analysis-root", str(analysis)],
            env,
            args.dry_run,
        )

    if "compare" in stages:
        run(
            [
                args.python,
                str(repo / "analysis/compare_recomputed_to_locked.py"),
                "--analysis-root", str(analysis),
                "--repo-root", str(repo),
            ],
            env,
            args.dry_run,
        )

    print(f"Completed stages: {', '.join(stages)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
