#!/usr/bin/env python3
"""Compare recomputed scientific registries with the locked release tables."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Comparison:
    name: str
    recomputed: str
    locked: str
    keys: tuple[str, ...]
    rename: dict[str, str] = field(default_factory=dict)
    source_base: str = "analysis"


COMPARISONS = [
    Comparison("CPTAC primary", "results/cptac/cptac_primary_concordance_registry_v2.csv", "tables/supplementary/S3_Primary.csv", ("cancer", "program", "comparison_layer")),
    Comparison("CPTAC feature sensitivity", "results/cptac/cptac_feature_row_sensitivity_registry_v2.csv", "tables/supplementary/S4a_FeatureSens.csv", ("cancer", "program", "comparison_layer")),
    Comparison("CPTAC threshold sensitivity", "results/cptac/cptac_threshold1_sensitivity_registry_v2.csv", "tables/supplementary/S4b_Threshold1.csv", ("cancer", "program", "comparison_layer")),
    Comparison("CPTAC pooled sensitivity", "results/cptac/cptac_within_cancer_standardized_pooled_sensitivity_v2.csv", "tables/supplementary/S4c_Pooled.csv", ("comparison_layer", "program")),
    Comparison("CPTAC residual diagnostics", "results/cptac/cptac_residual_model_diagnostics_v2.csv", "tables/supplementary/S4d_ResidualDiag.csv", ("cancer", "program", "comparison_layer")),
    Comparison("CPTAC residual context", "results/cptac/cptac_residual_context_registry_v2.csv", "tables/supplementary/S5_ResidualContext.csv", ("cancer", "program", "comparison_layer", "context_variable")),
    Comparison(
        "TCGA full atlas",
        "results/tcga/tcga_xena_68immune_within_cancer_correlations_v0_5.csv",
        "tables/supplementary/S6a_TCGA_Full.csv",
        ("cancer_type", "program", "correlated_program_or_gene"),
        {
            "signature": "program",
            "signature_label": "program_label",
            "cancer": "cancer_type",
            "immune_signature_68": "correlated_program_or_gene",
            "n": "sample_size_n",
            "fdr_bh": "q_value_bh",
        },
    ),
    Comparison("TCGA program summary", "results/tcga/tcga_xena_signature_summary_v0_5.csv", "tables/supplementary/S6b_TCGA_Summary.csv", ("signature",)),
    Comparison("TCGA fixed anchors", "analysis/definitions/TCGA_ANCHOR_SELECTION_CURRENT_P2.csv", "tables/supplementary/S6c_TCGA_Anchors.csv", ("signature", "immune_signature_68"), source_base="repo"),
    Comparison("TCGA immune subtypes", "results/tcga/tcga_xena_immune_subtype_summary_v0_5.csv", "tables/supplementary/S6d_TCGA_Subtypes.csv", ("signature", "immune_subtype")),
    Comparison("APOLLO", "results/apollo/apollo_concordance_registry_v2.csv", "tables/supplementary/S7_APOLLO.csv", ("cohort_label", "cancer_code", "comparison_layer", "program")),
    Comparison("single-cell global", "results/scrna_combined/public_scrna_global_wilcoxon_registry_v2.csv", "tables/supplementary/S8a_scRNA_Global.csv", ("program", "contrast")),
    Comparison("single-cell source", "results/scrna_combined/public_scrna_within_source_wilcoxon_registry_v2.csv", "tables/supplementary/S8b_scRNA_Source.csv", ("source_id", "cancer_code", "program", "contrast")),
    Comparison("spatial BRCA", "results/spatial/spatial_brca_pairwise_contrasts_v0_7.csv", "tables/supplementary/S9_Spatial_BRCA.csv", ("signature_label", "contrast")),
]


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument("--repo-root", default=str(repo))
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--output", default="")
    return parser.parse_args()


def normalized_text(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.replace(r"^(true|false)$", lambda match: match.group(1).lower(), regex=True)
    )


def numeric_series(series: pd.Series) -> tuple[bool, pd.Series]:
    nonempty = normalized_text(series).ne("")
    parsed = pd.to_numeric(series, errors="coerce")
    return bool(nonempty.any() and parsed[nonempty].notna().all()), parsed


def compare_frames(spec: Comparison, observed: pd.DataFrame, expected: pd.DataFrame, rtol: float, atol: float) -> dict[str, object]:
    observed = observed.rename(columns=spec.rename)
    missing_columns = sorted(set(expected.columns) - set(observed.columns))
    if missing_columns:
        return {"status": "FAIL", "reason": "missing_columns", "missing_columns": missing_columns}
    if any(key not in expected.columns or key not in observed.columns for key in spec.keys):
        return {"status": "FAIL", "reason": "missing_key_column", "keys": list(spec.keys)}
    if expected.duplicated(list(spec.keys)).any() or observed.duplicated(list(spec.keys)).any():
        return {"status": "FAIL", "reason": "duplicate_keys", "keys": list(spec.keys)}

    expected = expected.sort_values(list(spec.keys), kind="mergesort").reset_index(drop=True)
    observed = observed.sort_values(list(spec.keys), kind="mergesort").reset_index(drop=True)
    expected_keys = expected[list(spec.keys)].fillna("").astype(str)
    observed_keys = observed[list(spec.keys)].fillna("").astype(str)
    if expected_keys.shape != observed_keys.shape or not expected_keys.equals(observed_keys):
        return {
            "status": "FAIL",
            "reason": "key_set_mismatch",
            "expected_rows": int(len(expected)),
            "observed_rows": int(len(observed)),
        }

    mismatches: list[dict[str, object]] = []
    for column in expected.columns:
        left_is_numeric, left = numeric_series(expected[column])
        right_is_numeric, right = numeric_series(observed[column])
        if left_is_numeric and right_is_numeric:
            both_missing = left.isna() & right.isna()
            equal = np.isclose(left.fillna(0).to_numpy(float), right.fillna(0).to_numpy(float), rtol=rtol, atol=atol) | both_missing.to_numpy()
        else:
            equal = normalized_text(expected[column]).eq(normalized_text(observed[column])).to_numpy()
        if not bool(np.all(equal)):
            rows = np.flatnonzero(~equal)[:5]
            mismatches.append(
                {
                    "column": column,
                    "count": int(np.count_nonzero(~equal)),
                    "examples": [
                        {
                            "key": {key: str(expected.loc[index, key]) for key in spec.keys},
                            "locked": None if pd.isna(expected.loc[index, column]) else str(expected.loc[index, column]),
                            "recomputed": None if pd.isna(observed.loc[index, column]) else str(observed.loc[index, column]),
                        }
                        for index in rows
                    ],
                }
            )
    return {
        "status": "PASS" if not mismatches else "FAIL",
        "expected_rows": int(len(expected)),
        "observed_rows": int(len(observed)),
        "compared_columns": int(len(expected.columns)),
        "mismatches": mismatches,
    }


def main() -> int:
    args = parse_args()
    analysis = Path(args.analysis_root).expanduser().resolve()
    repo = Path(args.repo_root).expanduser().resolve()
    output = Path(args.output).expanduser().resolve() if args.output else analysis / "validation/analysis_parity.json"
    rows: list[dict[str, object]] = []
    for spec in COMPARISONS:
        observed_path = (repo if spec.source_base == "repo" else analysis) / spec.recomputed
        expected_path = repo / spec.locked
        if not observed_path.is_file():
            result: dict[str, object] = {"status": "FAIL", "reason": "recomputed_file_missing"}
        elif not expected_path.is_file():
            result = {"status": "FAIL", "reason": "locked_file_missing"}
        else:
            result = compare_frames(
                spec,
                pd.read_csv(observed_path),
                pd.read_csv(expected_path),
                args.rtol,
                args.atol,
            )
        rows.append(
            {
                "name": spec.name,
                "recomputed": spec.recomputed,
                "locked": spec.locked,
                **result,
            }
        )
    passed = sum(row["status"] == "PASS" for row in rows)
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if passed == len(rows) else "FAIL",
        "comparisons_passed": passed,
        "comparisons_total": len(rows),
        "rtol": args.rtol,
        "atol": args.atol,
        "comparisons": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: payload[key] for key in ["status", "comparisons_passed", "comparisons_total"]}, indent=2))
    print(f"Report: {output}")
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
