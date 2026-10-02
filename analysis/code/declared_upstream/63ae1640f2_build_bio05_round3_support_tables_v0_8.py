#!/usr/bin/env python3
from __future__ import annotations

import os as _workspace_os


def _workspace_path(value: str) -> str:
    """Resolve portable workspace markers without embedding a local user path."""
    mapping = {
        "${DATA_WORKSPACE}": "DATA_WORKSPACE",
        "${RESEARCH_WORKSPACE}": "RESEARCH_WORKSPACE",
        "${HOME}": "HOME",
    }
    resolved = value
    for marker, variable in mapping.items():
        if marker not in resolved:
            continue
        root = _workspace_os.environ.get(variable, "").rstrip("/")
        if not root:
            raise RuntimeError(f"Set {variable} before running this script")
        resolved = resolved.replace(marker, root)
    return resolved


import importlib.util
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(
    _workspace_path('${DATA_WORKSPACE}/14_pan_cancer_proteogenomic_immune_evasion')
)
VERSION = "v0_8_round3"
MIN_SIGNATURE_GENES = 3
APOLLO_BOOTSTRAP_N = 1000
APOLLO_RANDOM_SEED = 20260730


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


APOLLO = load_module(PROJECT_ROOT / "scripts/run_bio05_apollo_external_cohort_v0_7.py", "apollo_v07")


OUT_DIR = PROJECT_ROOT / "04_processed/discovery_v0_8_public_validation/round3_support"


def fisher_ci(rho: float, n: int) -> tuple[float | None, float | None]:
    if not np.isfinite(rho) or n <= 3 or abs(rho) >= 1:
        return None, None
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3)
    low = np.tanh(z - 1.96 * se)
    high = np.tanh(z + 1.96 * se)
    return float(low), float(high)


def norm_public_id(value: Any) -> str:
    text = str(value).upper().strip()
    m = re.match(r"^(CID\d+)[A-Z]$", text)
    if m:
        return m.group(1)
    return text


def build_cptac_within_cancer_table() -> pd.DataFrame:
    corr = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_4/bcm_matched_correlation_summary_v0_4.csv")
    corr = corr[
        corr["comparison_scope"].eq("within_cancer")
        & corr["comparison"].isin(["rna_vs_proteomics", "rna_vs_phosphoproteomics"])
    ].copy()
    corr["ci_lower"] = np.nan
    corr["ci_upper"] = np.nan
    for idx, row in corr.iterrows():
        ci_low, ci_high = fisher_ci(float(row["spearman_rho"]), int(row["n"]))
        corr.loc[idx, "ci_lower"] = ci_low
        corr.loc[idx, "ci_upper"] = ci_high
    corr["status"] = "analyzable"

    cancers = sorted(corr["cancer"].dropna().unique().tolist())
    layers = ["proteomics", "phosphoproteomics"]
    sigs = (
        corr[["signature", "signature_label"]]
        .drop_duplicates()
        .sort_values(["signature_label", "signature"])
        .to_dict("records")
    )
    grid_rows = []
    for cancer in cancers:
        for layer in layers:
            for sig in sigs:
                grid_rows.append(
                    {
                        "cancer": cancer,
                        "comparison_layer": layer,
                        "signature": sig["signature"],
                        "signature_label": sig["signature_label"],
                    }
                )
    grid = pd.DataFrame(grid_rows)
    out = grid.merge(
        corr[
            [
                "comparison",
                "comparison_layer",
                "cancer",
                "signature",
                "signature_label",
                "n",
                "spearman_rho",
                "p_value",
                "fdr_family",
                "fdr_bh",
                "ci_lower",
                "ci_upper",
                "status",
            ]
        ],
        on=["comparison_layer", "cancer", "signature", "signature_label"],
        how="left",
    )
    out["comparison"] = out["comparison"].fillna(
        out["comparison_layer"].map(
            {
                "proteomics": "rna_vs_proteomics",
                "phosphoproteomics": "rna_vs_phosphoproteomics",
            }
        )
    )
    out["status"] = out["status"].fillna("not_analyzable")
    return out.sort_values(["comparison_layer", "cancer", "signature_label"]).reset_index(drop=True)


def build_public_overlap_audit() -> pd.DataFrame:
    group_df = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_8_public_validation/public_source_group_counts_v0_8.csv")
    public_units = (
        group_df[["cancer_code", "source_id", "dataset_title", "donor_id_resolved"]]
        .drop_duplicates()
        .rename(columns={"donor_id_resolved": "public_unit_id"})
    )

    rows: list[dict[str, Any]] = []
    for cancer, sub in public_units.groupby("cancer_code", dropna=False):
        datasets = sub[["source_id", "dataset_title"]].drop_duplicates().to_dict("records")
        if len(datasets) < 2:
            continue
        for i in range(len(datasets)):
            for j in range(i + 1, len(datasets)):
                left = datasets[i]
                right = datasets[j]
                left_ids = sorted(
                    sub.loc[sub["source_id"].eq(left["source_id"]), "public_unit_id"].astype(str).unique().tolist()
                )
                right_ids = sorted(
                    sub.loc[sub["source_id"].eq(right["source_id"]), "public_unit_id"].astype(str).unique().tolist()
                )
                exact = sorted(set(left_ids) & set(right_ids))
                root = sorted(set(map(norm_public_id, left_ids)) & set(map(norm_public_id, right_ids)))
                rows.append(
                    {
                        "audit_scope": "same_cancer_public_scrna_pair",
                        "cancer_code": cancer,
                        "left_source_id": left["source_id"],
                        "left_dataset_title": left["dataset_title"],
                        "right_source_id": right["source_id"],
                        "right_dataset_title": right["dataset_title"],
                        "left_unit_type": "donor",
                        "right_unit_type": "donor",
                        "left_unit_n": len(left_ids),
                        "right_unit_n": len(right_ids),
                        "exact_overlap_n": len(exact),
                        "exact_overlap_examples": ";".join(exact[:10]),
                        "root_overlap_n": len(root),
                        "root_overlap_examples": ";".join(root[:10]),
                        "same_public_collection": "no",
                        "direct_id_comparable": "yes",
                        "overlap_test_status": "exact_public_donor_id_audit",
                        "evidence_label_after_audit": (
                            "separate_public_dataset_no_exact_donor_overlap"
                            if len(exact) == 0
                            else "same_or_overlapping_public_dataset"
                        ),
                        "note": "public donor IDs directly comparable within same-cancer source-mapping layer",
                    }
                )

    spatial = pd.read_csv(
        PROJECT_ROOT / "04_processed/discovery_v0_7_public_validation/spatial_brca_sample_signature_summary_v0_7.csv"
    )
    spatial_ids = sorted(spatial["sample_id"].astype(str).unique().tolist())
    brca_sc = public_units[public_units["cancer_code"].eq("brca")].copy()
    for source_id, same_collection in [
        ("dea97145-f712-431c-a223-6b5f565f362a", "yes"),
        ("3f7c572c-cd73-4b51-a313-207c7f20f188", "no"),
    ]:
        sub = brca_sc[brca_sc["source_id"].eq(source_id)].copy()
        if sub.empty:
            continue
        donor_ids = sorted(sub["public_unit_id"].astype(str).unique().tolist())
        exact = sorted(set(map(str.upper, donor_ids)) & set(map(str.upper, spatial_ids)))
        root = sorted(set(map(norm_public_id, donor_ids)) & set(map(norm_public_id, spatial_ids)))
        rows.append(
            {
                "audit_scope": "brca_scrna_vs_spatial",
                "cancer_code": "brca",
                "left_source_id": source_id,
                "left_dataset_title": sub["dataset_title"].iloc[0],
                "right_source_id": "BRCA_VISIUM_v0_7",
                "right_dataset_title": "Four official BRCA Visium samples",
                "left_unit_type": "donor",
                "right_unit_type": "sample",
                "left_unit_n": len(donor_ids),
                "right_unit_n": len(spatial_ids),
                "exact_overlap_n": len(exact),
                "exact_overlap_examples": ";".join(exact[:10]),
                "root_overlap_n": len(root),
                "root_overlap_examples": ";".join(root[:10]),
                "same_public_collection": same_collection,
                "direct_id_comparable": "yes",
                "overlap_test_status": "spatial_sample_vs_scrna_donor_id_audit",
                "evidence_label_after_audit": (
                    "same_collection_spatial_triangulation"
                    if same_collection == "yes"
                    else "separate_public_dataset"
                ),
                "note": (
                    "same atlas collection; not independent spatial replication"
                    if same_collection == "yes"
                    else "different public scRNA route; no exact overlap detected"
                ),
            }
        )

    discovery = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_4/bcm_matched_rna_protein_phosphoprotein_scores_v0_4.csv")
    apollo_scores = pd.read_csv(
        PROJECT_ROOT / "04_processed/discovery_v0_7_public_validation/apollo_signature_scores_v0_7.csv"
    )
    for cancer, cohort in [("luad", "APOLLO-LUAD"), ("ov", "APOLLO-OV")]:
        cptac_cases = sorted(discovery.loc[discovery["cancer"].eq(cancer), "normalized_case_id"].astype(str).unique())
        cptac_protein_samples = sorted(
            discovery.loc[discovery["cancer"].eq(cancer), "protein_layer_sample_id"].astype(str).unique()
        )
        apollo_ids = sorted(
            apollo_scores.loc[
                apollo_scores["cohort_label"].eq(cohort) & apollo_scores["comparison_layer"].eq("rna"),
                "sample_id",
            ]
            .astype(str)
            .unique()
            .tolist()
        )
        exact = sorted(set(cptac_cases) & set(apollo_ids))
        exact_protein = sorted(set(cptac_protein_samples) & set(apollo_ids))
        rows.append(
            {
                "audit_scope": "apollo_vs_cptac_case_overlap",
                "cancer_code": cancer,
                "left_source_id": f"CPTAC_{cancer.upper()}_discovery",
                "left_dataset_title": "CPTAC matched discovery",
                "right_source_id": cohort,
                "right_dataset_title": cohort,
                "left_unit_type": "case",
                "right_unit_type": "case",
                "left_unit_n": len(cptac_cases),
                "right_unit_n": len(apollo_ids),
                "exact_overlap_n": len(exact),
                "exact_overlap_examples": ";".join(exact[:10]),
                "root_overlap_n": len(exact_protein),
                "root_overlap_examples": ";".join(exact_protein[:10]),
                "same_public_collection": "no",
                "direct_id_comparable": "limited_no_harmonized_crosswalk",
                "overlap_test_status": "exact_public_id_check_only",
                "evidence_label_after_audit": "separate_public_cohort_by_provenance",
                "note": "exact public IDs do not overlap, but donor-level harmonized crosswalk is unavailable; use support rather than proven independent replication",
            }
        )

    return pd.DataFrame(rows)


def build_apollo_complete_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    signature_genes = APOLLO.load_signature_genes(PROJECT_ROOT / "04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv")
    score_df = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_7_public_validation/apollo_signature_scores_v0_7.csv")
    recovery_df = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_7_public_validation/apollo_signature_recovery_v0_7.csv")
    discovery_corr = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_4/bcm_matched_correlation_summary_v0_4.csv")
    discovery_corr = discovery_corr[
        discovery_corr["comparison_scope"].eq("within_cancer")
        & discovery_corr["comparison"].isin(["rna_vs_proteomics", "rna_vs_phosphoproteomics"])
    ].copy()
    discovery_corr["comparison_layer"] = discovery_corr["comparison_layer"].astype(str)

    signature_meta = []
    for sig, genes in signature_genes.items():
        signature_meta.append(
            {
                "signature": sig,
                "signature_label": APOLLO.SIGNATURE_LABELS.get(sig, sig),
                "requested_gene_count": len(genes),
            }
        )
    signature_meta_df = pd.DataFrame(signature_meta)
    cohorts = [
        ("APOLLO-LUAD", "luad", "proteomics"),
        ("APOLLO-LUAD", "luad", "phosphoproteomics"),
        ("APOLLO-OV", "ov", "proteomics"),
    ]
    flow_rows: list[dict[str, Any]] = []
    result_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(APOLLO_RANDOM_SEED)

    for cohort_label, cancer_code, layer in cohorts:
        rna_scores = score_df[
            score_df["cohort_label"].eq(cohort_label)
            & score_df["cancer_code"].eq(cancer_code)
            & score_df["comparison_layer"].eq("rna")
        ].copy()
        omic_scores = score_df[
            score_df["cohort_label"].eq(cohort_label)
            & score_df["cancer_code"].eq(cancer_code)
            & score_df["comparison_layer"].eq(layer)
        ].copy()
        rna_samples = sorted(rna_scores["sample_id"].astype(str).unique().tolist())
        omic_samples = sorted(omic_scores["sample_id"].astype(str).unique().tolist())
        matched_total = sorted(set(rna_samples) & set(omic_samples))

        rna_recovery = recovery_df[
            recovery_df["cohort_label"].eq(cohort_label) & recovery_df["comparison_layer"].eq("rna")
        ][["signature", "recovered_gene_count"]].rename(columns={"recovered_gene_count": "rna_recovered_gene_count"})
        omic_recovery = recovery_df[
            recovery_df["cohort_label"].eq(cohort_label) & recovery_df["comparison_layer"].eq(layer)
        ][["signature", "recovered_gene_count"]].rename(columns={"recovered_gene_count": "omic_recovered_feature_count"})

        merged_scores = rna_scores.merge(
            omic_scores,
            on=["signature", "signature_label", "sample_id"],
            how="outer",
            suffixes=("_rna", "_omic"),
        )

        layer_rows = []
        for sig_row in signature_meta_df.itertuples(index=False):
            sub = merged_scores[merged_scores["signature"].eq(sig_row.signature)].copy()
            x = pd.to_numeric(sub.get("signature_score_rna"), errors="coerce").to_numpy(dtype=float)
            y = pd.to_numeric(sub.get("signature_score_omic"), errors="coerce").to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            x = x[mask]
            y = y[mask]
            matched_n = int(mask.sum())
            omic_recovered = omic_recovery.loc[
                omic_recovery["signature"].eq(sig_row.signature), "omic_recovered_feature_count"
            ]
            omic_recovered_n = int(omic_recovered.iloc[0]) if not omic_recovered.empty else 0
            rna_recovered = rna_recovery.loc[
                rna_recovery["signature"].eq(sig_row.signature), "rna_recovered_gene_count"
            ]
            rna_recovered_n = int(rna_recovered.iloc[0]) if not rna_recovered.empty else 0

            row = {
                "cohort_label": cohort_label,
                "cancer_code": cancer_code,
                "comparison_layer": layer,
                "feature_space": (
                    "gene_collapsed_phosphosite_primary_gene"
                    if layer == "phosphoproteomics"
                    else "gene"
                ),
                "signature": sig_row.signature,
                "signature_label": sig_row.signature_label,
                "requested_gene_count": int(sig_row.requested_gene_count),
                "rna_recovered_gene_count": rna_recovered_n,
                "omic_recovered_feature_count": omic_recovered_n,
                "total_rna_samples": len(rna_samples),
                "total_omic_samples": len(omic_samples),
                "matched_sample_n_total": len(matched_total),
                "matched_analyzable_n": matched_n,
                "spearman_rho": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
                "spearman_p": np.nan,
                "bh_fdr_within_cohort_layer": np.nan,
                "direction_label": "",
                "status": "",
                "status_reason": "",
                "discovery_rho": np.nan,
                "sign_agreement_vs_discovery": np.nan,
            }
            disc = discovery_corr[
                discovery_corr["cancer"].eq(cancer_code)
                & discovery_corr["comparison_layer"].eq(layer)
                & discovery_corr["signature"].eq(sig_row.signature)
            ]
            if not disc.empty:
                row["discovery_rho"] = float(disc["spearman_rho"].iloc[0])

            if omic_recovered_n < MIN_SIGNATURE_GENES:
                row["status"] = "not_evaluable"
                row["status_reason"] = f"omic recovered feature count < {MIN_SIGNATURE_GENES}"
            elif matched_n < 10:
                row["status"] = "not_evaluable"
                row["status_reason"] = "matched analyzable n < 10"
            else:
                rho, ci_low, ci_high = APOLLO.bootstrap_spearman_ci(
                    x, y, bootstrap_n=APOLLO_BOOTSTRAP_N, rng=rng
                )
                p_value = float(stats.spearmanr(x, y, nan_policy="omit").pvalue)
                row["spearman_rho"] = float(rho)
                row["ci_lower"] = float(ci_low)
                row["ci_upper"] = float(ci_high)
                row["spearman_p"] = p_value
                row["direction_label"] = "positive" if rho > 0 else "negative" if rho < 0 else "zero"
                if np.isfinite(row["discovery_rho"]):
                    row["sign_agreement_vs_discovery"] = float(np.sign(rho) == np.sign(row["discovery_rho"]))
                layer_rows.append(row)
            result_rows.append(row)

        if layer_rows:
            layer_df = pd.DataFrame(layer_rows)
            fdr = stats.false_discovery_control(layer_df["spearman_p"].to_numpy(dtype=float), method="bh")
            for sig, value in zip(layer_df["signature"].tolist(), fdr):
                mask = (
                    pd.Series([r["cohort_label"] == cohort_label and r["comparison_layer"] == layer and r["signature"] == sig for r in result_rows])
                )
                idx = mask[mask].index[0]
                result_rows[idx]["bh_fdr_within_cohort_layer"] = float(value)
            for row in result_rows:
                if row["cohort_label"] == cohort_label and row["comparison_layer"] == layer and row["status"] == "":
                    fdr_val = float(row["bh_fdr_within_cohort_layer"])
                    rho = float(row["spearman_rho"])
                    if rho > 0 and fdr_val < 0.05:
                        row["status"] = "supported"
                        row["status_reason"] = "positive rho with cohort-layer BH-FDR < 0.05"
                    elif rho < 0:
                        row["status"] = "opposite"
                        row["status_reason"] = "negative-direction concordance relative to positive discovery expectation"
                    else:
                        row["status"] = "null_result"
                        row["status_reason"] = "positive or zero rho without cohort-layer BH-FDR support"

        result_subset = [
            r
            for r in result_rows
            if r["cohort_label"] == cohort_label and r["comparison_layer"] == layer
        ]
        flow_rows.append(
            {
                "cohort_label": cohort_label,
                "cancer_code": cancer_code,
                "comparison_layer": layer,
                "total_rna_samples": len(rna_samples),
                "total_omic_samples": len(omic_samples),
                "matched_sample_n_total": len(matched_total),
                "evaluable_signature_n": int(sum(r["status"] != "not_evaluable" for r in result_subset)),
                "not_evaluable_signature_n": int(sum(r["status"] == "not_evaluable" for r in result_subset)),
                "not_evaluable_signatures": ";".join(
                    sorted(r["signature_label"] for r in result_subset if r["status"] == "not_evaluable")
                ),
            }
        )

    result_df = pd.DataFrame(result_rows)
    result_df["status"] = result_df["status"].fillna("")
    result_df["status_reason"] = result_df["status_reason"].fillna("")
    evaluable_mask = result_df["status"].eq("")
    supported_mask = (
        evaluable_mask
        & result_df["spearman_rho"].gt(0)
        & result_df["bh_fdr_within_cohort_layer"].lt(0.05)
    )
    opposite_mask = evaluable_mask & result_df["spearman_rho"].lt(0)
    null_mask = evaluable_mask & ~(supported_mask | opposite_mask)

    result_df.loc[supported_mask, "status"] = "supported"
    result_df.loc[supported_mask, "status_reason"] = "positive rho with cohort-layer BH-FDR < 0.05"
    result_df.loc[opposite_mask, "status"] = "opposite"
    result_df.loc[opposite_mask, "status_reason"] = (
        "negative-direction concordance relative to positive discovery expectation"
    )
    result_df.loc[null_mask, "status"] = "null_result"
    result_df.loc[null_mask, "status_reason"] = (
        "positive or zero rho without cohort-layer BH-FDR support"
    )

    return result_df, pd.DataFrame(flow_rows)


def build_public_inference_claim_window() -> pd.DataFrame:
    inf = pd.read_csv(PROJECT_ROOT / "04_processed/discovery_v0_8_public_validation/public_source_pairwise_inference_v0_8.csv")
    rename = {
        "two_sided_sign_p": "sign_test_P",
        "two_sided_wilcoxon_p": "Wilcoxon_P",
        "sign_fdr_global": "global_sign_FDR",
        "wilcoxon_fdr_global": "global_Wilcoxon_FDR",
        "sign_fdr_within_source_cancer": "within_source_sign_FDR",
        "wilcoxon_fdr_within_source_cancer": "within_source_Wilcoxon_FDR",
    }
    out = inf.rename(columns=rename).copy()
    out["claim_window_n6"] = out["paired_donor_n"].astype(int) >= 6
    out["main_text_claim_used"] = (
        out["claim_window_n6"]
        & out["global_Wilcoxon_FDR"].astype(float).lt(0.05)
    )
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cptac_within = build_cptac_within_cancer_table()
    cptac_within.to_csv(OUT_DIR / f"cptac_within_cancer_concordance_full_{VERSION}.csv", index=False)

    overlap = build_public_overlap_audit()
    overlap.to_csv(OUT_DIR / f"public_source_overlap_audit_{VERSION}.csv", index=False)

    apollo_complete, apollo_flow = build_apollo_complete_tables()
    apollo_complete.to_csv(OUT_DIR / f"apollo_complete_concordance_{VERSION}.csv", index=False)
    apollo_flow.to_csv(OUT_DIR / f"apollo_sample_flow_{VERSION}.csv", index=False)

    public_inf = build_public_inference_claim_window()
    public_inf.to_csv(OUT_DIR / f"public_source_pairwise_inference_claim_window_{VERSION}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
