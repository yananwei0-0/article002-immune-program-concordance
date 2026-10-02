#!/usr/bin/env python3
"""Post-process real-source BIO-05 support layers under locked multiplicity rules."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


PROGRAMS = [
    "antigen_presentation_mhc_i",
    "antigen_presentation_mhc_ii",
    "ifn_gamma_response",
    "cytolytic_t_cell_context",
    "checkpoint_exhaustion_context",
    "myeloid_inflammatory_context",
    "tgfb_emt_exclusion_context",
    "proteasome_antigen_processing",
]
LABELS = {
    "antigen_presentation_mhc_i": "MHC-I",
    "antigen_presentation_mhc_ii": "MHC-II",
    "ifn_gamma_response": "IFN-gamma",
    "cytolytic_t_cell_context": "Cytolytic",
    "checkpoint_exhaustion_context": "Checkpoint",
    "myeloid_inflammatory_context": "Myeloid",
    "tgfb_emt_exclusion_context": "TGF/EMT",
    "proteasome_antigen_processing": "Proteasome",
}
CONTRASTS = [
    ("malignant", "immune"),
    ("malignant", "stromal"),
    ("immune", "stromal"),
]
SOURCE_PUBLICATION_IDENTITY = {
    "3f7c572c-cd73-4b51-a313-207c7f20f188": "doi:10.1038/s41467-024-49916-4; aggregated multi-study TME atlas",
    "dea97145-f712-431c-a223-6b5f565f362a": "doi:10.1038/s41588-021-00911-1",
    "1df8c90d-d299-4b2e-a54d-a5a80f36e780": "doi:10.1073/pnas.2103240118",
    "3f50314f-bdc9-40c6-8e4a-b0901ebfbe4c": "doi:10.1016/j.ccell.2021.03.007",
    "4d82bd8e-8827-410b-8fc7-685b7dd92585": "doi:10.1016/j.cell.2021.08.003",
    "0bebef1a-4607-4584-9070-dacf89a0d635": "doi:10.1186/s40164-025-00740-6",
    "bc7397a3-ea49-4d57-84b8-80bd6885d4c4": "CELLxGENE collection; publication DOI absent from local H5AD citation",
    "999f2a15-3d7e-440b-96ae-2c806799c08c": "doi:10.1101/2022.08.27.505439; aggregated GBmap",
    "3ca-peng2019-pdac": "3CA official Peng 2019 pancreas route",
    "3ca-regner2021-ucec": "3CA official Regner 2021 endometrial subset route",
}
SOURCE_DONOR_FIELD_NOTE = {
    "3f50314f-bdc9-40c6-8e4a-b0901ebfbe4c": "six public donor labels map one-to-one to six hca_data_portal_donor_uuid values in the retained H5AD",
    "3ca-regner2021-ucec": "sample identifier used as source-specific donor proxy because a separate patient field was unavailable",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True)
    return parser.parse_args()


def bh_fdr(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out


def spearman_ci(rho: float, n: int) -> tuple[float, float]:
    if not np.isfinite(rho) or n <= 3:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1 / math.sqrt(n - 3)
    return float(np.tanh(z - 1.9599639845 * se)), float(np.tanh(z + 1.9599639845 * se))


def apollo_registry(root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    source_dir = root / "results" / "apollo"
    scores = pd.read_csv(source_dir / "apollo_signature_scores_v0_7.csv")
    key = ["cohort_label", "cancer_code", "comparison_layer", "signature", "sample_id"]
    duplicate_n = int(scores.duplicated(key, keep=False).sum())
    if duplicate_n:
        raise RuntimeError(f"APOLLO score table contains {duplicate_n} duplicate inference keys")
    strata = [
        ("APOLLO-LUAD", "luad", "proteomics"),
        ("APOLLO-LUAD", "luad", "phosphoproteomics"),
        ("APOLLO-OV", "ov", "proteomics"),
    ]
    rows: list[dict[str, Any]] = []
    for cohort, cancer, layer in strata:
        rna = scores[
            scores["cohort_label"].eq(cohort)
            & scores["cancer_code"].eq(cancer)
            & scores["comparison_layer"].eq("rna")
        ]
        omic = scores[
            scores["cohort_label"].eq(cohort)
            & scores["cancer_code"].eq(cancer)
            & scores["comparison_layer"].eq(layer)
        ]
        for program in PROGRAMS:
            left = rna[rna["signature"].eq(program)][["sample_id", "signature_score"]].rename(columns={"signature_score": "rna_score"})
            right = omic[omic["signature"].eq(program)][["sample_id", "signature_score"]].rename(columns={"signature_score": "omic_score"})
            if left.duplicated("sample_id").any() or right.duplicated("sample_id").any():
                raise RuntimeError(f"APOLLO sample uniqueness failed for {cohort}/{layer}/{program}")
            matched = left.merge(right, on="sample_id", how="inner", validate="one_to_one").dropna()
            n = len(matched)
            if n >= 10 and matched["rna_score"].nunique() >= 3 and matched["omic_score"].nunique() >= 3:
                stat = spearmanr(matched["rna_score"], matched["omic_score"])
                rho, p = float(stat.statistic), float(stat.pvalue)
                ci_low, ci_high = spearman_ci(rho, n)
                eval_status, reason = "evaluable", ""
            else:
                rho = p = ci_low = ci_high = np.nan
                eval_status, reason = "not_evaluable", "missing_program_score_or_n_below_10"
            rows.append(
                {
                    "cohort_label": cohort,
                    "cancer_code": cancer,
                    "comparison_layer": layer,
                    "program": program,
                    "program_label": LABELS[program],
                    "matched_sample_n": n,
                    "status": eval_status,
                    "reason": reason,
                    "spearman_rho": rho,
                    "ci_lower": ci_low,
                    "ci_upper": ci_high,
                    "p_value": p,
                    "fdr_family": f"apollo_{cohort}_{layer}",
                }
            )
    registry = pd.DataFrame(rows)
    registry["q_value_bh"] = np.nan
    for _, idx in registry.groupby("fdr_family").groups.items():
        registry.loc[list(idx), "q_value_bh"] = bh_fdr(registry.loc[list(idx), "p_value"])
    registry["evidence_status"] = "not_evaluable"
    evaluable = registry["status"].eq("evaluable")
    registry.loc[evaluable & (registry["q_value_bh"] < 0.05) & (registry["spearman_rho"] > 0), "evidence_status"] = "supported"
    registry.loc[evaluable & (registry["q_value_bh"] < 0.05) & (registry["spearman_rho"] < 0), "evidence_status"] = "opposite"
    registry.loc[evaluable & ~(registry["q_value_bh"] < 0.05), "evidence_status"] = "null_result"
    registry["claim_boundary"] = "Same-cancer public proteomic support by provenance; no harmonized donor crosswalk to CPTAC."
    out = source_dir / "apollo_concordance_registry_v2.csv"
    registry.to_csv(out, index=False)
    summary = {
        "registered_rows": len(registry),
        "evaluable_rows": int(evaluable.sum()),
        "supported_rows": int(registry["evidence_status"].eq("supported").sum()),
        "opposite_rows": int(registry["evidence_status"].eq("opposite").sum()),
        "null_result_rows": int(registry["evidence_status"].eq("null_result").sum()),
        "not_evaluable_rows": int(registry["evidence_status"].eq("not_evaluable").sum()),
        "duplicate_inference_key_rows": duplicate_n,
    }
    return registry, summary


def load_scrna_scores(root: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    paths = [
        root / "results" / "scrna" / "scrna_signature_scores_v0_7.csv",
        root / "results" / "threeca" / "threeca_signature_scores_v0_8.csv",
    ]
    for path in paths:
        if path.exists() and path.stat().st_size > 0:
            frame = pd.read_csv(path)
            if not frame.empty:
                frame["source_file"] = str(path)
                frames.append(frame)
    if not frames:
        return pd.DataFrame()
    scores = pd.concat(frames, ignore_index=True, sort=False)
    if "signature" not in scores.columns:
        label_to_program = {v: k for k, v in LABELS.items()}
        scores["signature"] = scores["signature_label"].map(label_to_program)
    key = ["source_id", "cancer_code", "donor_id", "signature", "compartment"]
    duplicate_n = int(scores.duplicated(key, keep=False).sum())
    if duplicate_n:
        raise RuntimeError(f"Single-cell donor-compartment score keys are duplicated: {duplicate_n}")
    return scores


def paired_delta(group: pd.DataFrame, left: str, right: str) -> pd.DataFrame:
    wide = group.pivot(index="donor_unit", columns="compartment", values="signature_score")
    if left not in wide or right not in wide:
        return pd.DataFrame(columns=["donor_unit", "delta"])
    delta = pd.to_numeric(wide[left] - wide[right], errors="coerce").dropna()
    return delta.rename("delta").reset_index()


def donor_provenance_audit(scores: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Audit local donor namespaces without assuming cross-repository identity."""
    source_donors = scores[["source_id", "dataset_title", "cancer_code", "donor_id"]].drop_duplicates()
    raw_id_source_n = source_donors.groupby("donor_id")["source_id"].nunique()
    rows: list[dict[str, Any]] = []
    for source_id, group in source_donors.groupby("source_id", sort=True):
        raw_ids = set(group["donor_id"].astype(str))
        exact_overlap = sum(int(raw_id_source_n.get(donor, 0) > 1) for donor in raw_ids)
        rows.append(
            {
                "source_id": source_id,
                "dataset_title": "; ".join(sorted(set(group["dataset_title"].astype(str)))),
                "publication_or_collection_identity": SOURCE_PUBLICATION_IDENTITY.get(source_id, "not resolved from retained local metadata"),
                "cancer_codes": ";".join(sorted(set(group["cancer_code"].astype(str)))),
                "donor_identifier_field": "source-provided donor_id or participant/sample field documented in source schema audit",
                "donor_field_integrity_note": SOURCE_DONOR_FIELD_NOTE.get(source_id, "source schema audit retained; no cross-source identity inference"),
                "donor_identifier_namespace": f"source:{source_id}",
                "observed_unique_source_specific_donor_ids": len(raw_ids),
                "exact_raw_donor_labels_reused_in_other_sources": exact_overlap,
                "exact_label_overlap_result": "none_observed" if exact_overlap == 0 else "observed",
                "cross_source_biological_identity_crosswalk": "unavailable",
                "global_biological_donor_independence": "not_established",
                "pooled_donor_level_inference_eligible": False,
                "adjudication": "retain source-specific donor-paired analysis; global pooled donor inference removed",
            }
        )
    audit = pd.DataFrame(rows)
    summary = {
        "audited_source_n": int(audit["source_id"].nunique()),
        "source_specific_donor_unit_n": int(scores["donor_unit"].nunique()),
        "exact_raw_donor_labels_reused_across_sources_n": int((raw_id_source_n > 1).sum()),
        "globally_unique_biological_donor_crosswalk_available": False,
        "global_pooled_donor_inference_retained": False,
        "adjudication": "No exact donor-label reuse was observed, but labels are source-specific and cannot establish biological independence across repositories or aggregated collections.",
    }
    return audit, summary


def signed_rank(delta: pd.Series, min_n: int = 6) -> tuple[str, str, float, float, float]:
    values = pd.to_numeric(delta, errors="coerce").dropna().to_numpy(dtype=float)
    n = len(values)
    if n < min_n:
        return "not_evaluable", f"paired_donor_n_below_{min_n}", np.nan, np.nan, np.nan
    if np.allclose(values, 0):
        return "evaluable", "all_paired_differences_zero", 0.0, 1.0, 0.5
    result = wilcoxon(values, alternative="two-sided", zero_method="wilcox", method="auto")
    return "evaluable", "", float(np.median(values)), float(result.pvalue), float(np.mean(values > 0))


def scrna_registries(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    scores = load_scrna_scores(root)
    if scores.empty:
        empty = pd.DataFrame()
        return empty, empty, {"status": "not_run_or_no_scores"}
    scores = scores.copy()
    scores["donor_unit"] = scores["source_id"].astype(str) + "::" + scores["cancer_code"].astype(str) + "::" + scores["donor_id"].astype(str)
    provenance, provenance_summary = donor_provenance_audit(scores)
    within_rows: list[dict[str, Any]] = []
    source_cancers = scores[["source_id", "cancer_code"]].drop_duplicates().itertuples(index=False)
    for source_id, cancer_code in source_cancers:
        subset = scores[scores["source_id"].eq(source_id) & scores["cancer_code"].eq(cancer_code)]
        for program in PROGRAMS:
            group = subset[subset["signature"].eq(program)]
            for left, right in CONTRASTS:
                delta = paired_delta(group, left, right)
                status, reason, median, p, positive_fraction = signed_rank(delta["delta"])
                within_rows.append(
                    {
                        "source_id": source_id,
                        "cancer_code": cancer_code,
                        "program": program,
                        "program_label": LABELS[program],
                        "contrast": f"{left}_minus_{right}",
                        "paired_donor_n": len(delta),
                        "median_delta": median,
                        "positive_fraction": positive_fraction,
                        "wilcoxon_p_value": p,
                        "status": status,
                        "reason": reason,
                        "fdr_family": f"single_cell_within_source_{source_id}_{cancer_code}_8_programs_x_3_contrasts",
                    }
                )
    within = pd.DataFrame(within_rows)
    within["q_value_bh"] = np.nan
    for _, idx in within.groupby("fdr_family").groups.items():
        positions = list(idx)
        within.loc[positions, "q_value_bh"] = bh_fdr(within.loc[positions, "wilcoxon_p_value"])
    within["within_source_support_window"] = (
        within["status"].eq("evaluable") & (within["paired_donor_n"] >= 6) & (within["q_value_bh"] < 0.05)
    )
    within["n_ge_10_sensitivity_window"] = (
        within["status"].eq("evaluable") & (within["paired_donor_n"] >= 10) & (within["q_value_bh"] < 0.05)
    )

    global_rows: list[dict[str, Any]] = []
    for program in PROGRAMS:
        group = scores[scores["signature"].eq(program)]
        for left, right in CONTRASTS:
            delta = paired_delta(group, left, right)
            values = pd.to_numeric(delta["delta"], errors="coerce").dropna()
            median = float(values.median()) if len(values) else np.nan
            positive_fraction = float((values > 0).mean()) if len(values) else np.nan
            global_rows.append(
                {
                    "program": program,
                    "program_label": LABELS[program],
                    "contrast": f"{left}_minus_{right}",
                    "paired_donor_n": len(delta),
                    "contributing_source_n": int(delta["donor_unit"].str.split("::").str[0].nunique()) if not delta.empty else 0,
                    "median_delta": median,
                    "positive_fraction": positive_fraction,
                    "wilcoxon_p_value": np.nan,
                    "status": "descriptive_only" if len(delta) else "not_evaluable",
                    "reason": "cross_source_biological_donor_independence_not_established" if len(delta) else "no_paired_source_specific_donor_units",
                    "fdr_family": "not_applicable_global_pooling_removed_after_donor_provenance_audit",
                }
            )
    global_registry = pd.DataFrame(global_rows)
    global_registry["q_value_bh"] = np.nan
    global_registry["main_claim_window"] = False
    global_registry["n_ge_10_sensitivity_window"] = False
    within["claim_boundary"] = "Donor-paired source mapping only; not proof of malignant-cell origin or mechanism."
    global_registry["claim_boundary"] = "Descriptive pooled source-specific donor-unit summary only; no P/q value or global inferential claim."
    out_dir = root / "results" / "scrna_combined"
    out_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(out_dir / "public_scrna_and_3ca_donor_compartment_scores_v2.csv", index=False)
    within.to_csv(out_dir / "public_scrna_within_source_wilcoxon_registry_v2.csv", index=False)
    global_registry.to_csv(out_dir / "public_scrna_global_wilcoxon_registry_v2.csv", index=False)
    provenance.to_csv(out_dir / "public_scrna_source_donor_provenance_audit_v2.csv", index=False)
    summary = {
        "score_rows": len(scores),
        "unique_donor_units": int(scores["donor_unit"].nunique()),
        "unique_source_ids": int(scores["source_id"].nunique()),
        "within_source_registered_rows": len(within),
        "within_source_evaluable_rows": int(within["status"].eq("evaluable").sum()),
        "within_source_support_rows": int(within["within_source_support_window"].sum()),
        "global_registered_rows": len(global_registry),
        "global_evaluable_rows": 0,
        "global_descriptive_rows": int(global_registry["status"].eq("descriptive_only").sum()),
        "global_main_claim_rows": int(global_registry["main_claim_window"].sum()),
        "global_n_ge_10_rows": int(global_registry["n_ge_10_sensitivity_window"].sum()),
        "donor_provenance_audit": provenance_summary,
    }
    return within, global_registry, summary


def spatial_summary(root: Path) -> dict[str, Any]:
    source_dir = root / "results" / "spatial"
    path = source_dir / "spatial_brca_sample_signature_summary_v0_7.csv"
    if not path.exists():
        return {"status": "not_run"}
    table = pd.read_csv(path)
    sample_col = "sample_id"
    return {
        "status": "executed",
        "summary_rows": len(table),
        "unique_samples": int(table[sample_col].nunique()) if sample_col in table else 0,
        "paired_sample_n_values": sorted(pd.to_numeric(table.get("paired_sample_n", pd.Series(dtype=float)), errors="coerce").dropna().astype(int).unique().tolist()),
        "inference_boundary": "Sample-level medians only; spots are not independent patients and no spot-level P value is used.",
    }


def main() -> int:
    args = parse_args()
    root = Path(args.analysis_root).resolve()
    apollo, apollo_summary = apollo_registry(root)
    within, global_registry, scrna_summary = scrna_registries(root)
    spatial = spatial_summary(root)
    tcga_summary = json.loads((root / "results" / "tcga" / "tcga_source_execution_summary.json").read_text())
    tcga_scores = pd.read_csv(root / "results" / "tcga" / "tcga_xena_signature_scores_v0_5.csv")
    tcga_duplicate = int(tcga_scores.duplicated(["patient_id", "signature"], keep=False).sum())
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tcga": {
            "primary_tumor_patients": int(tcga_summary["tcga_primary_tumor_patients"]),
            "scored_samples": int(tcga_summary["tcga_scored_samples"]),
            "patient_program_duplicate_rows": tcga_duplicate,
            "boundary": "RNA-only orthogonal sensitivity; not matched cross-omic validation.",
        },
        "apollo": apollo_summary,
        "single_cell": scrna_summary,
        "spatial": spatial,
    }
    out = root / "results" / "support_layers_summary_v2.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    critical = tcga_duplicate != 0 or apollo_summary["duplicate_inference_key_rows"] != 0
    return 2 if critical else 0


if __name__ == "__main__":
    raise SystemExit(main())
