#!/usr/bin/env python3
from __future__ import annotations
from sci27_portability import expand_legacy_paths

import argparse
import json
import platform
import re
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path.cwd()
VERSION = "v0_7"
BOOTSTRAP_N = 1000
RANDOM_SEED = 20260730
MIN_SIGNATURE_GENES = 3

SIGNATURE_LABELS = {
    "antigen_presentation_mhc_i": "MHC-I",
    "antigen_presentation_mhc_ii": "MHC-II",
    "checkpoint_exhaustion_context": "Checkpoint",
    "cytolytic_t_cell_context": "Cytolytic",
    "ifn_gamma_response": "IFN-gamma",
    "myeloid_inflammatory_context": "Myeloid",
    "proteasome_antigen_processing": "Proteasome",
    "tgfb_emt_exclusion_context": "TGF/EMT",
}

DISCOVERY_COMPARISON = {
    "proteomics": "rna_vs_proteomics",
    "phosphoproteomics": "rna_vs_phosphoproteomics",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="BIO-05 APOLLO independent cohort signature concordance audit."
    )
    parser.add_argument(
        "--project-root",
        default=str(PROJECT_ROOT),
        help="Writable output root. Relative output/report paths are resolved below this directory.",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the five APOLLO source files listed in the public runbook.",
    )
    parser.add_argument(
        "--signature-table",
        default="04_processed/discovery_v0_1/signature_gene_sets_v0_1.csv",
    )
    parser.add_argument(
        "--discovery-corr",
        default="04_processed/discovery_v0_4/bcm_matched_correlation_summary_v0_4.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="04_processed/discovery_v0_7_public_validation",
    )
    parser.add_argument(
        "--report",
        default="05_reports/public_apollo_external_cohort_v0_7_report.md",
    )
    parser.add_argument("--bootstrap-n", type=int, default=BOOTSTRAP_N)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--min-signature-genes", type=int, default=MIN_SIGNATURE_GENES)
    return parser.parse_args()


def safe_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value)


def load_signature_genes(path: Path) -> dict[str, list[str]]:
    table = pd.read_csv(path)
    return {
        signature: sorted({safe_str(g).strip().upper() for g in group["gene_symbol"] if safe_str(g).strip()})
        for signature, group in table.groupby("signature", dropna=False)
    }


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for col in out.columns:
        values = pd.to_numeric(out[col], errors="coerce")
        std = float(values.std(ddof=0))
        if not np.isfinite(std) or std == 0:
            out[col] = np.nan
        else:
            out[col] = (values - float(values.mean())) / std
    return out


def collapse_gene_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    clean = frame.copy()
    clean.index = clean.index.astype(str).str.upper()
    clean = clean[clean.index.astype(str) != ""]
    clean = clean.groupby(level=0).mean()
    return clean


def score_signature_matrix(
    expr: pd.DataFrame,
    signature_genes: dict[str, list[str]],
    min_signature_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_z = zscore_columns(expr.T).T
    score_rows: list[dict[str, Any]] = []
    recovery_rows: list[dict[str, Any]] = []
    for signature, genes in signature_genes.items():
        recovered = [g for g in genes if g in expr_z.index]
        recovery_rows.append(
            {
                "signature": signature,
                "signature_label": SIGNATURE_LABELS.get(signature, signature),
                "requested_gene_count": len(genes),
                "recovered_gene_count": len(recovered),
            }
        )
        if len(recovered) < min_signature_genes:
            continue
        values = expr_z.loc[recovered].mean(axis=0)
        for sample_id, score in values.items():
            score_rows.append(
                {
                    "signature": signature,
                    "signature_label": SIGNATURE_LABELS.get(signature, signature),
                    "sample_id": sample_id,
                    "signature_score": float(score),
                    "n_genes_used": len(recovered),
                }
            )
    return pd.DataFrame(score_rows), pd.DataFrame(recovery_rows)


def bootstrap_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    bootstrap_n: int,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    obs = float(spearmanr(x, y, nan_policy="omit").statistic)
    n = len(x)
    vals = []
    for _ in range(bootstrap_n):
        idx = rng.integers(0, n, size=n)
        rho = spearmanr(x[idx], y[idx], nan_policy="omit").statistic
        if np.isfinite(rho):
            vals.append(float(rho))
    if not vals:
        return obs, np.nan, np.nan
    arr = np.asarray(vals, dtype=float)
    return obs, float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))


def standardize_luad_proteome_sample(col: str) -> str | None:
    match = re.match(r"AP\.([A-Z0-9]+)\.P\d+$", safe_str(col))
    return f"AP-{match.group(1)}" if match else None


def standardize_luad_phospho_sample(col: str) -> tuple[str | None, str | None]:
    sample_match = re.match(r"AP\.([A-Z0-9]+)\.\.(phospho|relative)\.$", safe_str(col))
    if not sample_match:
        return None, None
    return f"AP-{sample_match.group(1)}", sample_match.group(2)


def load_apollo_luad_rna(path: Path) -> pd.DataFrame:
    matrices = []
    with tarfile.open(path, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".counts.quant.txt")]
        for member in members:
            sample_id = Path(member.name).name.replace(".counts.quant.txt", "")
            handle = tar.extractfile(member)
            if handle is None:
                continue
            df = pd.read_csv(handle, sep="\t")
            df = df[["Gene", "count-UQ"]].copy()
            df.columns = ["gene", sample_id]
            matrices.append(df)
    merged = matrices[0]
    for df in matrices[1:]:
        merged = merged.merge(df, on="gene", how="outer")
    merged = merged.fillna(0.0)
    merged = merged.set_index("gene")
    merged = np.log2(merged.astype(float) + 1.0)
    return collapse_gene_matrix(merged)


def load_apollo_luad_proteome(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    sample_cols = [c for c in df.columns if standardize_luad_proteome_sample(c)]
    rename = {c: standardize_luad_proteome_sample(c) for c in sample_cols}
    mat = df[["Gene", *sample_cols]].copy()
    mat = mat.rename(columns={"Gene": "gene", **rename})
    mat = mat.set_index("gene")
    return collapse_gene_matrix(mat.astype(float))


def load_apollo_luad_phospho(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    keep_cols: list[str] = []
    rename: dict[str, str] = {}
    for col in df.columns:
        sample_id, layer = standardize_luad_phospho_sample(col)
        if sample_id and layer == "relative":
            keep_cols.append(col)
            rename[col] = sample_id
    mat = df[["Primary_Gene", *keep_cols]].copy()
    mat = mat.rename(columns={"Primary_Gene": "gene", **rename})
    mat = mat.set_index("gene")
    return collapse_gene_matrix(mat.astype(float))


def load_apollo_ov_rna(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "gene"}).set_index("gene")
    return collapse_gene_matrix(df.astype(float))


def load_apollo_ov_proteome(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    sample_cols = [c for c in df.columns if re.fullmatch(r"A\d+", safe_str(c))]
    mat = df[["Gene", *sample_cols]].copy().rename(columns={"Gene": "gene"}).set_index("gene")
    return collapse_gene_matrix(mat.astype(float))


def build_concordance(
    cohort_label: str,
    cancer_code: str,
    layer: str,
    rna_scores: pd.DataFrame,
    omic_scores: pd.DataFrame,
    bootstrap_n: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    merged = rna_scores.merge(
        omic_scores,
        on=["signature", "signature_label", "sample_id"],
        how="inner",
        suffixes=("_rna", "_omic"),
    )
    rows: list[dict[str, Any]] = []
    for (signature, signature_label), group in merged.groupby(["signature", "signature_label"], dropna=False):
        x = pd.to_numeric(group["signature_score_rna"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(group["signature_score_omic"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        if x.size < 10:
            continue
        rho, ci_low, ci_high = bootstrap_spearman_ci(x, y, bootstrap_n=bootstrap_n, rng=rng)
        rows.append(
            {
                "cohort_label": cohort_label,
                "cancer_code": cancer_code,
                "comparison_layer": layer,
                "signature": signature,
                "signature_label": signature_label,
                "n": int(x.size),
                "spearman_rho": rho,
                "ci_lower": ci_low,
                "ci_upper": ci_high,
            }
        )
    return pd.DataFrame(rows)


def add_discovery_consistency(
    concordance: pd.DataFrame,
    discovery_corr: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if concordance.empty:
        return concordance, pd.DataFrame()
    discovery = discovery_corr[
        discovery_corr["comparison"].isin(DISCOVERY_COMPARISON.values())
        & discovery_corr["comparison_scope"].eq("within_cancer")
    ].copy()
    discovery = discovery.rename(columns={"comparison_layer": "comparison_layer_discovery"})
    layer_map = {
        "proteomics": "proteomics",
        "phosphoproteomics": "phosphoproteomics",
    }
    merged = concordance.merge(
        discovery[
            [
                "cancer",
                "comparison_layer_discovery",
                "signature",
                "signature_label",
                "spearman_rho",
            ]
        ].rename(columns={"cancer": "cancer_code", "spearman_rho": "discovery_rho"}),
        left_on=["cancer_code", "comparison_layer", "signature", "signature_label"],
        right_on=["cancer_code", "comparison_layer_discovery", "signature", "signature_label"],
        how="left",
    )
    merged["sign_agreement_vs_discovery"] = np.where(
        np.isfinite(merged["spearman_rho"]) & np.isfinite(merged["discovery_rho"]),
        np.sign(merged["spearman_rho"]) == np.sign(merged["discovery_rho"]),
        np.nan,
    )

    summary_rows: list[dict[str, Any]] = []
    for (cohort_label, cancer_code, layer), group in merged.groupby(
        ["cohort_label", "cancer_code", "comparison_layer"], dropna=False
    ):
        valid = group[np.isfinite(group["spearman_rho"]) & np.isfinite(group["discovery_rho"])].copy()
        if valid.empty:
            continue
        rho_vec = float(
            spearmanr(valid["spearman_rho"], valid["discovery_rho"], nan_policy="omit").statistic
        )
        summary_rows.append(
            {
                "cohort_label": cohort_label,
                "cancer_code": cancer_code,
                "comparison_layer": layer,
                "signature_overlap_n": int(valid.shape[0]),
                "rank_correlation_vs_discovery_same_cancer": rho_vec,
                "positive_sign_agreement_fraction": float(valid["sign_agreement_vs_discovery"].mean()),
                "median_abs_rho_diff_vs_discovery": float(
                    (valid["spearman_rho"] - valid["discovery_rho"]).abs().median()
                ),
            }
        )
    return merged, pd.DataFrame(summary_rows)


def write_report(
    path: Path,
    concordance: pd.DataFrame,
    consistency: pd.DataFrame,
) -> None:
    lines = [
        "# BIO-05 APOLLO External Cohort v0.7",
        "",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        f"- Platform: {platform.platform()}",
        "",
        "## Cohorts",
        "",
        "- `APOLLO-LUAD-2022`: matched RNA, global proteome, and gene-collapsed phosphosite relative matrix across 87 public samples.",
        "- `APOLLO-OV-2024`: matched whole-tumor RNA and whole-tumor imputed global proteome across public ovarian samples.",
        "",
        "## Within-Cancer Concordance",
        "",
    ]
    if concordance.empty:
        lines.append("- No concordance result was generated.")
    else:
        for row in concordance.sort_values(["cohort_label", "comparison_layer", "signature_label"]).itertuples():
            lines.append(
                f"- `{row.cohort_label}` `{row.cancer_code}` `{row.comparison_layer}` `{row.signature_label}`: n={row.n}, rho={row.spearman_rho:.3f}, 95% CI [{row.ci_lower:.3f}, {row.ci_upper:.3f}]"
            )
    lines.extend(["", "## Discovery Consistency", ""])
    if consistency.empty:
        lines.append("- No same-cancer discovery consistency summary was available.")
    else:
        for row in consistency.sort_values(["cohort_label", "comparison_layer"]).itertuples():
            lines.append(
                f"- `{row.cohort_label}` `{row.cancer_code}` `{row.comparison_layer}`: overlap={row.signature_overlap_n}, rank_rho_vs_discovery={row.rank_correlation_vs_discovery_same_cancer:.3f}, sign_agreement={row.positive_sign_agreement_fraction:.3f}, median_abs_diff={row.median_abs_rho_diff_vs_discovery:.3f}"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).expanduser().resolve()
    input_dir = Path(args.input_dir).expanduser().resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = root / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)

    signature_genes = load_signature_genes(root / args.signature_table)
    rng = np.random.default_rng(args.seed)

    required_inputs = {
        "luad_rna": input_dir / "APOLLO_LUAD_RNA_quant_files_20181226.tar.gz",
        "luad_prot": input_dir / "APOLLO1_Level3_Gprot_v061418_n87.csv",
        "luad_phos": input_dir / "APOLLO1_Phospho_Level3_n87.csv",
        "ov_rna": input_dir / "APOLLO_OV_WholeTumor_RNASeq_NormalizedCounts.csv",
        "ov_prot": input_dir / "APOLLO_OV_WholeTumor_GlobalProteomics_imputed.csv",
    }
    missing = [str(path) for path in required_inputs.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing APOLLO input file(s): " + ", ".join(missing))

    luad_rna = load_apollo_luad_rna(required_inputs["luad_rna"])
    luad_prot = load_apollo_luad_proteome(required_inputs["luad_prot"])
    luad_phos = load_apollo_luad_phospho(required_inputs["luad_phos"])
    ov_rna = load_apollo_ov_rna(required_inputs["ov_rna"])
    ov_prot = load_apollo_ov_proteome(required_inputs["ov_prot"])

    matrices = {
        ("APOLLO-LUAD", "luad", "rna"): luad_rna,
        ("APOLLO-LUAD", "luad", "proteomics"): luad_prot,
        ("APOLLO-LUAD", "luad", "phosphoproteomics"): luad_phos,
        ("APOLLO-OV", "ov", "rna"): ov_rna,
        ("APOLLO-OV", "ov", "proteomics"): ov_prot,
    }

    score_frames: list[pd.DataFrame] = []
    recovery_frames: list[pd.DataFrame] = []
    for (cohort_label, cancer_code, layer), matrix in matrices.items():
        score_df, recovery_df = score_signature_matrix(
            matrix,
            signature_genes=signature_genes,
            min_signature_genes=args.min_signature_genes,
        )
        score_df["cohort_label"] = cohort_label
        score_df["cancer_code"] = cancer_code
        score_df["comparison_layer"] = layer
        recovery_df["cohort_label"] = cohort_label
        recovery_df["cancer_code"] = cancer_code
        recovery_df["comparison_layer"] = layer
        score_frames.append(score_df)
        recovery_frames.append(recovery_df)

    scores = pd.concat(score_frames, ignore_index=True)
    recovery = pd.concat(recovery_frames, ignore_index=True)
    scores.to_csv(out_dir / f"apollo_signature_scores_{VERSION}.csv", index=False)
    recovery.to_csv(out_dir / f"apollo_signature_recovery_{VERSION}.csv", index=False)

    concordance_frames: list[pd.DataFrame] = []
    for cohort_label, cancer_code in [("APOLLO-LUAD", "luad"), ("APOLLO-OV", "ov")]:
        rna_scores = scores[
            (scores["cohort_label"] == cohort_label)
            & (scores["cancer_code"] == cancer_code)
            & (scores["comparison_layer"] == "rna")
        ].copy()
        for layer in ["proteomics", "phosphoproteomics"]:
            omic_scores = scores[
                (scores["cohort_label"] == cohort_label)
                & (scores["cancer_code"] == cancer_code)
                & (scores["comparison_layer"] == layer)
            ].copy()
            if omic_scores.empty:
                continue
            concordance_frames.append(
                build_concordance(
                    cohort_label=cohort_label,
                    cancer_code=cancer_code,
                    layer=layer,
                    rna_scores=rna_scores,
                    omic_scores=omic_scores,
                    bootstrap_n=args.bootstrap_n,
                    rng=rng,
                )
            )

    concordance = pd.concat(concordance_frames, ignore_index=True) if concordance_frames else pd.DataFrame()
    discovery_corr = pd.read_csv(root / args.discovery_corr)
    concordance_annotated, consistency = add_discovery_consistency(concordance, discovery_corr)
    concordance_annotated.to_csv(out_dir / f"apollo_concordance_{VERSION}.csv", index=False)
    consistency.to_csv(out_dir / f"apollo_consistency_vs_discovery_{VERSION}.csv", index=False)

    write_report(report_path, concordance_annotated, consistency)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
