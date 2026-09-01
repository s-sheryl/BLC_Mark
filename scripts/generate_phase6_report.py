from pathlib import Path
from datetime import datetime, timezone

import pandas as pd


COHORTS = {
    "TCGA-BRCA": {
        "name": "Breast Cancer",
        "prefix": "brca",
    },
    "TCGA-LUAD": {
        "name": "Lung Adenocarcinoma",
        "prefix": "luad",
    },
    "TCGA-COAD": {
        "name": "Colorectal Cancer",
        "prefix": "coad",
    },
}


def load_phase3_results(
    repo_root: Path,
    cohort: str,
    prefix: str,
) -> pd.DataFrame:
    path = (
        repo_root
        / "results"
        / "phase3"
        / cohort
        / f"{prefix}_v1_de_differential_expression_results.csv"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing Phase 3 results: {path}"
        )

    return pd.read_csv(path)


def load_phase5_results(
    repo_root: Path,
    cohort: str,
) -> pd.DataFrame:
    path = (
        repo_root
        / "results"
        / "phase5"
        / cohort
        / "prioritization_results.csv"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing Phase 5 results: {path}"
        )

    return pd.read_csv(path)


def load_top_table(
    repo_root: Path,
    prefix: str,
) -> pd.DataFrame:
    path = (
        repo_root
        / "results"
        / "phase6"
        / "tables"
        / f"{prefix}_top_25_biomarkers.csv"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing Phase 6 top table: {path}"
        )

    return pd.read_csv(path)


def format_scientific(value: float) -> str:
    return f"{value:.3e}"


def format_score(value: float) -> str:
    return f"{value:.4f}"


def build_cohort_section(
    repo_root: Path,
    cohort: str,
    cancer_name: str,
    prefix: str,
) -> str:

    phase3 = load_phase3_results(
        repo_root,
        cohort,
        prefix,
    )

    phase5 = load_phase5_results(
        repo_root,
        cohort,
    )

    top = load_top_table(
        repo_root,
        prefix,
    )

    significant_count = int(
        phase3["significant"]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    scored_count = int(
        phase5["final_score"].notna().sum()
    )

    unavailable_count = int(
        phase5["final_score"].isna().sum()
    )

    top_gene = top.iloc[0]

    lines = [
        f"### {cohort} — {cancer_name}",
        "",
        f"- Significant differential-expression candidates: "
        f"**{significant_count:,}**",
        f"- Candidates retained for prioritization: "
        f"**{len(phase5):,}**",
        f"- Candidates with complete final scores: "
        f"**{scored_count:,}**",
        f"- Candidates retained with unavailable final scores: "
        f"**{unavailable_count:,}**",
        f"- Highest-ranked candidate: "
        f"**{top_gene['gene_id']}**",
        f"- Highest final prioritization score: "
        f"**{format_score(float(top_gene['final_score']))}**",
        "",
        "| Rank | Gene | Effect size | Adjusted p-value | "
        "Clinical category | Cross-cancer cohorts | Pathways | Final score |",
        "|---:|---|---:|---:|---|---:|---:|---:|",
    ]

    for _, row in top.head(10).iterrows():
        lines.append(
            f"| {int(row['rank'])} "
            f"| {row['gene_id']} "
            f"| {float(row['effect_size']):.4f} "
            f"| {format_scientific(float(row['adjusted_p_value']))} "
            f"| {row['clinical_category']} "
            f"| {int(row['cross_cancer_cohort_count'])} "
            f"| {int(row['pathway_count'])} "
            f"| {format_score(float(row['final_score']))} |"
        )

    lines.extend(
        [
            "",
            f"![{cohort} top biomarkers]"
            f"(../figures/{prefix}_top_15_biomarkers.png)",
            "",
        ]
    )

    return "\n".join(lines)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    output_directory = (
        repo_root
        / "results"
        / "phase6"
        / "report"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_path = (
        output_directory
        / "BLC_Mark_V1_scientific_report.md"
    )

    generated_at = datetime.now(
        timezone.utc
    ).isoformat()

    sections = []

    for cohort, config in COHORTS.items():
        sections.append(
            build_cohort_section(
                repo_root=repo_root,
                cohort=cohort,
                cancer_name=config["name"],
                prefix=config["prefix"],
            )
        )

    report = f"""# BLC Mark V1 — Scientific Results Report

## Overview

BLC Mark V1 is a reproducible cancer biomarker discovery,
evidence-integration, and prioritization workflow applied to three
TCGA cancer cohorts:

- TCGA-BRCA — Breast Cancer
- TCGA-LUAD — Lung Adenocarcinoma
- TCGA-COAD — Colorectal Cancer

Version 1 uses RNA-seq gene-expression data only.

The implemented workflow consists of six phases:

1. Scientific & Engineering Foundation
2. Data Acquisition & Preparation
3. Biomarker Discovery
4. Evidence Integration
5. Biomarker Prioritization
6. Scientific Outputs & Reproducibility

Generated at: `{generated_at}`

---

## Scientific Question

Can publicly available transcriptomic cancer datasets be integrated
into a transparent and reproducible evidence framework to systematically
prioritize candidate cancer biomarkers?

---

## Analysis Design

Differential-expression analysis compared Primary Tumor samples with
Solid Tissue Normal samples.

The implemented Version 1 differential-expression workflow used:

- normalized log2 expression data
- Welch's unequal-variance two-sample t-test
- Benjamini-Hochberg multiple-testing correction
- adjusted p-value significance threshold of 0.05
- no additional gene-expression filtering

Significant genes were passed to the evidence-integration phase.

Evidence categories used in prioritization included:

- differential-expression evidence
- cancer-association evidence
- clinical/prognostic evidence
- cross-cancer evidence
- functional and pathway context

The Version 1 prioritization model used four equally weighted scoring
components:

- differential-expression score: 0.25
- cancer-association score: 0.25
- clinical score: 0.25
- cross-cancer score: 0.25

A high final BLC Mark score indicates stronger prioritization under the
defined evidence framework. It does **not** constitute independent
clinical validation or proof of biomarker utility.

---

## Cohort Results

{chr(10).join(sections)}

---

## Cross-Cohort Prioritization

![Cross-cohort score profiles](../figures/cross_cohort_top_15_score_profiles.png)

The cross-cohort figure compares the final prioritization-score profiles
of the fifteen highest-ranked candidates in each cohort. Scores are
interpreted within the BLC Mark Version 1 prioritization framework and
should not be interpreted as clinical effect sizes or probabilities.

---

## Interpretation

The three cohorts produced distinct prioritized biomarker profiles.

TCGA-BRCA contained a mixture of validated, potential, favorable,
unfavorable, and unprognostic evidence categories among its highest-ranked
candidates.

TCGA-LUAD showed a particularly strong concentration of candidates
classified by the integrated evidence layer as validated prognostic and
unfavorable among the highest-ranked results.

TCGA-COAD produced a more heterogeneous top-ranked profile containing
both favorable and unfavorable prognostic evidence categories.

These observations describe the evidence classifications captured by
the BLC Mark workflow. They do not establish causal roles, diagnostic
performance, treatment response, or prospective clinical validity.

---

## Missing and Unavailable Evidence

BLC Mark distinguishes unavailable evidence from negative evidence.

Candidates lacking sufficient evidence to calculate every required
scoring component are retained rather than silently discarded.
Their final score and rank remain unavailable instead of being replaced
with zero.

This distinction preserves traceability between:

- absence of supporting evidence,
- evidence explicitly representing no support, and
- evidence that could not be resolved or retrieved.

---

## Reproducibility and Traceability

BLC Mark preserves machine-readable outputs for each analytical stage,
including:

- differential-expression results
- analysis metadata
- QC reports
- integrated evidence records
- prioritization results
- final summary tables
- scientific figures

Phase-specific metadata records analysis configuration, cohort identity,
software or scoring versions, and relevant upstream file paths.

The Phase 6 tables and figures are generated programmatically from
the frozen Phase 5 outputs rather than being manually reconstructed.

---

## Limitations

Version 1 has several important limitations.

First, the analysis is restricted to RNA-seq gene-expression data and
does not integrate genomic variants, methylation, proteomics, or other
omics layers.

Second, differential expression was evaluated using the implemented
Welch t-test workflow on normalized log2 expression data rather than a
raw-count negative-binomial model.

Third, external evidence availability is not uniform across genes.
Some candidate identifiers could not be resolved, and some candidates
therefore remain without complete final scores.

Fourth, prioritization reflects the predefined BLC Mark scoring system
and its Version 1 equal-weight configuration. Rankings should therefore
be interpreted as transparent computational prioritization, not as
clinical validation.

Finally, the analyzed Version 1 release contains three cohorts:
TCGA-BRCA, TCGA-LUAD, and TCGA-COAD.

---

## Conclusion

BLC Mark V1 demonstrates an end-to-end, reproducible workflow for
transcriptomic cancer biomarker discovery, external evidence integration,
and transparent candidate prioritization across three TCGA cancer
cohorts.

The system retains significant candidates through the evidence and
ranking stages, explicitly represents unavailable evidence, and produces
auditable machine-readable outputs alongside researcher-facing tables
and figures.

The resulting ranked candidates provide hypotheses for further
biological investigation and independent experimental or clinical
validation rather than definitive clinical biomarkers.
"""

    report_path.write_text(
        report,
        encoding="utf-8",
    )

    print(
        f"Wrote scientific report to {report_path}"
    )


if __name__ == "__main__":
    main()