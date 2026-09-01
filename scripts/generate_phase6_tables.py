from pathlib import Path

import pandas as pd


COHORTS = (
    "TCGA-BRCA",
    "TCGA-LUAD",
    "TCGA-COAD",
)

TOP_N = 25

INPUT_FILENAME = "prioritization_results.csv"

OUTPUT_COLUMNS = [
    "rank",
    "gene_id",
    "effect_size",
    "adjusted_p_value",
    "de_score",
    "cancer_association_score",
    "clinical_category",
    "clinical_score",
    "cross_cancer_cohort_count",
    "cross_cancer_score",
    "pathway_count",
    "final_score",
]


def load_prioritization_results(repo_root: Path, cohort: str) -> pd.DataFrame:
    input_path = (
        repo_root
        / "results"
        / "phase5"
        / cohort
        / INPUT_FILENAME
    )

    if not input_path.exists():
        raise FileNotFoundError(
            f"Phase 5 prioritization results not found: {input_path}"
        )

    dataframe = pd.read_csv(input_path)

    missing_columns = [
        column
        for column in OUTPUT_COLUMNS
        if column not in dataframe.columns
    ]

    if missing_columns:
        raise ValueError(
            f"{cohort} is missing required Phase 5 columns: "
            f"{missing_columns}"
        )

    return dataframe


def build_top_table(dataframe: pd.DataFrame, cohort: str) -> pd.DataFrame:
    ranked = dataframe.loc[
        dataframe["rank"].notna()
        & dataframe["final_score"].notna()
    ].copy()

    if ranked.empty:
        raise ValueError(
            f"{cohort} contains no scored and ranked candidates."
        )

    ranked["rank"] = pd.to_numeric(
        ranked["rank"],
        errors="raise",
    ).astype(int)

    ranked["final_score"] = pd.to_numeric(
        ranked["final_score"],
        errors="raise",
    )

    ranked = ranked.sort_values(
        by=["rank", "gene_id"],
        ascending=[True, True],
        kind="stable",
    )

    top_table = ranked.head(TOP_N)[OUTPUT_COLUMNS].copy()

    expected_ranks = list(
        range(1, len(top_table) + 1)
    )

    if top_table["rank"].tolist() != expected_ranks:
        raise ValueError(
            f"{cohort} top rankings are not sequential from rank 1."
        )

    top_table.insert(
        1,
        "cancer_cohort",
        cohort,
    )

    return top_table


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    output_directory = (
        repo_root
        / "results"
        / "phase6"
        / "tables"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    combined_tables = []

    for cohort in COHORTS:
        dataframe = load_prioritization_results(
            repo_root,
            cohort,
        )

        top_table = build_top_table(
            dataframe,
            cohort,
        )

        cohort_slug = cohort.lower().replace(
            "tcga-",
            "",
        )

        output_path = (
            output_directory
            / f"{cohort_slug}_top_{TOP_N}_biomarkers.csv"
        )

        top_table.to_csv(
            output_path,
            index=False,
        )

        combined_tables.append(top_table)

        print(
            f"{cohort}: wrote {len(top_table)} biomarkers "
            f"to {output_path}"
        )

    combined = pd.concat(
        combined_tables,
        ignore_index=True,
    )

    combined_output_path = (
        output_directory
        / f"combined_top_{TOP_N}_biomarkers.csv"
    )

    combined.to_csv(
        combined_output_path,
        index=False,
    )

    print(
        f"Combined table: wrote {len(combined)} rows "
        f"to {combined_output_path}"
    )


if __name__ == "__main__":
    main()