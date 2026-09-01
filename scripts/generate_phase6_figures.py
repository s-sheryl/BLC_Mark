from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


COHORTS = (
    "TCGA-BRCA",
    "TCGA-LUAD",
    "TCGA-COAD",
)

TOP_N = 15


def load_rankings(repo_root: Path, cohort: str) -> pd.DataFrame:
    path = (
        repo_root
        / "results"
        / "phase5"
        / cohort
        / "prioritization_results.csv"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing Phase 5 ranking file: {path}"
        )

    df = pd.read_csv(path)

    required_columns = {
        "rank",
        "gene_id",
        "final_score",
    }

    missing = required_columns.difference(df.columns)

    if missing:
        raise ValueError(
            f"{cohort} is missing required columns: "
            f"{sorted(missing)}"
        )

    ranked = df.loc[
        df["rank"].notna()
        & df["final_score"].notna()
    ].copy()

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

    return ranked.head(TOP_N)


def create_top_biomarker_plot(
    dataframe: pd.DataFrame,
    cohort: str,
    output_path: Path,
) -> None:

    plot_data = dataframe.sort_values(
        "final_score",
        ascending=True,
    )

    plt.figure(figsize=(9, 6))

    plt.barh(
        plot_data["gene_id"],
        plot_data["final_score"],
    )

    plt.xlabel("BLC Mark Final Prioritization Score")
    plt.ylabel("Gene")

    plt.title(
        f"{cohort}: Top {TOP_N} Prioritized Biomarkers"
    )

    plt.xlim(0, 1)

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def create_cross_cohort_summary(
    repo_root: Path,
    output_directory: Path,
) -> None:

    records = []

    for cohort in COHORTS:
        rankings = load_rankings(
            repo_root,
            cohort,
        ).copy()

        rankings["cancer_cohort"] = cohort

        records.append(
            rankings[
                [
                    "cancer_cohort",
                    "gene_id",
                    "rank",
                    "final_score",
                ]
            ]
        )

    combined = pd.concat(
        records,
        ignore_index=True,
    )

    combined = combined.sort_values(
        by=["cancer_cohort", "rank"],
        kind="stable",
    )

    plt.figure(figsize=(10, 7))

    for cohort in COHORTS:
        subset = combined.loc[
            combined["cancer_cohort"] == cohort
        ]

        plt.plot(
            subset["rank"],
            subset["final_score"],
            marker="o",
            label=cohort,
        )

    plt.xlabel("Biomarker Rank")
    plt.ylabel("BLC Mark Final Prioritization Score")

    plt.title(
        f"Top {TOP_N} Biomarker Score Profiles Across Cohorts"
    )

    plt.xlim(1, TOP_N)
    plt.ylim(0, 1)

    plt.xticks(
        range(1, TOP_N + 1)
    )

    plt.legend()

    plt.tight_layout()

    output_path = (
        output_directory
        / f"cross_cohort_top_{TOP_N}_score_profiles.png"
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(
        "Cross-cohort summary: wrote figure to "
        f"{output_path}"
    )


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    output_directory = (
        repo_root
        / "results"
        / "phase6"
        / "figures"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    for cohort in COHORTS:

        rankings = load_rankings(
            repo_root,
            cohort,
        )

        cohort_slug = cohort.lower().replace(
            "tcga-",
            "",
        )

        output_path = (
            output_directory
            / f"{cohort_slug}_top_{TOP_N}_biomarkers.png"
        )

        create_top_biomarker_plot(
            rankings,
            cohort,
            output_path,
        )

        print(
            f"{cohort}: wrote figure to "
            f"{output_path}"
        )

    create_cross_cohort_summary(
        repo_root,
        output_directory,
    )


if __name__ == "__main__":
    main()