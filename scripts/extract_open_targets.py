"""Extract BLC Mark cancer associations from Open Targets 26.06."""

import re
from pathlib import Path

import duckdb


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INDEX_FILE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "open_targets"
    / "association_index.html"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "open_targets"
    / "blc_mark_cancer_associations_26.06.parquet"
)

BASE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/"
    "opentargets/platform/26.06/output/"
    "association_overall_direct/"
)

DISEASE_IDS = (
    "MONDO_0007254",  # breast cancer
    "MONDO_0005061",  # lung adenocarcinoma
    "MONDO_0002271",  # colon adenocarcinoma
)


def main() -> None:
    html = INDEX_FILE.read_text(
        encoding="utf-8"
    )

    filenames = sorted(
        set(
            re.findall(
                r'href="([^"]+\.parquet)"',
                html,
            )
        )
    )

    if len(filenames) != 14:
        raise RuntimeError(
            "Expected 14 Open Targets association shards, "
            f"found {len(filenames)}."
        )

    urls = [
        BASE_URL + filename
        for filename in filenames
    ]

    url_sql = ",\n".join(
        f"'{url}'"
        for url in urls
    )

    disease_sql = ", ".join(
        f"'{disease_id}'"
        for disease_id in DISEASE_IDS
    )

    con = duckdb.connect()

    con.execute(
        "INSTALL httpfs;"
    )

    con.execute(
        "LOAD httpfs;"
    )

    con.execute(
        "SET http_retries = 10;"
    )

    con.execute(
        "SET http_retry_wait_ms = 5000;"
    )

    query = f"""
    COPY (
        SELECT
            diseaseId,
            targetId,
            associationScore,
            evidenceCount
        FROM read_parquet([
            {url_sql}
        ])
        WHERE diseaseId IN ({disease_sql})
        ORDER BY
            diseaseId,
            targetId
    )
    TO '{OUTPUT_FILE.as_posix()}'
    (
        FORMAT PARQUET,
        COMPRESSION ZSTD
    );
    """

    print(
        "Extracting only BLC Mark cancer associations..."
    )

    con.execute(query)

    result = con.execute(
        f"""
        SELECT
            diseaseId,
            COUNT(*) AS rows,
            COUNT(DISTINCT targetId) AS targets
        FROM read_parquet(
            '{OUTPUT_FILE.as_posix()}'
        )
        GROUP BY diseaseId
        ORDER BY diseaseId
        """
    ).fetchall()

    print()
    print("Extraction complete:")
    print(OUTPUT_FILE)

    for disease_id, rows, targets in result:
        print(
            disease_id,
            "| rows:",
            rows,
            "| targets:",
            targets,
        )


if __name__ == "__main__":
    main()