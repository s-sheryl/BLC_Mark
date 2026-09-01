"""Open Targets cancer-association evidence for BLC Mark Phase 4."""

import json
from urllib.request import Request, urlopen

from .models import EvidenceRecord, EvidenceType


OPEN_TARGETS_GRAPHQL_URL = (
    "https://api.platform.opentargets.org/api/v4/graphql"
)


BLC_MARK_CANCER_TERMS = {
    "TCGA-BRCA": (
        "breast cancer",
        "breast carcinoma",
        "breast invasive carcinoma",
    ),
    "TCGA-LUAD": (
        "lung adenocarcinoma",
    ),
    "TCGA-COAD": (
        "colon adenocarcinoma",
        "colorectal cancer",
        "colon cancer",
    ),
}


def _post_graphql(query: str) -> dict:
    """Send a GraphQL request to Open Targets and return decoded JSON."""
    if not isinstance(query, str) or not query.strip():
        raise ValueError("query must be a non-empty string.")

    payload = json.dumps(
        {"query": query}
    ).encode("utf-8")

    request = Request(
        OPEN_TARGETS_GRAPHQL_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "BLC-Mark/1.0",
        },
        method="POST",
    )

    with urlopen(request, timeout=30) as response:
        result = json.loads(
            response.read().decode("utf-8")
        )

    if result.get("errors"):
        raise RuntimeError(
            f"Open Targets GraphQL error: {result['errors']}"
        )

    return result


def fetch_target_disease_associations(
    ensembl_id: str,
    *,
    page_size: int = 100,
) -> dict | None:
    """
    Retrieve disease associations for one Open Targets target.

    Parameters
    ----------
    ensembl_id:
        Ensembl gene identifier used by Open Targets.

    page_size:
        Number of disease associations requested.

    Returns
    -------
    dict | None
        Target record, or None when Open Targets has no target.
    """
    if not isinstance(ensembl_id, str) or not ensembl_id.strip():
        raise ValueError(
            "ensembl_id must be a non-empty string."
        )

    if not isinstance(page_size, int) or page_size < 1:
        raise ValueError(
            "page_size must be a positive integer."
        )

    identifier = ensembl_id.strip()

    query = f"""
    query {{
      target(ensemblId: "{identifier}") {{
        id
        approvedSymbol
        associatedDiseases(
          page: {{index: 0, size: {page_size}}}
        ) {{
          count
          rows {{
            score
            disease {{
              id
              name
            }}
          }}
        }}
      }}
    }}
    """

    payload = _post_graphql(query)

    target = payload.get("data", {}).get("target")

    if not isinstance(target, dict):
        return None

    return target


def collect_cancer_association_evidence(
    gene_symbol: str,
    ensembl_id: str,
    cancer_cohort: str,
    *,
    retrieved_at: str,
    source_version: str,
    page_size: int = 100,
) -> list[EvidenceRecord]:
    """
    Convert Open Targets disease associations into cohort-relevant
    cancer evidence records.

    Only disease names matching the configured BLC Mark cancer cohort
    are retained.

    The Open Targets association score is preserved for provenance and
    interpretation only. No Phase 4 ranking or weighting is performed.
    """
    if not isinstance(gene_symbol, str) or not gene_symbol.strip():
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    if not isinstance(ensembl_id, str) or not ensembl_id.strip():
        raise ValueError(
            "ensembl_id must be a non-empty string."
        )

    if not isinstance(source_version, str) or not source_version.strip():
        raise ValueError(
            "source_version must be a non-empty string."
        )

    if cancer_cohort not in BLC_MARK_CANCER_TERMS:
        raise ValueError(
            f"Unsupported BLC Mark cancer cohort: {cancer_cohort}"
        )

    symbol = gene_symbol.strip()
    identifier = ensembl_id.strip()
    accepted_terms = BLC_MARK_CANCER_TERMS[cancer_cohort]

    target = fetch_target_disease_associations(
        identifier,
        page_size=page_size,
    )

    if target is None:
        return []

    approved_symbol = str(
        target.get("approvedSymbol", "")
    ).strip()

    if approved_symbol and approved_symbol != symbol:
        return []

    associated = target.get("associatedDiseases", {})

    if not isinstance(associated, dict):
        return []

    rows = associated.get("rows", [])

    if not isinstance(rows, list):
        return []

    records: list[EvidenceRecord] = []

    for row in rows:
        if not isinstance(row, dict):
            continue

        disease = row.get("disease", {})

        if not isinstance(disease, dict):
            continue

        disease_id = str(
            disease.get("id", "")
        ).strip()

        disease_name = str(
            disease.get("name", "")
        ).strip()

        score = row.get("score")

        if not disease_id or not disease_name:
            continue

        if not isinstance(score, (int, float)):
            continue

        normalized_name = disease_name.lower()

        if not any(
            term in normalized_name
            for term in accepted_terms
        ):
            continue

        records.append(
            EvidenceRecord(
                gene_id=symbol,
                evidence_type=EvidenceType.CANCER_ASSOCIATION,
                source="Open Targets Platform",
                source_version=source_version.strip(),
                evidence_id=disease_id,
                description=(
                    f"{disease_name}; "
                    f"Open Targets association score={score:.6f}"
                ),
                cancer_cohort=cancer_cohort,
                retrieved_at=retrieved_at,
                source_url=(
                    "https://platform.opentargets.org/"
                    f"disease/{disease_id}"
                ),
            )
        )

    records.sort(
        key=lambda record: (
            record.evidence_id or "",
            record.description,
        )
    )

    return records