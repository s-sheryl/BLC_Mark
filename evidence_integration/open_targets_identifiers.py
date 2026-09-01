"""Open Targets identifier resolution for BLC Mark Phase 4."""

import json
from urllib.request import Request, urlopen


OPEN_TARGETS_SEARCH_URL = (
    "https://api.platform.opentargets.org/api/v4/graphql"
)


def _post_graphql(query: str) -> dict:
    """Send a GraphQL request and return decoded JSON."""
    if not isinstance(query, str) or not query.strip():
        raise ValueError("query must be a non-empty string.")

    payload = json.dumps(
        {"query": query}
    ).encode("utf-8")

    request = Request(
        OPEN_TARGETS_SEARCH_URL,
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


def resolve_gene_symbol_to_ensembl(
    gene_symbol: str,
) -> str | None:
    """
    Resolve an exact human gene symbol to an Ensembl gene identifier.

    Returns None when no exact target match is found.
    """
    if not isinstance(gene_symbol, str) or not gene_symbol.strip():
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    symbol = gene_symbol.strip()

    query = f"""
    query {{
      search(
        queryString: "{symbol}"
        entityNames: ["target"]
        page: {{index: 0, size: 20}}
      ) {{
        hits {{
          id
          entity
          name
        }}
      }}
    }}
    """

    payload = _post_graphql(query)

    hits = (
        payload
        .get("data", {})
        .get("search", {})
        .get("hits", [])
    )

    if not isinstance(hits, list):
        return None

    exact_matches: list[str] = []

    for hit in hits:
        if not isinstance(hit, dict):
            continue

        entity = str(
            hit.get("entity", "")
        ).strip()

        name = str(
            hit.get("name", "")
        ).strip()

        identifier = str(
            hit.get("id", "")
        ).strip()

        if entity != "target":
            continue

        if name != symbol:
            continue

        if not identifier.startswith("ENSG"):
            continue

        exact_matches.append(identifier)

    exact_matches = sorted(set(exact_matches))

    if len(exact_matches) == 1:
        return exact_matches[0]

    return None