"""NCBI Gene evidence retrieval for BLC Mark Phase 4."""

import json
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .models import EvidenceRecord, EvidenceType


ESEARCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
)

ESUMMARY_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
)

HUMAN_TAXONOMY_ID = "9606"

# NCBI permits up to 3 requests/second without an API key.
# 0.40 seconds between requests keeps BLC Mark below that limit.
NCBI_REQUEST_INTERVAL_SECONDS = 0.40

# Transient network failures should not abort an entire Phase 4 run.
NCBI_MAX_RETRIES = 5
NCBI_INITIAL_BACKOFF_SECONDS = 1.0


_last_request_time = 0.0


def _respect_rate_limit() -> None:
    """Ensure NCBI requests remain below the unauthenticated rate limit."""
    global _last_request_time

    now = time.monotonic()

    elapsed = now - _last_request_time

    remaining = (
        NCBI_REQUEST_INTERVAL_SECONDS
        - elapsed
    )

    if remaining > 0:
        time.sleep(remaining)

    _last_request_time = time.monotonic()


def _get_json(
    url: str,
    params: dict,
) -> dict:
    """
    Perform a rate-limited NCBI HTTP GET request.

    Transient HTTP/network errors are retried with exponential backoff.
    Permanent failures are raised explicitly.
    """
    query = urlencode(params)

    request_url = f"{url}?{query}"

    backoff = NCBI_INITIAL_BACKOFF_SECONDS

    for attempt in range(
        1,
        NCBI_MAX_RETRIES + 1,
    ):
        _respect_rate_limit()

        request = Request(
            request_url,
            headers={
                "Accept": "application/json",
                "User-Agent": "BLC-Mark/1.0",
            },
            method="GET",
        )

        try:
            with urlopen(
                request,
                timeout=30,
            ) as response:
                return json.loads(
                    response
                    .read()
                    .decode("utf-8")
                )

        except HTTPError as exc:
            retryable = (
                exc.code == 429
                or 500 <= exc.code < 600
            )

            if (
                not retryable
                or attempt == NCBI_MAX_RETRIES
            ):
                raise

        except (
            URLError,
            ConnectionResetError,
            TimeoutError,
            OSError,
        ):
            if attempt == NCBI_MAX_RETRIES:
                raise

        time.sleep(backoff)

        backoff *= 2

    raise RuntimeError(
        "NCBI request failed after retry loop."
    )


def find_human_gene_id(
    gene_symbol: str,
) -> str | None:
    """
    Resolve an exact human gene symbol to an NCBI Gene ID.

    Returns None when no matching human Gene record is found.
    """
    if (
        not isinstance(gene_symbol, str)
        or not gene_symbol.strip()
    ):
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    symbol = gene_symbol.strip()

    payload = _get_json(
        ESEARCH_URL,
        {
            "db": "gene",
            "term": (
                f"{symbol}[sym] "
                f"AND {HUMAN_TAXONOMY_ID}[taxid]"
            ),
            "retmode": "json",
            "retmax": 5,
            "sort": "relevance",
            "tool": "BLC_Mark",
        },
    )

    ids = (
        payload
        .get("esearchresult", {})
        .get("idlist", [])
    )

    if not ids:
        return None

    return str(ids[0])


def fetch_gene_summary(
    gene_id: str,
) -> dict | None:
    """Retrieve the NCBI Gene document summary for a Gene ID."""
    if (
        not isinstance(gene_id, str)
        or not gene_id.strip()
    ):
        raise ValueError(
            "gene_id must be a non-empty string."
        )

    uid = gene_id.strip()

    payload = _get_json(
        ESUMMARY_URL,
        {
            "db": "gene",
            "id": uid,
            "retmode": "json",
            "tool": "BLC_Mark",
        },
    )

    result = payload.get(
        "result",
        {},
    )

    summary = result.get(uid)

    if not isinstance(summary, dict):
        return None

    return summary


def collect_functional_evidence(
    gene_symbol: str,
    *,
    retrieved_at: str,
) -> EvidenceRecord | None:
    """
    Collect NCBI Gene functional evidence for one human candidate gene.

    No evidence record is fabricated when the gene cannot be resolved,
    the returned record does not match the requested symbol, or NCBI
    provides no descriptive annotation.
    """
    if (
        not isinstance(gene_symbol, str)
        or not gene_symbol.strip()
    ):
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    symbol = gene_symbol.strip()

    gene_id = find_human_gene_id(
        symbol
    )

    if gene_id is None:
        return None

    summary = fetch_gene_summary(
        gene_id
    )

    if summary is None:
        return None

    ncbi_symbol = str(
        summary.get(
            "name",
            "",
        )
    ).strip()

    organism_data = summary.get(
        "organism",
        {},
    )

    if isinstance(
        organism_data,
        dict,
    ):
        organism = str(
            organism_data.get(
                "scientificname",
                "",
            )
        ).strip()
    else:
        organism = ""

    if (
        organism
        and organism != "Homo sapiens"
    ):
        return None

    if (
        ncbi_symbol
        and ncbi_symbol != symbol
    ):
        return None

    description = str(
        summary.get(
            "summary",
            "",
        )
    ).strip()

    if not description:
        description = str(
            summary.get(
                "description",
                "",
            )
        ).strip()

    if not description:
        return None

    return EvidenceRecord(
        gene_id=symbol,
        evidence_type=(
            EvidenceType.FUNCTIONAL
        ),
        source="NCBI Gene",
        source_version=None,
        evidence_id=(
            f"NCBI_GENE:{gene_id}"
        ),
        description=description,
        retrieved_at=retrieved_at,
        source_url=(
            "https://www.ncbi.nlm.nih.gov/"
            f"gene/{gene_id}"
        ),
    )