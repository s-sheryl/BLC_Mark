"""
BLC Mark - UCSC Xena Client

Purpose:
    Provide the small, explicit interface BLC Mark uses to
    communicate with a UCSC Xena hub.

Responsibilities:
    - Establish a reusable HTTP session.
    - Verify that the configured Xena hub is reachable and responds
      correctly to a basic Xena query.
    - Verify that a configured dataset exists on the hub.
    - Resolve a verified direct-download URL for a dataset.

Scope:
    This module performs network discovery and verification only.

    It does NOT:
        - download complete datasets;
        - write dataset files;
        - preprocess expression matrices;
        - infer biological metadata;
        - perform differential expression analysis.

    Actual file downloading belongs to src.download_manager.
"""

from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import quote

import requests


XENA_CLIENT_VERSION = "1.0"

DEFAULT_XENA_HOST = "https://tcga.xenahubs.net"
DEFAULT_TIMEOUT_SECONDS = 30

__all__ = [
    "XenaClient",
    "XenaClientConfig",
    "DEFAULT_XENA_HOST",
    "DEFAULT_TIMEOUT_SECONDS",
    "XENA_CLIENT_VERSION",
]


@dataclass(frozen=True)
class XenaClientConfig:
    """Connection configuration for one UCSC Xena hub."""

    host: str = DEFAULT_XENA_HOST
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS

    def __post_init__(self) -> None:
        if not isinstance(self.host, str) or not self.host.strip():
            raise ValueError("'host' must be a non-empty string.")

        if not self.host.startswith(("http://", "https://")):
            raise ValueError(
                "'host' must begin with 'http://' or 'https://'."
            )

        if (
            isinstance(self.timeout_seconds, bool)
            or not isinstance(self.timeout_seconds, int)
            or self.timeout_seconds <= 0
        ):
            raise ValueError(
                "'timeout_seconds' must be a positive integer."
            )


class XenaClient:
    """Client for verified interaction with a UCSC Xena hub."""

    def __init__(self, config: XenaClientConfig) -> None:
        if not isinstance(config, XenaClientConfig):
            raise TypeError(
                "'config' must be an XenaClientConfig instance."
            )

        self.config = config
        self.session: requests.Session | None = None
        self.connected = False

    def connect(self) -> None:
        """Create the reusable HTTP session.

        Creating the session itself does not prove that the Xena hub
        is reachable. Use check_connection() for that verification.
        """
        if self.session is None:
            session = requests.Session()
            session.headers.update(
                {
                    "User-Agent": (
                        f"BLC Mark/{XENA_CLIENT_VERSION} "
                        "UCSC-Xena-Client"
                    )
                }
            )
            self.session = session

    def close(self) -> None:
        """Close the reusable HTTP session."""
        if self.session is not None:
            self.session.close()
            self.session = None

        self.connected = False

    def __enter__(self) -> "XenaClient":
        self.connect()
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def _ensure_session(self) -> requests.Session:
        """Return an initialized HTTP session."""
        if self.session is None:
            self.connect()

        assert self.session is not None
        return self.session

    def check_connection(self) -> bool:
        """Verify that the configured Xena hub responds to a Xena query.

        UCSC Xena hubs accept query expressions through POST requests
        to the /data/ endpoint. A simple arithmetic expression is used
        here because it verifies the Xena query service itself rather
        than merely checking whether the web server returns HTTP 200.

        Returns:
            True when the hub returns the expected Xena response.

        Raises:
            ConnectionError:
                If the hub cannot be reached or returns an unexpected
                response.
        """
        session = self._ensure_session()

        endpoint = f"{self.config.host.rstrip('/')}/data/"

        try:
            response = session.post(
                endpoint,
                data="(+ 1 2)",
                headers={"Content-Type": "text/plain"},
                timeout=self.config.timeout_seconds,
            )
            response.raise_for_status()
        except requests.RequestException as error:
            self.connected = False
            raise ConnectionError(
                "Could not connect to UCSC Xena hub "
                f"'{self.config.host}': {error}"
            ) from error

        result = response.text.strip()

        if result not in {"3", "3.0"}:
            self.connected = False
            raise ConnectionError(
                "The configured server responded, but did not behave "
                "like the expected UCSC Xena query service. "
                f"Received response: {result!r}."
            )

        self.connected = True
        return True

    @staticmethod
    def _validate_dataset_id(dataset_id: str) -> str:
        """Validate and normalize a caller-supplied Xena dataset ID."""
        if not isinstance(dataset_id, str):
            raise TypeError(
                "'dataset_id' must be a string, "
                f"got {type(dataset_id).__name__}."
            )

        dataset_id = dataset_id.strip()

        if not dataset_id:
            raise ValueError("'dataset_id' must not be empty.")

        if dataset_id.startswith(("/", "\\")):
            raise ValueError(
                "'dataset_id' must be a Xena dataset identifier, "
                "not an absolute path."
            )

        if ".." in dataset_id.split("/"):
            raise ValueError(
                "'dataset_id' must not contain parent-directory "
                "components ('..')."
            )

        return dataset_id

    def _candidate_download_url(self, dataset_id: str) -> str:
        """Construct the canonical Xena download endpoint.

        The returned URL is only a candidate. Public callers never
        receive it without search_dataset() first verifying that the
        remote resource actually exists.
        """
        validated_id = self._validate_dataset_id(dataset_id)

        # Preserve "/" because Xena dataset identifiers use path-like
        # components such as TCGA.BRCA.sampleMap/HiSeqV2.
        encoded_id = quote(validated_id, safe="/._-")

        return (
            f"{self.config.host.rstrip('/')}"
            f"/download/{encoded_id}.gz"
        )

    def _verify_remote_file(self, url: str) -> requests.Response:
        """Verify that a remote downloadable resource exists.

        HEAD is attempted first to avoid downloading dataset content.
        Some HTTP servers do not support HEAD correctly, so a one-byte
        ranged GET is used as a conservative fallback.
        """
        session = self._ensure_session()

        try:
            response = session.head(
                url,
                allow_redirects=True,
                timeout=self.config.timeout_seconds,
            )

            if response.status_code in {405, 501}:
                response.close()
                response = session.get(
                    url,
                    headers={"Range": "bytes=0-0"},
                    stream=True,
                    allow_redirects=True,
                    timeout=self.config.timeout_seconds,
                )

            response.raise_for_status()
            return response

        except requests.RequestException as error:
            raise LookupError(
                f"Xena dataset resource could not be verified at "
                f"'{url}': {error}"
            ) from error

    def search_dataset(self, dataset_id: str) -> dict:
        """Verify that a dataset exists on the configured Xena hub.

        This method deliberately performs a real remote check rather
        than trusting the registry entry or fabricating metadata.

        Args:
            dataset_id:
                Xena dataset identifier, for example
                "TCGA.BRCA.sampleMap/HiSeqV2".

        Returns:
            Dictionary containing the verified dataset identifier,
            hub, resolved URL, and available HTTP metadata.

        Raises:
            LookupError:
                If the dataset resource cannot be verified.
        """
        validated_id = self._validate_dataset_id(dataset_id)

        if not self.connected:
            self.check_connection()

        download_url = self._candidate_download_url(validated_id)
        response = self._verify_remote_file(download_url)

        try:
            content_length_header = response.headers.get(
                "Content-Length"
            )

            content_length = (
                int(content_length_header)
                if content_length_header
                and content_length_header.isdigit()
                else None
            )

            return {
                "dataset_id": validated_id,
                "xena_host": self.config.host.rstrip("/"),
                "download_url": str(response.url),
                "verified": True,
                "status_code": response.status_code,
                "content_type": response.headers.get("Content-Type"),
                "content_length": content_length,
                "content_encoding": response.headers.get(
                    "Content-Encoding"
                ),
            }
        finally:
            response.close()

    def resolve_download_url(self, dataset_id: str) -> str:
        """Return a verified direct-download URL for a Xena dataset.

        The URL is not returned solely because it can be constructed.
        search_dataset() performs a real network verification first.

        Args:
            dataset_id:
                Xena dataset identifier.

        Returns:
            Verified directly downloadable URL.

        Raises:
            LookupError:
                If the resource cannot be verified.
        """
        metadata = self.search_dataset(dataset_id)

        download_url = metadata.get("download_url")

        if not isinstance(download_url, str) or not download_url:
            raise LookupError(
                "Xena dataset verification succeeded but no valid "
                f"download URL was returned for '{dataset_id}'."
            )

        return download_url


if __name__ == "__main__":
    print("BLC Mark Xena Client")
    print("-" * 40)

    config = XenaClientConfig()
    client = XenaClient(config)

    print(f"Version : {XENA_CLIENT_VERSION}")
    print(f"Host    : {client.config.host}")
    print()

    try:
        print("Checking Xena connection...")
        client.check_connection()
        print("Connection: OK")

        # Verification only. This does not download the dataset.
        dataset_id = "TCGA.BRCA.sampleMap/HiSeqV2"

        print(f"Checking dataset: {dataset_id}")
        metadata = client.search_dataset(dataset_id)

        print("Dataset: VERIFIED")
        print(f"Download URL: {metadata['download_url']}")

    except (ConnectionError, LookupError) as error:
        print(f"Verification failed: {error}")

    finally:
        client.close()