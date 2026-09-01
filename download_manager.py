"""
Purpose:
    Download files from HTTP or HTTPS URLs to local disk, safely,
    efficiently, atomically, and resumably -- nothing else.

Responsibilities:
    - Stream a file from a URL to a user-specified local path, in
      fixed-size chunks, over a reused HTTP connection.
    - Write downloads atomically: to a ``.part`` file first, verified,
      then renamed into place -- an interrupted download can never be
      mistaken for a complete one.
    - Resume an interrupted download via HTTP Range requests when the
      server supports it, falling back to a full restart when it does
      not.
    - Verify a downloaded file is non-empty, matches the server's
      reported Content-Length (when meaningful), and optionally
      matches a caller-supplied SHA-256 checksum.
    - Retry only genuinely transient failures (timeouts, connection
      errors, and a defined set of temporary HTTP status codes) with
      exponential backoff and jitter -- never retry a 404, 403, or a
      malformed URL.
    - Report progress via an optional callback. This module never
      prints anything.
    - Return structured DownloadResult objects rather than raising for
      ordinary download failures, so batch callers can inspect every
      outcome without wrapping every call in a try/except.

Scope:
    This module only downloads files. It does not inspect expression
    matrices, does not validate dataset structure, does not register
    datasets, does not preprocess or normalize data, and has no
    knowledge of TCGA, GEO, Xena, cancer types, or genes of any kind.
    It never creates JSON, never persists metadata, and never talks
    to DatasetManager or any other module.

Version:
    This is a frozen Version 1.0 implementation. Its public interface
    (DownloadResult, DownloadStatus, the exception hierarchy, and
    DownloadManager's public methods) is intended to remain stable.
"""

import hashlib
import logging
import random
import re
import tempfile
import threading
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable

import requests

logger = logging.getLogger(__name__)

DOWNLOAD_MANAGER_VERSION: str = "1.0"

DEFAULT_TIMEOUT_SECONDS: int = 30
DEFAULT_MAX_RETRIES: int = 3
DEFAULT_CHUNK_SIZE_BYTES: int = 1024 * 1024  # 1 MB
DEFAULT_BACKOFF_BASE_SECONDS: float = 1.0
MAX_BACKOFF_JITTER_SECONDS: float = 0.5
PARTIAL_FILE_SUFFIX: str = ".part"

# HTTP status codes considered temporary and therefore worth retrying.
# Anything outside this set (404, 403, 401, 400, 410, ...) is treated
# as a definitive answer from the server that retrying will not change.
from http import HTTPStatus
TRANSIENT_HTTP_STATUS_CODES: frozenset[int] = frozenset(
    {
        HTTPStatus.REQUEST_TIMEOUT,
        HTTPStatus.TOO_EARLY,
        HTTPStatus.TOO_MANY_REQUESTS,
        HTTPStatus.INTERNAL_SERVER_ERROR,
        HTTPStatus.BAD_GATEWAY,
        HTTPStatus.SERVICE_UNAVAILABLE,
        HTTPStatus.GATEWAY_TIMEOUT,
    }
)

_CONTENT_RANGE_TOTAL_PATTERN: re.Pattern[str] = re.compile(r"/(\d+)\s*$")

ProgressCallback = Callable[[int, "int | None"], None]

__all__ = [
    "DOWNLOAD_MANAGER_VERSION",
    "DownloadError",
    "InvalidURLError",
    "DownloadVerificationError",
    "ChecksumMismatchError",
    "DownloadStatus",
    "DownloadResult",
    "DownloadManager",
]


class DownloadError(Exception):
    """Base exception for all download failures raised by this module."""


class InvalidURLError(DownloadError):
    """Raised when a URL is malformed or uses an unsupported scheme."""


class DownloadVerificationError(DownloadError):
    """Raised when a downloaded file fails post-download verification.

    Covers a missing file, an empty file, or a size that does not
    match the server-reported Content-Length.
    """


class ChecksumMismatchError(DownloadVerificationError):
    """Raised when a downloaded file's SHA-256 does not match expectations."""


class DownloadStatus(Enum):
    """Outcome of a single download attempt.

    Members:
        SUCCESS: The file was downloaded (or already existed and was
            verified) and passed every requested verification step.
        SKIPPED: The destination file already existed, overwrite was
            not requested, and (if a checksum was supplied) it
            matched -- no network request was made.
        FAILED: The download could not be completed, or the result
            failed verification. See DownloadResult.error_message for
            details.
    """

    SUCCESS = "success"
    SKIPPED = "skipped"
    FAILED = "failed"


@dataclass(frozen=True)
class DownloadResult:
    """Structured outcome of a single download_file() call.

    Attributes:
        url: The source URL that was requested.
        destination: The local path the file was (or would have been)
            written to.
        status: The overall outcome of this download attempt.
        skipped: True if no network request was made because the
            destination already existed and overwrite was not
            requested.
        file_size: Size of the resulting file in bytes, or None if the
            download failed before any file existed.
        sha256_checksum: The file's SHA-256 hex digest, computed only
            when a `checksum` argument was supplied to download_file()
            (so a caller who never asked for verification never pays
            the cost of hashing). None otherwise.
        download_time_seconds: Wall-clock time spent inside
            download_file(), including any retries and backoff delay.
        attempts_used: Number of network attempts made. Zero for a
            skipped download or a download that failed before any
            request was sent (e.g. an invalid URL).
        error_message: Human-readable description of what went wrong,
            or None if status is SUCCESS or SKIPPED.
    """

    url: str
    destination: Path
    status: DownloadStatus
    skipped: bool = False
    file_size: int | None = None
    sha256_checksum: str | None = None
    download_time_seconds: float = 0.0
    attempts_used: int = 0
    error_message: str | None = None


class DownloadManager:
    """Downloads files from HTTP/HTTPS URLs to local disk.

    This class has exactly one job: get bytes from a URL onto disk
    reliably, atomically, and resumably. It knows nothing about
    datasets, cancer types, or file formats -- every URL and
    destination path is treated as opaque.

    Thread safety:
        A DownloadManager instance holds no mutable state that is
        written during a download beyond the shared `requests.Session`
        -- every other value used while downloading (attempt counters,
        byte counts, partial-file paths) is a local variable scoped to
        one download_file() call. `requests.Session` is designed to
        support concurrent requests from multiple threads for the
        common case of independent GET requests, which is the only
        way this class uses it. Multiple threads may therefore safely
        call download_file() or download_multiple() on the same
        DownloadManager instance concurrently, provided they target
        different destination paths.

    Attributes:
        timeout_seconds: Per-request timeout, in seconds.
        max_retries: Maximum number of attempts for a single download
            before it is treated as failed.
        chunk_size_bytes: Size, in bytes, of each chunk read from the
            response stream, written to disk, and hashed during
            checksum verification.
        backoff_base_seconds: Base delay, in seconds, for exponential
            backoff between retries.
        session: The persistent `requests.Session` used for every HTTP
            request made by this instance.
    """

    def __init__(
        self,
        timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
        max_retries: int = DEFAULT_MAX_RETRIES,
        chunk_size_bytes: int = DEFAULT_CHUNK_SIZE_BYTES,
        backoff_base_seconds: float = DEFAULT_BACKOFF_BASE_SECONDS,
    ) -> None:
        """Store download configuration and open a persistent session.

        Args:
            timeout_seconds: Per-request timeout, in seconds.
            max_retries: Maximum number of attempts for a single
                download before it is treated as failed.
            chunk_size_bytes: Size, in bytes, of each chunk read from
                the response stream and written to disk.
            backoff_base_seconds: Base delay, in seconds, used to
                compute exponential backoff between retries.

        Raises:
            ValueError: If any argument is not a positive value.
        """
        if timeout_seconds <= 0:
            raise ValueError(f"timeout_seconds must be positive, got {timeout_seconds}.")
        if max_retries <= 0:
            raise ValueError(f"max_retries must be positive, got {max_retries}.")
        if chunk_size_bytes <= 0:
            raise ValueError(f"chunk_size_bytes must be positive, got {chunk_size_bytes}.")
        if backoff_base_seconds <= 0:
            raise ValueError(
                f"backoff_base_seconds must be positive, got {backoff_base_seconds}."
            )

        self.timeout_seconds: int = timeout_seconds
        self.max_retries: int = max_retries
        self.chunk_size_bytes: int = chunk_size_bytes
        self.backoff_base_seconds: float = backoff_base_seconds
        self.session: requests.Session = requests.Session()
        self._log_lock: threading.Lock = threading.Lock()

    def __enter__(self) -> "DownloadManager":
        """Support use as a context manager.

        Returns:
            This DownloadManager instance.
        """
        return self

    def __exit__(self, *_exc_info: object) -> None:
        """Close the underlying session when used as a context manager."""
        self.close()

    def close(self) -> None:
        """Close the persistent HTTP session.

        Raises:
            None.
        """
        self.session.close()

    # -- Public downloading API --------------------------------------------------

    def download_file(
        self,
        url: str,
        destination: Path,
        *,
        overwrite: bool = False,
        checksum: str | None = None,
        progress_callback: ProgressCallback | None = None,
    ) -> DownloadResult:
        """Download a single file from a URL to a local destination.

        The file is downloaded atomically: bytes are streamed to a
        temporary ``destination.part`` file, verified, and only then
        renamed to `destination`. If a download attempt fails for a
        transient reason (timeout, connection error, or a temporary
        HTTP status), the partial file is kept and the next attempt
        resumes from where it left off via an HTTP Range request, if
        the server honors it; if the server does not, the download
        restarts from the beginning automatically. If the download
        fails for a non-transient reason (404, 403, invalid URL, or
        exhausted retries), the partial file is removed and
        `destination` is left untouched.

        This method does not raise for ordinary download failures --
        every outcome, success or failure, is reported through the
        returned DownloadResult so callers (including
        download_multiple()) can inspect results uniformly. It raises
        only for programming errors in how it was called.

        Args:
            url: The HTTP or HTTPS URL to download from.
            destination: Local path to write the downloaded file to.
                Parent directories are created automatically.
            overwrite: If False (default), an existing file at
                `destination` is left in place and treated as a
                SKIPPED result. If True, any existing file is
                replaced.
            checksum: If provided, the downloaded (or existing, when
                skipped) file's SHA-256 hex digest must match this
                value exactly, or the result's status is FAILED. If
                None, no checksum is computed or verified.
            progress_callback: Optional callable invoked as
                `progress_callback(downloaded_bytes, total_bytes)`
                after each chunk is written. `total_bytes` is the
                server-reported total size as an int when known, or
                None otherwise. This module never prints progress
                itself.

        Returns:
            A DownloadResult describing the outcome.

        Raises:
            TypeError: If `url` or `destination` is not the expected
                type. This is the only case this method raises rather
                than reporting failure through the returned result,
                since it indicates a caller programming error rather
                than a download failure.
        """
        if not isinstance(url, str):
            raise TypeError(f"url must be a str, got {type(url).__name__}.")
        if not isinstance(destination, Path):
            raise TypeError(f"destination must be a Path, got {type(destination).__name__}.")

        start_time = time.perf_counter()
        self._log_event(logging.INFO, "download_started", url=url, destination=str(destination))

        try:
            self._validate_url_scheme(url)
        except InvalidURLError as error:
            return self._build_failed_result(url, destination, 0, start_time, str(error))

        if destination.exists() and not overwrite:
            return self._handle_skip(url, destination, checksum, start_time)

        try:
            self._verify_destination_writable(destination.parent)
        except DownloadError as error:
            return self._build_failed_result(url, destination, 0, start_time, str(error))

        partial_path = self._partial_path_for(destination)

        try:
            attempts_used = self._stream_with_retries(url, partial_path, progress_callback)
            self._verify_download(partial_path)

            sha256_checksum: str | None = None
            if checksum is not None:
                sha256_checksum = self._calculate_sha256(partial_path)
                self._verify_checksum(sha256_checksum, checksum, str(partial_path))

            partial_path.replace(destination)
            file_size = destination.stat().st_size
            elapsed = time.perf_counter() - start_time

            self._log_event(
                logging.INFO,
                "download_succeeded",
                url=url,
                destination=str(destination),
                attempts=attempts_used,
                file_size=file_size,
                elapsed_seconds=round(elapsed, 3),
            )

            return DownloadResult(
                url=url,
                destination=destination,
                status=DownloadStatus.SUCCESS,
                skipped=False,
                file_size=file_size,
                sha256_checksum=sha256_checksum,
                download_time_seconds=elapsed,
                attempts_used=attempts_used,
            )

        except DownloadError as error:
            self.remove_partial_download(partial_path)
            attempts_used = getattr(error, "attempts_used", self.max_retries)
            return self._build_failed_result(
                url, destination, attempts_used, start_time, str(error)
            )
        except OSError as error:
            self.remove_partial_download(partial_path)
            message = f"Disk write error while downloading '{url}': {error}"
            self._log_event(logging.ERROR, "download_failed", url=url, error=message)
            return self._build_failed_result(url, destination, 0, start_time, message)

    def download_multiple(
        self,
        downloads: list[tuple[str, Path]],
        *,
        overwrite: bool = False,
    ) -> dict[str, DownloadResult]:
        """Download several files, tolerating individual failures.

        Each (url, destination) pair is downloaded independently via
        download_file(). Since download_file() never raises for
        ordinary failures, every pair is always attempted regardless
        of how earlier ones turned out.

        Args:
            downloads: A list of (url, destination) tuples to
                download.
            overwrite: Passed through to download_file() for every
                pair.

        Returns:
            A dictionary mapping each URL to its DownloadResult.

        Raises:
            None.
        """
        results: dict[str, DownloadResult] = {}

        for url, destination in downloads:
            results[url] = self.download_file(url, destination, overwrite=overwrite)

        succeeded = sum(1 for result in results.values() if result.status == DownloadStatus.SUCCESS)
        self._log_event(
            logging.INFO,
            "download_multiple_complete",
            succeeded=succeeded,
            total=len(downloads),
        )

        return results
   
    def _verify_download(self, path: Path) -> bool:
        """Verify a downloaded file exists and is non-empty."""
        if not path.exists():
            raise DownloadVerificationError(
                f"Downloaded file not found: {path}"
            )

        file_size = path.stat().st_size

        if file_size == 0:
            raise DownloadVerificationError(
                f"Downloaded file is empty: {path}"
            )

        self._log_event(
            logging.INFO,
            "verification_passed",
            path=str(path),
            file_size=file_size,
        )

        return True

    def remove_partial_download(self, path: Path) -> None:
        """Delete a partially written or failed download, if present."""
        if path.exists():
            path.unlink()
            self._log_event(
                logging.INFO,
                "partial_removed",
                path=str(path),
            )

    def _handle_skip(
        self,
        url: str,
        destination: Path,
        checksum: str | None,
        start_time: float,
    ) -> DownloadResult:
        """Build the result for an already-existing destination file."""
        self._log_event(
            logging.INFO,
            "download_skipped",
            url=url,
            destination=str(destination),
        )

        file_size = destination.stat().st_size
        sha256_checksum: str | None = None

        if checksum is not None:
            sha256_checksum = self._calculate_sha256(destination)

            try:
                self._verify_checksum(
                    sha256_checksum,
                    checksum,
                    str(destination),
                )
            except ChecksumMismatchError as error:
                return self._build_failed_result(
                    url,
                    destination,
                    0,
                    start_time,
                    str(error),
                )

        return DownloadResult(
            url=url,
            destination=destination,
            status=DownloadStatus.SKIPPED,
            skipped=True,
            file_size=file_size,
            sha256_checksum=sha256_checksum,
            download_time_seconds=time.perf_counter() - start_time,
            attempts_used=0,
        )

    def _build_failed_result(
        self,
        url: str,
        destination: Path,
        attempts_used: int,
        start_time: float,
        error_message: str,
    ) -> DownloadResult:
        """Build a FAILED DownloadResult."""
        return DownloadResult(
            url=url,
            destination=destination,
            status=DownloadStatus.FAILED,
            skipped=False,
            download_time_seconds=time.perf_counter() - start_time,
            attempts_used=attempts_used,
            error_message=error_message,
        )

    def _validate_url_scheme(self, url: str) -> None:
        """Confirm that a URL uses HTTP or HTTPS."""
        if not url.lower().startswith(("http://", "https://")):
            raise InvalidURLError(
                f"Unsupported URL scheme for '{url}'. "
                "Only http:// and https:// URLs are supported."
            )

    def _verify_destination_writable(
        self,
        directory: Path,
    ) -> None:
        """Confirm the destination directory exists and is writable."""
        try:
            directory.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            message = (
                f"Cannot create destination directory "
                f"'{directory}': {error}"
            )
            self._log_event(
                logging.ERROR,
                "directory_creation_failed",
                error=message,
            )
            raise DownloadError(message) from error

        try:
            with tempfile.NamedTemporaryFile(
                dir=directory,
                delete=True,
            ):
                pass
        except OSError as error:
            message = (
                f"Destination directory '{directory}' "
                f"is not writable: {error}"
            )
            self._log_event(
                logging.ERROR,
                "directory_not_writable",
                error=message,
            )
            raise DownloadError(message) from error

    def _partial_path_for(self, destination: Path) -> Path:
        """Return the temporary .part path used during download."""
        return destination.with_name(
            destination.name + PARTIAL_FILE_SUFFIX
        )

    def _stream_with_retries(
        self,
        url: str,
        partial_path: Path,
        progress_callback: ProgressCallback | None,
    ) -> int:
        """Stream a file with retry and resume support."""
        last_error: Exception | None = None

        for attempt in range(1, self.max_retries + 1):
            resume_from = (
                partial_path.stat().st_size
                if partial_path.exists()
                else 0
            )

            self._log_event(
                logging.INFO,
                "download_attempt",
                url=url,
                attempt=attempt,
                max_attempts=self.max_retries,
                resume_from=resume_from,
            )

            try:
                self._download_once(
                    url,
                    partial_path,
                    progress_callback,
                    resume_from,
                )
                return attempt

            except requests.exceptions.RequestException as error:
                last_error = error

                if not self._is_retryable(error):
                    self.remove_partial_download(partial_path)

                    message = (
                        f"Non-retryable download failure for "
                        f"'{url}': {error}"
                    )

                    self._log_event(
                        logging.ERROR,
                        "download_failed_permanent",
                        url=url,
                        error=message,
                    )

                    raise self._error_with_attempts(
                        DownloadError(message),
                        attempt,
                    ) from error

                self._log_event(
                    logging.WARNING,
                    "download_attempt_failed",
                    url=url,
                    attempt=attempt,
                    max_attempts=self.max_retries,
                    error=str(error),
                )

                if attempt < self.max_retries:
                    delay_seconds = self._compute_backoff_delay(
                        attempt
                    )

                    self._log_event(
                        logging.INFO,
                        "retry_scheduled",
                        url=url,
                        delay_seconds=round(delay_seconds, 3),
                    )

                    time.sleep(delay_seconds)

            except DownloadVerificationError as error:
                last_error = error

                self._log_event(
                    logging.WARNING,
                    "download_attempt_failed",
                    url=url,
                    attempt=attempt,
                    max_attempts=self.max_retries,
                    error=str(error),
                )

                if attempt < self.max_retries:
                    delay_seconds = self._compute_backoff_delay(
                        attempt
                    )
                    time.sleep(delay_seconds)

        self.remove_partial_download(partial_path)

        message = (
            f"Download failed after {self.max_retries} "
            f"attempt(s) for '{url}': {last_error}"
        )

        self._log_event(
            logging.ERROR,
            "download_failed_exhausted",
            url=url,
            error=message,
        )

        raise self._error_with_attempts(
            DownloadError(message),
            self.max_retries,
        ) from last_error

    def _error_with_attempts(
        self,
        error: DownloadError,
        attempts_used: int,
    ) -> DownloadError:
        """Attach attempts_used to a DownloadError."""
        error.attempts_used = attempts_used  # type: ignore[attr-defined]
        return error

    def _is_retryable(
        self,
        error: requests.exceptions.RequestException,
    ) -> bool:
        """Return True when retrying a request may succeed."""
        if isinstance(
            error,
            (
                requests.exceptions.Timeout,
                requests.exceptions.ConnectionError,
            ),
        ):
            return True

        if isinstance(
            error,
            requests.exceptions.ChunkedEncodingError,
        ):
            return True

        if isinstance(
            error,
            requests.exceptions.HTTPError,
        ):
            response = error.response

            if (
                response is not None
                and response.status_code
                in TRANSIENT_HTTP_STATUS_CODES
            ):
                return True

            return False

        return False

    def _compute_backoff_delay(
        self,
        attempt: int,
    ) -> float:
        """Calculate exponential retry delay plus jitter."""
        exponential_delay = (
            self.backoff_base_seconds
            * (2 ** (attempt - 1))
        )

        jitter = random.uniform(
            0,
            MAX_BACKOFF_JITTER_SECONDS,
        )

        return exponential_delay + jitter

    def _download_once(
        self,
        url: str,
        partial_path: Path,
        progress_callback: ProgressCallback | None,
        resume_from: int,
    ) -> None:
        """Perform one streaming download attempt."""
        headers = (
            {"Range": f"bytes={resume_from}-"}
            if resume_from > 0
            else {}
        )

        with self.session.get(
            url,
            headers=headers,
            stream=True,
            timeout=self.timeout_seconds,
        ) as response:
            response.raise_for_status()

            is_resumed_response = (
                resume_from > 0
                and response.status_code == 206
            )

            write_mode = (
                "ab"
                if is_resumed_response
                else "wb"
            )

            downloaded_bytes = (
                resume_from
                if is_resumed_response
                else 0
            )

            if (
                resume_from > 0
                and not is_resumed_response
            ):
                self._log_event(
                    logging.INFO,
                    "resume_not_supported",
                    url=url,
                    detail=(
                        "Server did not honor Range request; "
                        "restarting from the beginning."
                    ),
                )

            total_bytes = self._determine_total_bytes(
                response,
                is_resumed_response,
            )

            with open(
                partial_path,
                write_mode,
            ) as file_handle:

                for chunk in response.iter_content(
                    chunk_size=self.chunk_size_bytes
                ):
                    if not chunk:
                        continue

                    file_handle.write(chunk)
                    downloaded_bytes += len(chunk)

                    if progress_callback is not None:
                        progress_callback(
                            downloaded_bytes,
                            total_bytes,
                        )

        actual_size = partial_path.stat().st_size

        if (
            total_bytes is not None
            and actual_size != total_bytes
        ):
            raise DownloadVerificationError(
                f"Downloaded size ({actual_size} bytes) "
                f"does not match expected total size "
                f"({total_bytes} bytes) for '{url}'."
            )

    def _determine_total_bytes(
        self,
        response: requests.Response,
        is_resumed_response: bool,
    ) -> int | None:
        """Determine expected full file size from HTTP headers."""
        if is_resumed_response:
            content_range = response.headers.get(
                "Content-Range"
            )

            if content_range:
                match = _CONTENT_RANGE_TOTAL_PATTERN.search(
                    content_range
                )

                if match:
                    return int(match.group(1))

            return None

        if response.headers.get("Content-Encoding"):
            return None

        content_length = response.headers.get(
            "Content-Length"
        )

        return (
            int(content_length)
            if content_length is not None
            else None
        )

    def _verify_checksum(
        self,
        actual_sha256: str,
        expected_sha256: str,
        path_label: str,
    ) -> None:
        """Verify SHA-256 checksum equality."""
        self._log_event(
            logging.INFO,
            "checksum_verification",
            path=path_label,
            expected=expected_sha256,
            actual=actual_sha256,
        )

        if (
            actual_sha256.lower()
            != expected_sha256.lower()
        ):
            raise ChecksumMismatchError(
                f"Checksum mismatch for '{path_label}': "
                f"expected {expected_sha256}, "
                f"got {actual_sha256}."
            )

    def _calculate_sha256(
        self,
        path: Path,
    ) -> str:
        """Calculate SHA-256 for a file."""
        digest = hashlib.sha256()

        with open(path, "rb") as file_handle:
            for chunk in iter(
                lambda: file_handle.read(
                    self.chunk_size_bytes
                ),
                b"",
            ):
                digest.update(chunk)

        return digest.hexdigest()

    def _log_event(
        self,
        level: int,
        event: str,
        **fields: object,
    ) -> None:
        """Emit a structured log entry."""
        with self._log_lock:
            field_text = " ".join(
                f"{key}={value}"
                for key, value in fields.items()
            )

            logger.log(
                level,
                "event=%s %s",
                event,
                field_text,
            )


if __name__ == "__main__":
    import tempfile as tempfile_module

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    print("BLC Mark Download Manager")
    print(f"Version: {DOWNLOAD_MANAGER_VERSION}\n")

    # Demonstration only: downloads a small, stable public text file
    # to show streaming, atomic write, retry, and verification
    # mechanics end to end. This is a generic file, not a
    # transcriptomic dataset -- no dataset source, cancer type, or
    # dataset name is implied, consistent with this module staying
    # fully generic.
    demonstration_url = "https://raw.githubusercontent.com/pandas-dev/pandas/main/LICENSE"

    def _report_progress(downloaded: int, total: int | None) -> None:
        if total is not None:
            print(f"  progress: {downloaded}/{total} bytes")
        else:
            print(f"  progress: {downloaded} bytes (total unknown)")

    with DownloadManager(timeout_seconds=30, max_retries=3) as manager:
        with tempfile_module.TemporaryDirectory() as scratch_dir:
            destination_path = Path(scratch_dir) / "downloads" / "demo_file.txt"

            result = manager.download_file(
                demonstration_url,
                destination_path,
                progress_callback=_report_progress,
            )

            print(f"\nStatus: {result.status.value}")
            print(f"File size: {result.file_size} bytes")
            print(f"Attempts used: {result.attempts_used}")
            print(f"Elapsed: {result.download_time_seconds:.3f}s")
            if result.error_message:
                print(f"Error: {result.error_message}")