"""
Purpose:
    Provide the single source of truth for BLC Mark's project
    configuration -- paths, network defaults, and behavioral defaults
    that other modules read from instead of defining duplicate
    constants.

Responsibilities:
    - Define ProjectConfig, an immutable dataclass containing the
      project's directory paths and shared network/behavior defaults.
    - Validate configuration values at construction time.
    - Derive sensible default paths from the project root.
    - Expose stable Version 1.0 configuration constants.

Scope:
    This module represents configuration in memory only.

    It does not:
    - create directories,
    - read datasets,
    - download files,
    - configure logging,
    - make network requests,
    - validate biological data,
    - perform preprocessing,
    - perform quality control,
    - read environment variables,
    - parse configuration files,
    - perform file hashing.

    Constructing a ProjectConfig has no filesystem or network side
    effects.

Version:
    Frozen Version 1.0 configuration API.

    Field names, public constants, and validation behavior are intended
    to remain stable so dependent modules can rely on them without
    future breaking changes.
"""

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

from src.exceptions import ConfigurationError


CONFIGURATION_MANAGER_VERSION = "1.0"

DEFAULT_PROJECT_ROOT: Path = Path(__file__).resolve().parents[1]

DEFAULT_XENA_HOST: str = "https://tcga.xenahubs.net"
DEFAULT_TIMEOUT_SECONDS: int = 30
DEFAULT_MAX_RETRIES: int = 3
DEFAULT_CHUNK_SIZE_BYTES: int = 1024 * 1024
DEFAULT_BACKOFF_BASE_SECONDS: float = 1.0
DEFAULT_VERIFY_CHECKSUMS_BY_DEFAULT: bool = False


__all__ = [
    "CONFIGURATION_MANAGER_VERSION",
    "DEFAULT_PROJECT_ROOT",
    "DEFAULT_XENA_HOST",
    "DEFAULT_TIMEOUT_SECONDS",
    "DEFAULT_MAX_RETRIES",
    "DEFAULT_CHUNK_SIZE_BYTES",
    "DEFAULT_BACKOFF_BASE_SECONDS",
    "DEFAULT_VERIFY_CHECKSUMS_BY_DEFAULT",
    "ProjectConfig",
]


@dataclass(frozen=True)
class ProjectConfig:
    """Immutable, validated configuration shared across BLC Mark.

    Every field has a sensible default, so ``ProjectConfig()`` produces
    a complete usable configuration.

    Because the dataclass is frozen, configuration cannot be mutated
    after construction. To change configuration, construct a new
    ProjectConfig instance.

    Attributes:
        project_root:
            Root directory of the BLC Mark project.

        data_dir:
            Root directory for dataset-related data.

        downloads_dir:
            Directory for downloaded raw files.

        processed_dir:
            Directory for preprocessed, analysis-ready data.

        results_dir:
            Directory for analysis outputs and reports.

        logs_dir:
            Directory for log files.

        temp_dir:
            Directory for temporary or scratch files.

        xena_host:
            Default UCSC Xena data hub URL.

        timeout_seconds:
            Default network timeout in seconds.

        max_retries:
            Maximum number of retries for retryable network
            operations.

        chunk_size_bytes:
            Default streaming chunk size used by download and
            hashing infrastructure.

        backoff_base_seconds:
            Base delay used for exponential retry backoff.

        verify_checksums_by_default:
            Default policy for checksum verification.
    """

    project_root: Path = DEFAULT_PROJECT_ROOT

    data_dir: Path | None = None
    downloads_dir: Path | None = None
    processed_dir: Path | None = None
    results_dir: Path | None = None
    logs_dir: Path | None = None
    temp_dir: Path | None = None

    xena_host: str = DEFAULT_XENA_HOST

    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS
    max_retries: int = DEFAULT_MAX_RETRIES
    chunk_size_bytes: int = DEFAULT_CHUNK_SIZE_BYTES
    backoff_base_seconds: float = DEFAULT_BACKOFF_BASE_SECONDS

    verify_checksums_by_default: bool = DEFAULT_VERIFY_CHECKSUMS_BY_DEFAULT

    def __post_init__(self) -> None:
        """Validate and complete the configuration.

        Raises:
            TypeError:
                If a field has an invalid Python type.

            ConfigurationError:
                If a field has an invalid value.
        """

        self._require_path("project_root", self.project_root)
        self._validate_absolute_path(
            "project_root",
            self.project_root,
        )

        self._require_optional_path("data_dir", self.data_dir)
        self._require_optional_path("downloads_dir", self.downloads_dir)
        self._require_optional_path("processed_dir", self.processed_dir)
        self._require_optional_path("results_dir", self.results_dir)
        self._require_optional_path("logs_dir", self.logs_dir)
        self._require_optional_path("temp_dir", self.temp_dir)

        data_dir = (
            self.data_dir
            if self.data_dir is not None
            else self.project_root / "data"
        )

        downloads_dir = (
            self.downloads_dir
            if self.downloads_dir is not None
            else data_dir / "downloads"
        )

        processed_dir = (
            self.processed_dir
            if self.processed_dir is not None
            else data_dir / "processed"
        )

        results_dir = (
            self.results_dir
            if self.results_dir is not None
            else self.project_root / "results"
        )

        logs_dir = (
            self.logs_dir
            if self.logs_dir is not None
            else self.project_root / "logs"
        )

        temp_dir = (
            self.temp_dir
            if self.temp_dir is not None
            else self.project_root / "tmp"
        )

        object.__setattr__(self, "data_dir", data_dir)
        object.__setattr__(self, "downloads_dir", downloads_dir)
        object.__setattr__(self, "processed_dir", processed_dir)
        object.__setattr__(self, "results_dir", results_dir)
        object.__setattr__(self, "logs_dir", logs_dir)
        object.__setattr__(self, "temp_dir", temp_dir)

        for field_name in (
            "data_dir",
            "downloads_dir",
            "processed_dir",
            "results_dir",
            "logs_dir",
            "temp_dir",
        ):
            path_value = getattr(self, field_name)

            self._require_path(field_name, path_value)
            self._validate_absolute_path(
                field_name,
                path_value,
            )

        self._validate_xena_host()

        self._validate_positive_int(
            "timeout_seconds",
            self.timeout_seconds,
        )

        self._validate_non_negative_int(
            "max_retries",
            self.max_retries,
        )

        self._validate_positive_int(
            "chunk_size_bytes",
            self.chunk_size_bytes,
        )

        self._validate_positive_number(
            "backoff_base_seconds",
            self.backoff_base_seconds,
        )

        if not isinstance(
            self.verify_checksums_by_default,
            bool,
        ):
            raise TypeError(
                "'verify_checksums_by_default' must be a bool, "
                f"got {type(self.verify_checksums_by_default).__name__}."
            )

    @staticmethod
    def _require_path(
        field_name: str,
        value: object,
    ) -> None:
        """Require a value to be a pathlib.Path."""

        if not isinstance(value, Path):
            raise TypeError(
                f"'{field_name}' must be a pathlib.Path, "
                f"got {type(value).__name__}."
            )

    @staticmethod
    def _require_optional_path(
        field_name: str,
        value: object,
    ) -> None:
        """Require a value to be a pathlib.Path or None."""

        if value is not None and not isinstance(value, Path):
            raise TypeError(
                f"'{field_name}' must be a pathlib.Path or None, "
                f"got {type(value).__name__}."
            )

    @staticmethod
    def _validate_absolute_path(
        field_name: str,
        path_value: Path,
    ) -> None:
        """Require a path to be absolute."""

        if not path_value.is_absolute():
            raise ConfigurationError(
                f"'{field_name}' must be an absolute path, "
                f"got {path_value!r}."
            )

    @staticmethod
    def _validate_positive_int(
        field_name: str,
        value: object,
    ) -> None:
        """Require a strictly positive integer."""

        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(
                f"'{field_name}' must be an int, "
                f"got {type(value).__name__}."
            )

        if value <= 0:
            raise ConfigurationError(
                f"'{field_name}' must be positive, got {value!r}."
            )

    @staticmethod
    def _validate_non_negative_int(
        field_name: str,
        value: object,
    ) -> None:
        """Require a non-negative integer."""

        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(
                f"'{field_name}' must be an int, "
                f"got {type(value).__name__}."
            )

        if value < 0:
            raise ConfigurationError(
                f"'{field_name}' must be non-negative, got {value!r}."
            )

    @staticmethod
    def _validate_positive_number(
        field_name: str,
        value: object,
    ) -> None:
        """Require a strictly positive int or float."""

        if isinstance(value, bool) or not isinstance(
            value,
            (int, float),
        ):
            raise TypeError(
                f"'{field_name}' must be an int or float, "
                f"got {type(value).__name__}."
            )

        if value <= 0:
            raise ConfigurationError(
                f"'{field_name}' must be positive, got {value!r}."
            )

    def _validate_xena_host(self) -> None:
        """Validate the configured Xena host URL."""

        if not isinstance(self.xena_host, str):
            raise TypeError(
                f"'xena_host' must be a str, "
                f"got {type(self.xena_host).__name__}."
            )

        if not self.xena_host.strip():
            raise ConfigurationError(
                "'xena_host' must not be empty."
            )

        parsed = urlparse(self.xena_host)

        if parsed.scheme not in {"http", "https"}:
            raise ConfigurationError(
                "'xena_host' must use the http:// or https:// "
                f"scheme, got {self.xena_host!r}."
            )

        if not parsed.netloc:
            raise ConfigurationError(
                "'xena_host' must contain a valid host, "
                f"got {self.xena_host!r}."
            )


if __name__ == "__main__":
    print("BLC Mark Configuration Manager")
    print(f"Version: {CONFIGURATION_MANAGER_VERSION}\n")

    config = ProjectConfig()

    print("Paths:")
    print(f"  project_root:  {config.project_root}")
    print(f"  data_dir:      {config.data_dir}")
    print(f"  downloads_dir: {config.downloads_dir}")
    print(f"  processed_dir: {config.processed_dir}")
    print(f"  results_dir:   {config.results_dir}")
    print(f"  logs_dir:      {config.logs_dir}")
    print(f"  temp_dir:      {config.temp_dir}")

    print("\nNetwork:")
    print(f"  xena_host:            {config.xena_host}")
    print(f"  timeout_seconds:      {config.timeout_seconds}")
    print(f"  max_retries:          {config.max_retries}")
    print(f"  chunk_size_bytes:     {config.chunk_size_bytes}")
    print(f"  backoff_base_seconds: {config.backoff_base_seconds}")

    print("\nBehavior:")
    print(
        "  verify_checksums_by_default: "
        f"{config.verify_checksums_by_default}"
    )

    print("\nValidation demonstrations:")

    try:
        ProjectConfig(
            project_root=Path("relative/path"),
        )
    except ConfigurationError as error:
        print(f"  Relative project_root rejected: {error}")

    try:
        ProjectConfig(
            timeout_seconds=0,
        )
    except ConfigurationError as error:
        print(f"  Invalid timeout rejected: {error}")

    try:
        ProjectConfig(
            max_retries=-1,
        )
    except ConfigurationError as error:
        print(f"  Negative max_retries rejected: {error}")

    try:
        ProjectConfig(
            xena_host="not-a-url",
        )
    except ConfigurationError as error:
        print(f"  Invalid Xena host rejected: {error}")

    try:
        ProjectConfig(
            data_dir="data",
        )
    except TypeError as error:
        print(f"  Invalid data_dir type rejected: {error}")

    print("\nFrozen dataclass demonstration:")

    try:
        config.timeout_seconds = 60
    except (AttributeError, TypeError) as error:
        print(f"  Mutation correctly rejected: {error}")

    print("\nProjectConfig constructed and validated successfully.")