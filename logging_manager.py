
"""
Purpose:
    Provide the single, central logging configuration for
    BLC Mark, so every module obtains its logger from here
    instead of configuring the logging framework independently.

Responsibilities:
    - Define LoggingConfig, an immutable dataclass describing how
      logging should be set up: level, format, and whether to log to
      console, file, or both.
    - Configure Python's built-in logging framework via
      configure_logging(), using explicit handler setup rather than
      logging.basicConfig(), and without leaving duplicate handlers
      behind if called more than once.
    - Provide get_logger() as a thin, direct wrapper around
      logging.getLogger() for modules to obtain a named logger once
      configuration is in place.

Scope:
    This module only configures Python's logging framework. It does
    not perform downloads, does not read or validate datasets, does
    not perform preprocessing, does not make network requests, and
    does not create analysis reports. The one file-system action it
    performs is creating the parent directory of the configured log
    file when file logging is enabled -- a direct, necessary
    consequence of the file-logging responsibility this module
    explicitly owns, not a side effect incidental to something else.

Version:
    This is a frozen Version 1.0 logging API. Its field names,
    function signatures, and behavior are intended to remain stable
    so dependent modules (dataset_manager.py, download_manager.py,
    xena_client.py, preprocessing_manager.py, quality_control.py,
    metadata_manager.py) can rely on it without future breaking
    changes.
"""

import logging
import sys
from dataclasses import dataclass
from pathlib import Path

from .exceptions import ConfigurationError


LOGGING_MANAGER_VERSION = "1.0"

DEFAULT_LOG_LEVEL: int = logging.INFO
DEFAULT_LOG_FORMAT: str = (
    "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
)
DEFAULT_DATE_FORMAT: str = "%Y-%m-%d %H:%M:%S"

_VALID_LOG_LEVELS: frozenset[int] = frozenset(
    {
        logging.DEBUG,
        logging.INFO,
        logging.WARNING,
        logging.ERROR,
        logging.CRITICAL,
    }
)

# Attribute used to mark handlers this module has added to the root
# logger, so configure_logging() can find and remove exactly those
# handlers (and no others) before adding fresh ones on a later call.

_MANAGED_HANDLER_ATTRIBUTE: str = "_blc_mark_managed"


__all__ = [
    "LOGGING_MANAGER_VERSION",
    "DEFAULT_LOG_LEVEL",
    "DEFAULT_LOG_FORMAT",
    "DEFAULT_DATE_FORMAT",
    "LoggingConfig",
    "configure_logging",
    "get_logger",
]


@dataclass(frozen=True)
class LoggingConfig:
    """Immutable, validated description of how logging should be set up.

    Attributes:
        log_level: Minimum severity level the root logger will emit.
            Must be one of logging.DEBUG, logging.INFO,
            logging.WARNING, logging.ERROR, or logging.CRITICAL.

        log_to_console: If True, log records are written to standard
            error via a console handler.

        log_to_file: If True, log records are also written to
            `log_file` via a file handler. Requires `log_file` to be
            set.

        log_file: Path to the log file to write to, if `log_to_file`
            is True. Ignored if `log_to_file` is False.

        log_format: logging.Formatter format string applied to every
            handler this module configures.

        date_format: logging.Formatter date format string applied to
            every handler this module configures.

    Raises:
        ConfigurationError: If `log_level` is not a recognized
            logging level, if `log_to_file` is True but `log_file` is
            not provided or points at an existing directory, or if
            `log_format` or `date_format` is empty.
        TypeError: If a field has an invalid type.
    """

    log_level: int = DEFAULT_LOG_LEVEL
    log_to_console: bool = True
    log_to_file: bool = False
    log_file: Path | None = None
    log_format: str = DEFAULT_LOG_FORMAT
    date_format: str = DEFAULT_DATE_FORMAT

    def __post_init__(self) -> None:
        """Validate every field immediately after construction."""

        if (
            isinstance(self.log_level, bool)
            or not isinstance(self.log_level, int)
        ):
            raise TypeError(
                "'log_level' must be an int, "
                f"got {type(self.log_level).__name__}."
            )

        if self.log_level not in _VALID_LOG_LEVELS:
            raise ConfigurationError(
                f"'log_level' must be one of {sorted(_VALID_LOG_LEVELS)}, "
                f"got {self.log_level!r}."
            )

        if not isinstance(self.log_to_console, bool):
            raise TypeError(
                "'log_to_console' must be a bool, "
                f"got {type(self.log_to_console).__name__}."
            )

        if not isinstance(self.log_to_file, bool):
            raise TypeError(
                "'log_to_file' must be a bool, "
                f"got {type(self.log_to_file).__name__}."
            )

        if self.log_file is not None and not isinstance(self.log_file, Path):
            raise TypeError(
                "'log_file' must be a pathlib.Path or None, "
                f"got {type(self.log_file).__name__}."
            )

        if not isinstance(self.log_format, str):
            raise TypeError(
                "'log_format' must be a str, "
                f"got {type(self.log_format).__name__}."
            )

        if not isinstance(self.date_format, str):
            raise TypeError(
                "'date_format' must be a str, "
                f"got {type(self.date_format).__name__}."
            )

        if not self.log_format.strip():
            raise ConfigurationError("'log_format' must not be empty.")

        if not self.date_format.strip():
            raise ConfigurationError("'date_format' must not be empty.")

        if self.log_to_file:
            if self.log_file is None:
                raise ConfigurationError(
                    "'log_file' must be provided when 'log_to_file' is True."
                )

            if self.log_file.exists() and self.log_file.is_dir():
                raise ConfigurationError(
                    "'log_file' must be a file path, not a directory: "
                    f"{self.log_file}"
                )


def configure_logging(config: LoggingConfig) -> None:
    """Configure the root logger according to a LoggingConfig.

    Uses explicit handler configuration rather than
    logging.basicConfig(). Calling this function more than once
    is safe: any handlers previously added by this function are
    removed and replaced with fresh ones reflecting the current
    `config`. This prevents duplicate BLC Mark handlers.

    Handlers not added by this module are left untouched.

    If `config.log_to_file` is True, the parent directory of
    `config.log_file` is created if it does not already exist,
    since a FileHandler cannot otherwise be opened.

    Args:
        config: The logging configuration to apply.

    Returns:
        None.

    Raises:
        TypeError: If `config` is not a LoggingConfig instance.
        OSError: If `config.log_to_file` is True and the log file's
            parent directory cannot be created, or the log file
            cannot be opened for writing.
    """
    if not isinstance(config, LoggingConfig):
        raise TypeError(
            "'config' must be a LoggingConfig, "
            f"got {type(config).__name__}."
        )

    root_logger = logging.getLogger()

    _remove_managed_handlers(root_logger)

    root_logger.setLevel(config.log_level)

    formatter = logging.Formatter(
        fmt=config.log_format,
        datefmt=config.date_format,
    )

    if config.log_to_console:
        console_handler = logging.StreamHandler(
            stream=sys.stderr
        )
        console_handler.setFormatter(formatter)

        setattr(
            console_handler,
            _MANAGED_HANDLER_ATTRIBUTE,
            True,
        )

        root_logger.addHandler(console_handler)

    if config.log_to_file:
        # LoggingConfig validation guarantees this is not None.
        assert config.log_file is not None

        config.log_file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        file_handler = logging.FileHandler(
            config.log_file,
            mode="a",
        )

        file_handler.setFormatter(formatter)

        setattr(
            file_handler,
            _MANAGED_HANDLER_ATTRIBUTE,
            True,
        )

        root_logger.addHandler(file_handler)


def _remove_managed_handlers(
    root_logger: logging.Logger,
) -> None:
    """Remove and close every handler previously added by
    configure_logging(), leaving unrelated handlers untouched.

    Args:
        root_logger: The root logger to clean up.

    Returns:
        None.
    """
    for handler in list(root_logger.handlers):
        if getattr(
            handler,
            _MANAGED_HANDLER_ATTRIBUTE,
            False,
        ):
            root_logger.removeHandler(handler)
            handler.close()


def get_logger(name: str) -> logging.Logger:
    """Obtain a named logger.

    This is a thin wrapper around logging.getLogger().

    Args:
        name: The logger name, conventionally `__name__` of the
            calling module.

    Returns:
        The logger for `name`, as returned by logging.getLogger().

    Raises:
        TypeError: If `name` is not a string.
        ValueError: If `name` is empty or contains only whitespace.
    """
    if not isinstance(name, str):
        raise TypeError(
            "'name' must be a str, "
            f"got {type(name).__name__}."
        )

    if not name.strip():
        raise ValueError("'name' must not be empty.")

    return logging.getLogger(name)


if __name__ == "__main__":
    import tempfile

    print("BLC Mark Logging Manager")
    print(f"Version: {LOGGING_MANAGER_VERSION}\n")

    console_config = LoggingConfig(
        log_level=logging.INFO,
        log_to_console=True,
    )

    configure_logging(console_config)

    demo_logger = get_logger(
        "blc_mark.logging_manager.demo"
    )

    demo_logger.info(
        "This is a demonstration INFO message."
    )

    demo_logger.warning(
        "This is a demonstration WARNING message."
    )

    with tempfile.TemporaryDirectory() as scratch_dir:
        demo_log_file = Path(scratch_dir) / "logs" / "demo.log"

        file_config = LoggingConfig(
            log_level=logging.INFO,
            log_to_console=True,
            log_to_file=True,
            log_file=demo_log_file,
        )

        configure_logging(file_config)

        demo_logger.info(
            "This message should appear on console and "
            "in the log file."
        )

        print(
            f"Log file created: {demo_log_file.exists()}"
        )

        print(
            f"Log file contents:\n"
            f"{demo_log_file.read_text(encoding='utf-8')}"
        )

        # Explicitly remove all handlers managed by this module.
        # This closes the FileHandler before TemporaryDirectory
        # attempts to remove the temporary directory on Windows.
        _remove_managed_handlers(logging.getLogger())

        # Ensure Python's logging system has completed any pending
        # handler shutdown work before temporary-file cleanup.
        logging.shutdown()

    print(
        "Logging configured and verified successfully."
    )




