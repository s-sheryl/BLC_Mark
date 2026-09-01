"""
Purpose:
    Provide small, generic filesystem utility functions reused across
    BLC Mark -- checking existence, ensuring directories, safe
    deletion, size lookups, and read/write access checks.

Responsibilities:
    - Wrap common pathlib/os filesystem operations in a consistent,
      well-documented, type-validated interface.
    - Raise standard Python exceptions (FileNotFoundError,
      IsADirectoryError, NotADirectoryError, PermissionError, OSError)
      with clear messages rather than failing silently or returning
      ambiguous sentinel values.

Scope:
    This module contains only generic filesystem utilities. It does
    not read CSV/TSV files, does not parse JSON, does not download
    files, does not configure logging, does not validate datasets,
    does not preprocess data, does not communicate with Xena, and
    does not perform any biological analysis. Every function here
    would be equally at home in a project that has nothing to do with
    cancer genomics -- no function in this module knows anything
    about BLC Mark's domain.

Version:
    This is a frozen Version 1.0 filesystem utility API. Its function
    names and signatures are intended to remain stable so dependent
    modules (configuration_manager.py, logging_manager.py,
    download_manager.py, dataset_manager.py, metadata_manager.py,
    preprocessing_manager.py, quality_control.py) can rely on them
    without future breaking changes.
"""

import os
import shutil
from pathlib import Path


FILE_UTILS_VERSION = "1.0"


__all__ = [
    "FILE_UTILS_VERSION",
    "ensure_directory_exists",
    "file_exists",
    "directory_exists",
    "delete_file",
    "delete_directory",
    "get_file_size",
    "is_directory_empty",
    "resolve_path",
    "is_readable",
    "is_writable",
]


def _require_path(path: Path, argument_name: str = "path") -> None:
    """Confirm an argument is a pathlib.Path instance.

    Args:
        path: The value to check.
        argument_name: Name of the argument, used only in the error
            message.

    Raises:
        TypeError: If `path` is not a pathlib.Path instance.
    """
    if not isinstance(path, Path):
        raise TypeError(
            f"'{argument_name}' must be a pathlib.Path, "
            f"got {type(path).__name__}."
        )


def _require_bool(value: bool, argument_name: str) -> None:
    """Confirm an argument is a bool instance.

    Args:
        value: The value to check.
        argument_name: Name of the argument, used only in the error
            message.

    Raises:
        TypeError: If `value` is not a bool instance.
    """
    if not isinstance(value, bool):
        raise TypeError(
            f"'{argument_name}' must be a bool, "
            f"got {type(value).__name__}."
        )


def ensure_directory_exists(path: Path) -> None:
    """Create a directory, including any missing parent directories,
    if it does not already exist.

    Args:
        path: The directory path to ensure exists.

    Returns:
        None.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
        NotADirectoryError: If `path` already exists but is a file
            rather than a directory.
        PermissionError: If the directory cannot be created due to
            insufficient permissions.
        OSError: If the directory cannot be created for any other
            filesystem reason.
    """
    _require_path(path)

    if path.exists() and not path.is_dir():
        raise NotADirectoryError(
            f"Path exists but is not a directory: {path}"
        )

    path.mkdir(parents=True, exist_ok=True)


def file_exists(path: Path) -> bool:
    """Check whether a regular file exists at the given path.

    Args:
        path: The path to check.

    Returns:
        True if `path` exists and is a regular file. False if `path`
        does not exist or exists but is a directory.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
    """
    _require_path(path)
    return path.is_file()


def directory_exists(path: Path) -> bool:
    """Check whether a directory exists at the given path.

    Args:
        path: The path to check.

    Returns:
        True if `path` exists and is a directory. False if `path`
        does not exist or exists but is a file.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
    """
    _require_path(path)
    return path.is_dir()


def delete_file(path: Path, *, missing_ok: bool = True) -> None:
    """Delete a file safely.

    Args:
        path: The file to delete.
        missing_ok: If True (default), deleting a path that does not
            exist is a no-op. If False, a missing file raises
            FileNotFoundError.

    Returns:
        None.

    Raises:
        TypeError: If `path` is not a pathlib.Path or `missing_ok`
            is not a bool.
        FileNotFoundError: If `path` does not exist and `missing_ok`
            is False.
        IsADirectoryError: If `path` exists but is a directory, not a
            file.
        PermissionError: If the file cannot be deleted due to
            insufficient permissions.
        OSError: If the file cannot be deleted for any other
            filesystem reason.
    """
    _require_path(path)
    _require_bool(missing_ok, "missing_ok")

    if not path.exists():
        if missing_ok:
            return
        raise FileNotFoundError(
            f"Cannot delete; file does not exist: {path}"
        )

    if path.is_dir():
        raise IsADirectoryError(
            f"Cannot delete file; path is a directory: {path}"
        )

    path.unlink()


def delete_directory(path: Path, *, missing_ok: bool = True) -> None:
    """Delete a directory and everything inside it, safely.

    Args:
        path: The directory to delete.
        missing_ok: If True (default), deleting a path that does not
            exist is a no-op. If False, a missing directory raises
            FileNotFoundError.

    Returns:
        None.

    Raises:
        TypeError: If `path` is not a pathlib.Path or `missing_ok`
            is not a bool.
        FileNotFoundError: If `path` does not exist and `missing_ok`
            is False.
        NotADirectoryError: If `path` exists but is a file, not a
            directory.
        PermissionError: If the directory cannot be deleted due to
            insufficient permissions.
        OSError: If the directory cannot be deleted for any other
            filesystem reason.
    """
    _require_path(path)
    _require_bool(missing_ok, "missing_ok")

    if not path.exists():
        if missing_ok:
            return
        raise FileNotFoundError(
            f"Cannot delete; directory does not exist: {path}"
        )

    if not path.is_dir():
        raise NotADirectoryError(
            f"Cannot delete directory; path is a file: {path}"
        )

    shutil.rmtree(path)


def get_file_size(path: Path) -> int:
    """Get the size of a file in bytes.

    Args:
        path: The file to measure.

    Returns:
        The file size in bytes.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
        FileNotFoundError: If `path` does not exist.
        IsADirectoryError: If `path` exists but is a directory, not a
            file.
    """
    _require_path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Cannot get size; file does not exist: {path}"
        )

    if path.is_dir():
        raise IsADirectoryError(
            f"Cannot get size; path is a directory: {path}"
        )

    return path.stat().st_size


def is_directory_empty(path: Path) -> bool:
    """Check whether a directory contains no entries.

    Args:
        path: The directory to check.

    Returns:
        True if `path` is a directory containing no files or
        subdirectories.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
        FileNotFoundError: If `path` does not exist.
        NotADirectoryError: If `path` exists but is a file, not a
            directory.
    """
    _require_path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Directory does not exist: {path}"
        )

    if not path.is_dir():
        raise NotADirectoryError(
            f"Path is not a directory: {path}"
        )

    return not any(path.iterdir())


def resolve_path(path: Path) -> Path:
    """Resolve a path to its absolute form, expanding `~` and
    resolving any `..` components and symlinks.

    Does not require `path` to exist -- resolution is purely
    syntactic and symlink-aware, not a validity check.

    Args:
        path: The path to resolve.

    Returns:
        The absolute, resolved path.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
    """
    _require_path(path)
    return path.expanduser().resolve()


def is_readable(path: Path) -> bool:
    """Check whether a path can be read.

    Args:
        path: The path to check.

    Returns:
        True if `path` exists and the current process has read
        permission on it. False if `path` does not exist or is not
        readable.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
    """
    _require_path(path)

    if not path.exists():
        return False

    return os.access(path, os.R_OK)


def is_writable(path: Path) -> bool:
    """Check whether a path can be written to.

    If `path` already exists, this checks write permission on `path`
    itself. If `path` does not exist, this checks write permission on
    its parent directory instead, since the meaningful question for a
    non-existent path is whether it could be created there.

    Args:
        path: The path to check.

    Returns:
        True if `path` (or its parent directory, if `path` does not
        exist) is writable by the current process. False otherwise,
        including when the parent directory itself does not exist.

    Raises:
        TypeError: If `path` is not a pathlib.Path.
    """
    _require_path(path)

    if path.exists():
        return os.access(path, os.W_OK)

    if not path.parent.exists():
        return False

    return os.access(path.parent, os.W_OK)


if __name__ == "__main__":
    import tempfile

    print("BLC Mark File Utilities")
    print(f"Version: {FILE_UTILS_VERSION}\n")

    with tempfile.TemporaryDirectory() as scratch_dir:
        root = Path(scratch_dir)

        nested_dir = root / "nested" / "subdirectory"
        ensure_directory_exists(nested_dir)
        print(
            "ensure_directory_exists: "
            f"created {nested_dir} -> {directory_exists(nested_dir)}"
        )

        sample_file = nested_dir / "sample.txt"
        sample_file.write_text(
            "BLC Mark file_utils demonstration.\n",
            encoding="utf-8",
        )

        print(f"file_exists: {file_exists(sample_file)}")
        print(
            "directory_exists (on a file): "
            f"{directory_exists(sample_file)}"
        )

        print(f"get_file_size: {get_file_size(sample_file)} bytes")
        print(
            "is_directory_empty (nested_dir): "
            f"{is_directory_empty(nested_dir)}"
        )

        empty_dir = root / "empty"
        ensure_directory_exists(empty_dir)
        print(
            "is_directory_empty (empty_dir): "
            f"{is_directory_empty(empty_dir)}"
        )

        print(f"resolve_path: {resolve_path(Path('.'))}")
        print(f"is_readable: {is_readable(sample_file)}")
        print(
            "is_writable (existing file): "
            f"{is_writable(sample_file)}"
        )
        print(
            "is_writable (new path in existing dir): "
            f"{is_writable(nested_dir / 'new_file.txt')}"
        )

        delete_file(sample_file)
        print(
            "delete_file: file_exists after delete -> "
            f"{file_exists(sample_file)}"
        )

        delete_directory(root / "nested")
        print(
            "delete_directory: directory_exists after delete -> "
            f"{directory_exists(root / 'nested')}"
        )

        delete_file(
            root / "does_not_exist.txt",
            missing_ok=True,
        )
        print(
            "delete_file with missing_ok=True on absent file: "
            "no error raised"
        )

    print("\nAll file_utils functions demonstrated successfully.")

