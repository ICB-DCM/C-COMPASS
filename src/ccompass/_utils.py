"""Internal utility functions."""

import logging
import os
import platform
from collections.abc import Iterable
from pathlib import Path


def get_data_directory() -> Path:
    """Get the platform-specific data directory."""
    system = platform.system()

    if system == "Linux":
        return Path(
            os.getenv(
                "XDG_DATA_HOME",
                os.path.join(os.path.expanduser("~"), ".local", "share"),
            )
        )

    if system == "Windows":
        return Path(
            os.getenv(
                "APPDATA",
                os.path.join(os.path.expanduser("~"), "AppData", "Roaming"),
            )
        )

    if system == "Darwin":
        return Path(os.path.expanduser("~"), "Library", "Application Support")

    raise NotImplementedError(f"Unsupported platform: {system}")


def get_ccmps_data_directory() -> Path:
    """Get the platform-specific data directory for ccmps."""
    return get_data_directory() / "C-COMPASS"


def unique_preserve_order(seq: Iterable) -> list:
    """Return a list of unique elements from some iterable, keeping only the
    first occurrence of each element.

    :param seq: Sequence to prune
    :return: List of unique elements in ``seq``
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


class PrefixFilter(logging.Filter):
    """A logging filter that adds a prefix to the log message."""

    def __init__(self, prefix):
        super().__init__()
        self.prefix = prefix

    def filter(self, record):
        record.msg = f"{self.prefix} {record.msg}"
        return True
