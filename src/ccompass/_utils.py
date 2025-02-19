"""Internal utility functions."""

import logging
import multiprocessing
import os
import platform
import sys
from collections.abc import Iterable
from contextlib import contextmanager
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


class StreamToLogger:
    """Redirect a stream (e.g., stdout) to a logger."""

    def __init__(self, logger: logging.Logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ""

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass


@contextmanager
def stdout_to_logger(logger: logging.Logger, log_level=logging.INFO):
    """Context manager to redirect stdout to a logger."""
    old_stdout = sys.stdout
    try:
        sys.stdout = StreamToLogger(logger, log_level)
        yield
    finally:
        sys.stdout = old_stdout


def get_mp_ctx() -> multiprocessing.context.BaseContext:
    """Get the multiprocessing context."""
    if platform.system() == "Windows":
        return multiprocessing.get_context("spawn")
    else:
        return multiprocessing.get_context("forkserver")
