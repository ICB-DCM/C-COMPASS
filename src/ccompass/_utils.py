"""Internal utility functions."""

import os
import platform
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
