"""C-COMPASS"""

import os
from pathlib import Path

from ._utils import get_ccmps_data_directory

__all__ = []

# tensorflow logging
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

# the application settings file
config_filepath: Path = get_ccmps_data_directory() / "settings.yaml"
