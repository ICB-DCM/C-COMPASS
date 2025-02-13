"""C-COMPASS


Acronyms:

* CA: class abundance
* CC: class contribution (contribution of a compartment to a profile,
  CC ∈ [0, 1])
* DS: distance score
* RL: relocalization (difference between two class contributions, RL ∈ [-1, 1])
* RLS: relocalization score (sum of RL values across all compartments)
  RLS ∈ [0, 2] (no relocalization .. full relocalization)
* nCC: normalized class contribution
"""

from pathlib import Path

from ._utils import get_ccmps_data_directory

__all__ = []

# the application settings file
config_filepath: Path = get_ccmps_data_directory() / "settings.yaml"

# the repository URL
repository_url = "https://github.com/ICB-DCM/C-COMPASS/"
# the ReadTheDocs URL
readthedocs_url = "https://c-compass.readthedocs.io/en/latest/"

# name of the application
app_name = "C-COMPASS"
