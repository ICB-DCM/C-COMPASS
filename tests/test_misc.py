"""Miscellaneous tests for the ccompass package."""

import numpy as np
import pandas as pd

from ccompass.CCMPS import create_markerlist


def test_create_markerlist():
    marker_sets = {
        "somefile": {
            "class_col": "MarkerCompartment",
            "classes": ["PROTEIN - COMPLEX", "CYTOPLASM", "LYSOSOME"],
            "identifier_col": "Genename",
            "table": pd.DataFrame(
                {
                    "Genename": ["AAGAB", "AAK1", "AARS1"],
                    "MarkerCompartment": [
                        "CYTOPLASM",
                        "PROTEIN - COMPLEX",
                        "CYTOPLASM",
                    ],
                }
            ),
        }
    }
    marker_conv = {
        "PROTEIN - COMPLEX": "PROTEIN_COMPLEX",
        "CYTOPLASM": "CYTOPLASM",
        "LYSOSOME": "LYSOSOME",
        "ignored...": np.nan,
    }
    marker_params = {"how": "exclude", "what": "unite"}
    markerlist = create_markerlist(marker_sets, marker_conv, marker_params)
    assert markerlist.to_dict() == {
        "class": {
            "AAGAB": "CYTOPLASM",
            "AAK1": "PROTEIN_COMPLEX",
            "AARS1": "CYTOPLASM",
        }
    }
