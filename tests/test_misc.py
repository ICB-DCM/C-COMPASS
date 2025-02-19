"""Miscellaneous tests for the ccompass package."""

import numpy as np
import pandas as pd

from ccompass.core import MarkerSet, create_markerlist


def test_create_markerlist():
    marker_sets = {
        "somefile": MarkerSet(
            class_col="MarkerCompartment",
            identifier_col="Genename",
            df=pd.DataFrame(
                {
                    "Genename": [
                        "AAGAB",
                        "AAK1",
                        "AARS1",
                        "only_in_first",
                        "mismatch",
                    ],
                    "MarkerCompartment": [
                        "CYTOPLASM",
                        "PROTEIN - COMPLEX",
                        "CYTOPLASM",
                        "PROTEIN - COMPLEX",
                        "CYTOPLASM",
                    ],
                    "ignored...": [np.nan] * 5,
                }
            ),
        ),
        "somefile2": MarkerSet(
            class_col="MarkerCompartment",
            identifier_col="Genename",
            df=pd.DataFrame(
                {
                    "Genename": [
                        "AAGAB",
                        "AAK1",
                        "AARS1",
                        "only_in_second",
                        "mismatch",
                    ],
                    "MarkerCompartment": [
                        "CYTOPLASM",
                        "PROTEIN - COMPLEX",
                        "CYTOPLASM",
                        "PROTEIN - COMPLEX",
                        "PROTEIN - COMPLEX",
                    ],
                    "ignored...": [np.nan] * 5,
                }
            ),
        ),
    }
    marker_conv = {
        "PROTEIN - COMPLEX": "PROTEIN_COMPLEX",
        "CYTOPLASM": "CYTOPLASM",
        "LYSOSOME": "LYSOSOME",
        "ignored...": np.nan,
    }

    markerlist = create_markerlist(
        marker_sets, marker_conv, what="unite", how="exclude"
    )
    assert markerlist.to_dict() == {
        "class": {
            "AAGAB": "CYTOPLASM",
            "AAK1": "PROTEIN_COMPLEX",
            "AARS1": "CYTOPLASM",
            "only_in_first": "PROTEIN_COMPLEX",
            "only_in_second": "PROTEIN_COMPLEX",
        }
    }

    markerlist = create_markerlist(
        marker_sets, marker_conv, what="intersect", how="exclude"
    )
    assert markerlist.to_dict() == {
        "class": {
            "AAGAB": "CYTOPLASM",
            "AAK1": "PROTEIN_COMPLEX",
            "AARS1": "CYTOPLASM",
        }
    }

    # FIXME: what to do in case of a tie and "majority"?
    #  Currently, the first marker set wins.
    marker_sets["somefile3"] = marker_sets["somefile2"].model_copy(deep=True)
    markerlist = create_markerlist(
        marker_sets, marker_conv, what="intersect", how="majority"
    )
    assert markerlist.to_dict() == {
        "class": {
            "AAGAB": "CYTOPLASM",
            "AAK1": "PROTEIN_COMPLEX",
            "AARS1": "CYTOPLASM",
            "mismatch": "PROTEIN_COMPLEX",
        }
    }
