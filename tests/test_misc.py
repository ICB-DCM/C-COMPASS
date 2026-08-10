"""Miscellaneous tests for the ccompass package."""

import numpy as np
import pandas as pd
from pandas._testing import assert_frame_equal

from ccompass._testing.synthetic_data import (
    SyntheticDataConfig,
    create_profiles,
    total_proteome,
)
from ccompass.core import (
    ConditionPredictionModel,
    MarkerSet,
    TrainingRoundModel,
    create_markerlist,
)
from ccompass.MOA import stats_proteome


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


def test_synth_data_deterministic():
    """Test that the synthetic data generation is deterministic."""
    c = SyntheticDataConfig()
    fractionation_df1, marker_df1 = create_profiles(c=c)
    total_prot_df1 = total_proteome(
        proteins=list(fractionation_df1[c.protein_id_col]), c=c
    )

    c = SyntheticDataConfig()
    fractionation_df2, marker_df2 = create_profiles(c=c)
    total_prot_df2 = total_proteome(
        proteins=list(fractionation_df2[c.protein_id_col]), c=c
    )

    assert_frame_equal(fractionation_df1, fractionation_df2)
    assert_frame_equal(marker_df1, marker_df2)
    assert_frame_equal(total_prot_df1, total_prot_df2)


def test_stats_proteome_all_na_row():
    """`stats_proteome` must not crash if a protein has no valid class
    prediction in any replicate/round of a condition (i.e. its `CC_*`
    columns are all NaN).

    Previously raised ``ValueError: Encountered all NA values`` from
    `DataFrame.idxmax` (mirrors the `fCC`/`fNN_winner` fix from #143, which
    was never applied to the analogous `CC`/`NN_winner` computation).
    """
    # Two replicates of "Cond1", both agreeing that "P1" is class "A".
    # "P2" has no valid prediction in either replicate (e.g. it was filtered
    # out before NN training due to an incomplete profile), so it is only
    # known via `fract_data`/`fract_marker`, not via any
    # `z_full_df`/`w_full_prob_df`.
    subcons = ["Cond1_Rep.1", "Cond1_Rep.2"]
    class_predictions = {
        subcon: ConditionPredictionModel(
            condition_id=subcon,
            classes=["A", "B"],
            round_results={
                "round1": TrainingRoundModel(
                    w_full_prob_df=pd.DataFrame(
                        {"SVM_winner": ["A"], "SVM_prob": [0.9]},
                        index=["P1"],
                    ),
                    z_full_df=pd.DataFrame(
                        {"A": [0.7], "B": [0.3]}, index=["P1"]
                    ),
                )
            },
        )
        for subcon in subcons
    }
    fract_data = {
        "class": {
            subcon: pd.DataFrame(index=["P1", "P2"]) for subcon in subcons
        }
    }
    fract_marker = {
        subcon: pd.DataFrame({"class": ["A", "B"]}, index=["P1", "P2"])
        for subcon in subcons
    }

    results = stats_proteome(
        class_predictions,
        fract_data,
        fract_marker,
        fract_conditions=["Cond1"],
        reliability=50,
    )

    metrics = results["Cond1"].metrics
    assert metrics.loc["P1", "NN_winner"] == "A"
    assert metrics.loc["P2", "NN_winner"] == "nan"
    assert metrics.loc["P1", "fNN_winner"] == "A"
    assert metrics.loc["P2", "fNN_winner"] == "nan"
