import os
from pathlib import Path

import numpy as np
from test_session_model import assert_session_equal

from ccompass._testing.synthetic_data import (
    SyntheticDataConfig,
    create_profiles,
    fract_col_id_to_row,
    total_proteome,
    tp_col_id_to_row,
)
from ccompass.core import (
    IDENTIFIER,
    KEEP,
    FractDataset,
    MarkerSet,
    NeuralNetworkParametersModel,
    SessionModel,
    TotalProtDataset,
    create_fullprofiles,
    create_identity_conversion,
    create_marker_profiles,
    create_markerlist,
)
from ccompass.FDP import start_fract_data_processing
from ccompass.main_gui import (
    logger,
)
from ccompass.MOA import class_comparisons, global_comparisons, stats_proteome
from ccompass.TPP import start_total_proteome_processing


def test_full():
    """Check that we can run the full analysis.

    For now, no ground truth, only check that the code runs without crashing.
    """
    max_procs = os.cpu_count()

    # generate synthetic data
    c = SyntheticDataConfig()
    fractionation_df0, marker_df = create_profiles(c=c)
    total_prot_df = total_proteome(
        proteins=list(fractionation_df0[c.protein_id_col]), c=c
    )
    fractionation_df = fractionation_df0.drop(columns=[c.class_id_col])
    # uppercase is expected elsewhere
    marker_df = marker_df.apply(lambda x: x.astype(str).str.upper())

    # simulate user input
    fract_filepath = "bla/fract.csv"
    marker_filepath = "bla/marker.csv"
    total_prot_filepath = "bla/total_prot.csv"
    fract_dset = FractDataset(
        df=fractionation_df,
        table=[
            fract_col_id_to_row(col_id, c)
            for col_id in fractionation_df
            if not col_id.startswith("Amount_")
        ],
    )
    tp_dset = TotalProtDataset(
        df=total_prot_df,
        table=[
            tp_col_id_to_row(col_id, c=c)
            for col_id in total_prot_df
            if not col_id.startswith("RelativeRegulation")
        ],
    )
    sess = SessionModel(
        fract_input={fract_filepath: fract_dset},
    )

    # process fractionation data
    (
        sess.fract_data,
        sess.fract_std,
        sess.fract_info,
        sess.fract_conditions,
    ) = start_fract_data_processing(
        sess.fract_input,
        sess.fract_preparams,
    )

    # process marker data
    sess.marker_sets = {
        marker_filepath: MarkerSet(
            df=marker_df,
            identifier_col=c.gene_id_col,
            class_col=c.class_id_col,
        )
    }
    sess.marker_fractkey = c.gene_id_col
    sess.marker_conv = create_identity_conversion(sess.marker_sets.values())

    sess.marker_list = create_markerlist(
        sess.marker_sets,
        sess.marker_conv,
        **sess.marker_params,
    )

    logger.info("Marker list created")
    (
        sess.fract_marker,
        sess.fract_marker_vis,
        sess.fract_test,
    ) = create_marker_profiles(
        sess.fract_data,
        sess.marker_fractkey,
        sess.fract_info,
        sess.marker_list,
    )
    logger.info("Marker profiles created")
    sess.fract_full = create_fullprofiles(sess.fract_marker, sess.fract_test)
    logger.info("Full profiles created")

    # process total proteome data
    sess.tp_input = {total_prot_filepath: tp_dset}

    sess.tp_data, sess.tp_info, sess.tp_icorr = (
        start_total_proteome_processing(
            sess.tp_input,
            sess.tp_preparams,
            sess.tp_data,
            sess.tp_info,
            sess.tp_icorr,
        )
    )

    # train model
    from ccompass.MOP import multi_organelle_prediction

    sess.NN_params = NeuralNetworkParametersModel(
        rounds=1,
        subrounds=3,
        optimizers=["adam"],
        NN_epochs=10,
        NN_optimization="short",
    )
    sess.learning_xyz = multi_organelle_prediction(
        sess.fract_full,
        sess.fract_marker,
        sess.fract_test,
        sess.NN_params,
        max_procs,
    )

    # "static statistics"
    sess.results = stats_proteome(
        sess.learning_xyz,
        sess.fract_data,
        sess.fract_conditions,
        sess.NN_params.reliability,
    )
    for condition, result in sess.results.items():
        ccols = [f"CC_{cls}" for cls in result.classnames]
        assert np.isclose(result.metrics.loc[:, ccols].sum(axis=1), 1).all()

        for cls in result.classnames:
            assert (result.metrics[f"CC_{cls}"] >= 0).all()
            assert (result.metrics[f"CC_{cls}"] <= 1).all()
            assert (result.metrics[f"fCC_{cls}"] >= 0).all()
            assert (result.metrics[f"fCC_{cls}"] <= 1).all()
            assert (
                (result.metrics[f"fCC_{cls}"] >= result.metrics[f"CC_{cls}"])
                | (result.metrics[f"fCC_{cls}"] == 0.0)
            ).all(), result.metrics[[f"fCC_{cls}", f"CC_{cls}"]][
                ~(
                    (
                        result.metrics[f"fCC_{cls}"]
                        >= result.metrics[f"CC_{cls}"]
                    )
                    | (result.metrics[f"fCC_{cls}"] == 0.0)
                )
            ]

    # "global changes"
    sess.comparison = global_comparisons(
        sess.results,
        max_procs,
    )

    # "class-centric changes"
    class_comparisons(
        sess.tp_data,
        sess.results,
        sess.comparison,
    )

    # TODO test all generating all reports
    # TODO add some checks
    # check that we have results for all conditions
    conditions = {
        row[1]
        for row in next(iter(sess.fract_input.values())).table
        if row[1] not in (IDENTIFIER, KEEP)
    }
    assert set(sess.results.keys()) == conditions

    ...
    sess.to_numpy(Path(__file__).parent / "session_test_full.npy")
    sess.to_zip(Path(__file__).parent / "session_test_full.ccompass")

    sess2 = SessionModel.from_numpy(
        Path(__file__).parent / "session_test_full.npy"
    )
    assert_session_equal(sess, sess2)
    sess2 = SessionModel.from_zip(
        Path(__file__).parent / "session_test_full.ccompass"
    )
    assert_session_equal(sess, sess2)
