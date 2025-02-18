import os
import re
from pathlib import Path

from create_synthetic_data import (
    class_id_col,
    create_profiles,
    gene_id_col,
    protein_id_col,
    total_proteome,
)

from ccompass.core import (
    IDENTIFIER,
    KEEP,
    NA,
    MarkerSet,
    NeuralNetworkParametersModel,
    SessionModel,
)
from ccompass.FDP import start_fract_data_processing
from ccompass.main_gui import (
    create_fullprofiles,
    create_identity_conversion,
    create_marker_profiles,
    create_markerlist,
    logger,
)
from ccompass.MOA import class_comparisons, global_comparisons, stats_proteome
from ccompass.TPP import start_total_proteome_processing

# regexes to parse column IDs
fract_id_rx = re.compile(
    r"(?P<condition>Con\d+)_Rep(?P<replicate>\d+)_Fr(?P<fraction>\d+)"
)
tp_id_rx = re.compile(r"(?P<condition>Con\d+)_Rep(?P<replicate>\d+)")


def fract_col_id_to_row(col_id: str) -> list:
    """Convert fractionation data column id to fractionation table rows."""
    if col_id == protein_id_col:
        return [col_id, IDENTIFIER, NA, NA]
    if col_id == gene_id_col:
        return [col_id, KEEP, NA, NA]

    if not (match := fract_id_rx.match(col_id)):
        raise ValueError(f"Invalid fractionation ID: {col_id}")

    condition = match["condition"]
    replicate = int(match["replicate"])
    fraction = int(match["fraction"])
    return [col_id, condition, replicate, fraction]


def tp_col_id_to_row(col_id: str) -> list:
    """Convert total proteome data column id to total proteome table rows."""
    if col_id == protein_id_col:
        return [col_id, IDENTIFIER]

    if not (match := tp_id_rx.match(col_id)):
        raise ValueError(f"Invalid total proteome ID: {col_id}")

    condition = match["condition"]
    return [col_id, condition]


def test_full():
    """Check that we can run the full analysis.

    For now, no ground truth, only check that the code runs without crashing.
    """
    max_procs = os.cpu_count()

    # generate synthetic data
    fractionation_df0, marker_df = create_profiles()
    total_prot_df = total_proteome(
        proteins=list(fractionation_df0[protein_id_col])
    )
    fractionation_df = fractionation_df0.drop(columns=[class_id_col])
    # uppercase is expected elsewhere
    marker_df = marker_df.apply(lambda x: x.astype(str).str.upper())

    # simulate user input
    sess = SessionModel()
    fract_filepath = "bla/fract.csv"
    marker_filepath = "bla/marker.csv"
    total_prot_filepath = "bla/total_prot.csv"
    sess.fract_tables = {
        fract_filepath: [
            fract_col_id_to_row(col_id)
            for col_id in fractionation_df
            if not col_id.startswith("Amount_")
        ]
    }
    sess.fract_indata = {fract_filepath: fractionation_df}
    sess.fract_identifiers = {fract_filepath: protein_id_col}

    # process fractionation data
    (
        sess.fract_data,
        sess.fract_std,
        sess.fract_info,
        sess.fract_conditions,
    ) = start_fract_data_processing(
        sess.fract_tables,
        sess.fract_preparams,
        sess.fract_identifiers,
        sess.fract_indata,
    )

    # process marker data
    sess.marker_sets = {
        marker_filepath: MarkerSet(
            df=marker_df,
            identifier_col=gene_id_col,
            class_col=class_id_col,
        )
    }
    sess.marker_fractkey = gene_id_col
    sess.marker_conv = create_identity_conversion(sess.marker_sets)

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
    sess.tp_tables = {
        total_prot_filepath: [
            tp_col_id_to_row(col_id)
            for col_id in total_prot_df
            if not col_id.startswith("RelativeRegulation")
        ]
    }
    sess.tp_indata = {total_prot_filepath: total_prot_df}
    sess.tp_identifiers = {total_prot_filepath: protein_id_col}

    sess.tp_data, sess.tp_info, sess.tp_conditions, sess.tp_icorr = (
        start_total_proteome_processing(
            sess.tp_data,
            sess.tp_tables,
            sess.tp_preparams,
            sess.tp_identifiers,
            sess.tp_info,
            sess.tp_icorr,
            sess.tp_indata,
            sess.tp_conditions,
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
        sess.fract_std,
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

    sess.status.fractionation_data = True
    sess.status.tp_data = True
    sess.status.marker_file = True
    sess.status.marker_matched = True
    sess.status.training = True
    sess.status.proteome_prediction = True
    sess.status.comparison_global = True
    sess.status.comparison_class = True

    # TODO test all generating all reports
    # TODO add some checks
    # check that we have results for all conditions
    conditions = {
        row[1]
        for row in next(iter(sess.fract_tables.values()))
        if row[1] not in (IDENTIFIER, KEEP)
    }
    assert set(sess.results.keys()) == conditions

    ...
    sess.to_numpy(Path(__file__).parent / "session_test_full.npy")
