"""Fractionation data processing."""

import copy
import logging
from collections import Counter
from itertools import chain
from tkinter import messagebox
from typing import Literal

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler

from ._utils import unique_preserve_order
from .core import IDENTIFIER, KEEP, FractDataset

logger = logging.getLogger(__package__)


def create_dataset(
    fract_input: dict[str, FractDataset],
    conditions: list[str],
    window: sg.Window,
    progress: float,
) -> tuple[dict[str, dict[str, pd.DataFrame]], dict[str, pd.DataFrame], float]:
    """Consolidate fractionation data.

    :returns: A tuple containing the dataset, the data to keep and the updated progress.
    The dataset is a dictionary with the following structure:
    - condition_id -> replicate_id -> DataFrame
    - all datasets have the same index (protein identifiers)
    """
    conditions = [x for x in conditions if x not in ["", IDENTIFIER]]

    # collect all identifiers
    all_identifiers = list(
        set(
            chain.from_iterable(
                dset.df[dset.id_col] for dset in fract_input.values()
            )
        )
    )

    # progress increment per condition
    stepsize = 10.0 / len(conditions)
    # dataset to be created
    #  condition_id -> replicate_id -> DataFrame
    combined: dict[str, dict[str, pd.DataFrame]] = {}

    for condition in conditions:
        progress += stepsize
        if window:
            window["--status2--"].update(condition)
            window["--progress--"].update(progress)
            window.read(timeout=50)

        data_new = pd.DataFrame(index=all_identifiers)
        for dset in fract_input.values():
            for (
                samplename,
                sample_condition,
                sample_replicate,
                sample_fraction,
            ) in dset.table:
                if sample_condition != condition:
                    continue

                data = pd.DataFrame()
                data[samplename] = dset.df[samplename]
                data.set_index(dset.df[dset.id_col], inplace=True)
                data_new = pd.merge(
                    data_new,
                    data,
                    right_index=True,
                    left_index=True,
                    how="outer",
                )

                if condition == KEEP:
                    if samplename + "_x" in data_new.columns:
                        # handle merge conflicts
                        for element in list(data_new.index):
                            if pd.isnull(data_new[samplename + "_x"][element]):
                                data_new[samplename + "_x"][element] = (
                                    data_new[samplename + "_y"][element]
                                )
                            if pd.isnull(data_new[samplename + "_y"][element]):
                                data_new[samplename + "_y"][element] = (
                                    data_new[samplename + "_x"][element]
                                )
                    data_new = data_new.T.drop_duplicates().T
                    data_new.rename(
                        {samplename + "_x": samplename},
                        axis=1,
                        inplace=True,
                    )
                else:
                    data_new = data_new.rename(
                        columns={
                            samplename: f"Fr.{sample_fraction}_{samplename}_Rep.{sample_replicate}"
                        }
                    )
        # list of unique replicate numbers
        replist = []
        for sample in data_new.columns:
            # get replicate number
            rep = sample[sample.rfind("_") + 1 :]
            if rep not in replist:
                replist.append(rep)
        # collect all samples for each replicate
        repdata = {}
        for rep in replist:
            data = pd.DataFrame(index=data_new.index)
            for sample in data_new.columns:
                if sample[sample.rfind("_") + 1 :] == rep:
                    data = pd.merge(
                        data,
                        data_new[sample],
                        left_index=True,
                        right_index=True,
                    )
            repdata[rep] = data
        combined[condition] = repdata

    data_keep = {}
    if KEEP in combined:
        data_keep = combined[KEEP]
        del combined[KEEP]

    return combined, data_keep, progress


def pre_post_scaling(
    data: dict[str, dict[str, pd.DataFrame]],
    how: Literal["minmax", "area"],
    window: sg.Window,
    progress: int,
):
    """Scale data using MinMaxScaler or area normalization.

    Returns a tuple containing the scaled data and the updated progress.
    """
    for condition in data:
        stepsize = (5.0 / len(data)) / len(data[condition])
        for replicate in data[condition]:
            progress += stepsize
            if window:
                window["--progress--"].update(progress)
                window["--status2--"].update(" ".join([condition, replicate]))
                window.read(timeout=50)

            if how == "minmax":
                scaler = MinMaxScaler()
                data[condition][replicate] = pd.DataFrame(
                    scaler.fit_transform(data[condition][replicate].T).T,
                    columns=data[condition][replicate].columns,
                ).set_index(data[condition][replicate].index)
            elif how == "area":
                data[condition][replicate] = data[condition][replicate].div(
                    data[condition][replicate].sum(axis=1), axis=0
                )
            else:
                raise ValueError(f"Unknown scaling method: {how}")
    return data, progress


def filter_missing(
    data: dict[str, pd.DataFrame], mincount: int, window, progress
):
    """Remove rows with at least mincount non-zero values and set remaining
    N/A values to zero."""
    for condition in data:
        stepsize = (5.0 / len(data)) / len(data[condition])
        for replicate in data[condition]:
            progress += stepsize
            if window:
                window["--progress--"].update(progress)
                window["--status2--"].update(" ".join([condition, replicate]))
                window.read(timeout=50)

            data[condition][replicate].dropna(thresh=mincount, inplace=True)
            data[condition][replicate].replace(np.nan, 0.0, inplace=True)
    return data, progress


def filter_count(
    data: dict[str, dict[str, pd.DataFrame]],
    mincount: int,
    window: sg.Window,
    progress: float,
) -> tuple[dict[str, dict[str, pd.DataFrame]], dict[str, list[str]], float]:
    """Filter data by minimum number of replicates.

    Returns a tuple containing the filtered data, the remaining protein list
    and the updated progress.
    """
    protlist_remaining = {}
    for condition in data:
        stepsize = (5.0 / len(data)) / len(data[condition])
        peplist = []
        for replicate in data[condition]:
            peplist.extend(data[condition][replicate].index)
        peplist = list(set(remove_elements(peplist, mincount)))
        for replicate in data[condition]:
            progress += stepsize
            if window:
                window["--progress--"].update(progress)
                window["--status2--"].update(" ".join([condition, replicate]))
                window.read(timeout=50)

            # drop rows that are not in peplist
            for index in list(data[condition][replicate].index):
                if index not in peplist:
                    data[condition][replicate].drop(
                        index, axis=0, inplace=True
                    )
        protlist_remaining[condition] = peplist
    return data, protlist_remaining, progress


def remove_elements(lst: list, k: int) -> list:
    """Remove elements that occur less than k times in a list.

    Returns a list with elements that occur at least k times in the input list.
    """
    counted = Counter(lst)
    return [el for el in lst if counted[el] >= k]


def list_samples(data, window, progress):
    fracts_con = {}
    fracts_count = {}
    fracts_corr = {}
    for condition in data:
        fracts_count[condition] = {}
        fractions = []
        stepsize = (10.0 / len(data)) / len(data[condition])
        for replicate in data[condition]:
            progress += stepsize
            if window:
                window["--progress--"].update(progress)
                window["--status2--"].update(" ".join([condition, replicate]))
                window.read(timeout=50)

            for sample in list(data[condition][replicate].columns):
                prefix = sample[: sample.find("_")]
                fractnumber = int(prefix[3:])
                if fractnumber not in fractions:
                    fractions.append(fractnumber)
                    fracts_count[condition][fractnumber] = 1
                else:
                    fracts_count[condition][fractnumber] += 1
        fractions = sorted(fractions)
        fracts_con[condition] = fractions
        fracts_corr[condition] = [
            "Fr." + str(k)
            for k, v in fracts_count[condition].items()
            if v == max(fracts_count[condition].values())
        ]
    return fracts_con, fracts_count, fracts_corr, progress


def combine_median_std(
    data: dict[str, dict[str, pd.DataFrame]],
    fracts_con: dict[str, list[int]],
    window: sg.Window,
    progress: float,
) -> tuple[dict[str, dict[str, pd.DataFrame]], dict[str, pd.DataFrame], float]:
    """Compute median and standard deviation for each condition x fraction
    across the provided replicates."""
    data_median = {}
    data_std = {}
    stepsize = 5.0 / len(data)

    for condition in data:
        progress += stepsize
        if window:
            window["--progress--"].update(progress)
            window["--status2--"].update(condition)
            window.read(timeout=50)

        con_vals, con_std = combine_median_std_for_condition(
            data[condition], fracts_con[condition], condition
        )
        data_std[condition] = con_std
        # data_median[condition+'_median'] = con_vals
        data_median[condition] = {"median": con_vals}
    return data_median, data_std, progress


def combine_median_std_for_condition(
    data: dict[str, pd.DataFrame], fracts_con: list[int], condition: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute median and standard deviation for each fraction across the
    provided replicates.

    :param data: Data for one condition. ReplicateID -> DataFrame.
    :param fracts_con: List of fraction indices to consider.
    :param condition: Condition id for constructing column names.
    :return: A tuple containing the median and standard deviation dataframes.
    """
    con_vals = pd.DataFrame()
    con_std = pd.DataFrame()

    # TODO(performance): some merging can be done more efficiently
    for fract in fracts_con:
        fract_df = pd.DataFrame()
        frac_id = f"Fr.{fract}"

        for replicate_df in data.values():
            for sample in replicate_df:
                if sample[: sample.find("_")] == frac_id:
                    fract_df = pd.merge(
                        fract_df,
                        replicate_df[sample],
                        left_index=True,
                        right_index=True,
                        how="outer",
                    )

        col_median = f"{condition}_median_{frac_id}"
        col_std = f"{condition}_std_{frac_id}"

        con_vals = pd.merge(
            con_vals,
            fract_df.median(axis=1).rename(col_median),
            left_index=True,
            right_index=True,
            how="outer",
        )

        con_std = pd.merge(
            con_std,
            fract_df.std(axis=1).rename(col_std),
            left_index=True,
            right_index=True,
            how="outer",
        )

    return con_vals.fillna(0.0), con_std.fillna(0.0)


def combine_concat(data, window):
    for condition in data:
        if window:
            window["--status2--"].update(condition)
            window.read(timeout=50)

        con_vals = pd.DataFrame()
        for replicate in data[condition]:
            renamedict = {}
            for sample in data[condition][replicate]:
                oldname = sample
                newname = (
                    condition
                    + sample[sample.rfind("_") :]
                    + "_"
                    + sample[: sample.find("_")]
                )
                renamedict[oldname] = newname
            con_vals = pd.merge(
                con_vals,
                data[condition][replicate],
                left_index=True,
                right_index=True,
                how="outer",
            )
            con_vals.rename(renamedict, axis="columns", inplace=True)
        data[condition] = {"concat": con_vals.fillna(0.0)}
    return data


def remove_zeros(
    data: dict[str, dict[str, pd.DataFrame]],
) -> dict[str, dict[str, pd.DataFrame]]:
    """Remove rows with all zeros from the data."""
    for condition in data:
        for replicate, df in data[condition].items():
            df = df.apply(pd.to_numeric, errors="coerce")
            df = df.fillna(0)[~(df == 0).all(axis=1)]
            data[condition][replicate] = df
    return data


def calculate_outcorr(data, protlist_remaining, comb, window, progress):
    outer_corrs = pd.DataFrame()
    stepsize = 5.0 / len(data)

    for condition in data:
        progress += stepsize
        if window:
            window["--progress--"].update(progress)
            window["--status2--"].update(condition)
            window.read(timeout=50)

        outcorr = pd.DataFrame(index=protlist_remaining[condition])
        for con in data:
            if not con == condition:
                col_new = "OuterCorrelation_" + condition + "_" + con
                outcorr[col_new] = np.nan
                data_own = data[condition][comb].fillna(0.0)
                data_other = data[con][comb].fillna(0.0)
                fracts_own = []
                fracts_other = []
                for fract in data_own.columns:
                    fracts_own.append(fract[fract.rfind("_") + 1 :])
                for fract in data_other.columns:
                    fracts_other.append(fract[fract.rfind("_") + 1 :])
                fracts_both = [x for x in fracts_own if x in fracts_other]
                for fract in data_own.columns:
                    suffix = fract[fract.rfind("_") + 1 :]
                    if suffix not in fracts_both:
                        data_own = data_own.drop([fract], axis=1)
                for fract in data_other.columns:
                    suffix = fract[fract.rfind("_") + 1 :]
                    if suffix not in fracts_both:
                        data_other = data_other.drop([fract], axis=1)
                for ID in data_own.index:
                    if ID in data_other.index:
                        profile_own = data_own.loc[ID].tolist()
                        profile_other = data_other.loc[ID].tolist()
                        corr = pearsonr(profile_own, profile_other)
                        outcorr[col_new][ID] = corr[0]
        outer_corrs = pd.merge(
            outer_corrs,
            outcorr,
            left_index=True,
            right_index=True,
            how="outer",
        )
    return outer_corrs, progress


def modify_structure(data_in):
    data_out = {"class": {}, "vis": {}}
    for way in data_in:
        for condition in data_in[way]:
            for mode in data_in[way][condition]:
                data_out[way][f"{condition}_{mode}"] = data_in[way][condition][
                    mode
                ]
                # data_out[way][condition] = data_in[way][condition][mode]
    return data_out


def create_fract_processing_window() -> sg.Window:
    """Create fractionation data processing progress dialog."""
    layout = [
        [
            sg.Column(
                [
                    [
                        sg.ProgressBar(
                            100,
                            orientation="h",
                            size=(38, 25),
                            key="--progress--",
                        )
                    ],
                    [
                        sg.Text(
                            "-ready-",
                            font=("Arial", 9),
                            size=(60, 2),
                            key="--status1--",
                        )
                    ],
                    [
                        sg.Text(
                            "for run",
                            font=("Arial", 9),
                            size=(60, 2),
                            key="--status2--",
                        )
                    ],
                ],
                size=(420, 100),
            ),
            sg.Column(
                [
                    [
                        sg.Button(
                            "Start",
                            size=(15, 1),
                            key="--start--",
                            disabled=False,
                            enable_events=True,
                            button_color="dark green",
                        )
                    ],
                    [
                        sg.Button(
                            "Cancel",
                            size=(15, 1),
                            key="--cancel--",
                            disabled=False,
                            enable_events=True,
                            button_color="black",
                        )
                    ],
                ],
                size=(150, 70),
                vertical_alignment="top",
            ),
        ]
    ]
    return sg.Window(
        "Processing...",
        layout,
        size=(600, 120),
        modal=True,
    )


def start_fract_data_processing(
    fract_input: dict[str, FractDataset],
    preparams: dict[str, dict],
    window: sg.Window | None = None,
):
    """Start fractionation data processing."""
    # collect conditions (including [KEEP])
    conditions = unique_preserve_order(
        row[1]
        for dataset in fract_input.values()
        for row in dataset.table
        if row[1] != IDENTIFIER
    )

    # ---------------------------------------------------------------------
    logger.info("creating dataset...")
    progress = 0
    if window:
        window["--status1--"].update(value="creating dataset...")
        window.read(timeout=50)

    dataset, protein_info, progress = create_dataset(
        fract_input,
        conditions,
        window,
        progress,
    )
    data_ways = {
        "class": copy.deepcopy(dataset),
        "vis": copy.deepcopy(dataset),
    }
    shapes = ", ".join(
        [
            f"{cond_id}/{rep_id}: {df.shape}"
            for cond_id, cond_dict in dataset.items()
            for rep_id, df in cond_dict.items()
        ]
    )
    logger.debug(f"Fractionation dataset shapes are {shapes}")

    # ---------------------------------------------------------------------
    logger.info("converting dataset...")
    progress = 10
    if window:
        window["--status1--"].update(value="converting dataset...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    for way in data_ways:
        data_ways[way] = remove_zeros(data_ways[way])

    # ---------------------------------------------------------------------
    logger.info("pre-scaling...")
    progress = 20
    if window:
        window["--status1--"].update(value="pre-scaling...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    for way in data_ways:
        if preparams[way]["scale1"][0]:
            data_ways[way], progress = pre_post_scaling(
                data_ways[way],
                preparams[way]["scale1"][1],
                window,
                progress,
            )

    # ---------------------------------------------------------------------
    logger.info("filtering by missing fractions...")
    progress = 30
    if window:
        window["--status1--"].update(value="filtering by missing values...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    if preparams["global"]["missing"][0]:
        for way in data_ways:
            data_ways[way], progress = filter_missing(
                data_ways[way],
                int(preparams["global"]["missing"][1]),
                window,
                progress,
            )

    # ---------------------------------------------------------------------
    logger.info("finding IDs...")
    progress = 40
    if window:
        window["--status1--"].update(value="finding IDs...")
        window["--progress--"].update(40)
        window.read(timeout=50)

    for way in data_ways:
        data_ways[way], proteins_remaining, progress = filter_count(
            data_ways[way],
            int(preparams["global"]["minrep"][1]),
            window,
            progress,
        )

    # ---------------------------------------------------------------------
    logger.info("detecting samples...")
    progress = 50
    if window:
        window["--status1--"].update(value="detecting samples...")
        window["--progress--"].update(50)
        window.read(timeout=50)

    fracts_con, fracts_count, fracts_corr, progress = list_samples(
        data_ways["class"], window, progress
    )

    # ---------------------------------------------------------------------
    logger.info("combining data...")
    progress = 60
    if window:
        window["--status1--"].update(value="combining data...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    std_ways = {"class": [], "vis": []}

    for way in data_ways:
        data_combined, std_ways[way], progress = combine_median_std(
            data_ways[way], fracts_con, window, progress
        )
        if preparams[way]["combination"] == "median":
            data_ways[way] = data_combined
        elif preparams[way]["combination"] == "concat":
            data_ways[way] = combine_concat(data_ways[way], window)
        elif preparams[way]["combination"] == "separate":
            pass

    # ---------------------------------------------------------------------
    logger.info("post-scaling...")
    progress = 70
    if window:
        window["--status1--"].update(value="post-scaling...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    for way in data_ways:
        if preparams[way]["scale2"][0]:
            data_ways[way], progress = pre_post_scaling(
                data_ways[way],
                preparams[way]["scale2"][1],
                window,
                progress,
            )

    # ---------------------------------------------------------------------
    logger.info("removing zeros...")
    progress = 80
    if window:
        window["--status1--"].update(value="removing baseline profiles...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    for way in data_ways:
        if preparams[way]["zeros"]:
            data_ways[way] = remove_zeros(data_ways[way])

    # ---------------------------------------------------------------------
    logger.info("calculating outer correlations...")
    progress = 90
    if window:
        window["--status1--"].update(value="calculating outer correlations...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    if preparams["global"]["outcorr"]:
        outcorr, progress = calculate_outcorr(
            data_ways["vis"],
            proteins_remaining,
            preparams["vis"]["combination"],
            window,
            progress,
        )
        for column in outcorr.columns:
            protein_info[column] = outcorr[column]

    data_ways = modify_structure(data_ways)

    # ---------------------------------------------------------------------
    progress = 100
    if window:
        window["--status1--"].update(value="calculating outer correlations...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    logger.info("done!")

    return (
        data_ways,
        std_ways,
        protein_info,
        conditions,
    )


def sample_tables_are_valid(
    input_tables: dict[str, list[list[int | str]]],
    min_replicates: int,
) -> bool:
    """Check that the sample table is valid.

    And show error boxes if not.
    """

    # validate samples table
    if not all(
        any(sample[1] == IDENTIFIER for sample in input_table)
        for input_table in input_tables.values()
    ):
        messagebox.showerror(
            "Error",
            "At least one Identifier is missing.\n"
            "Please check for multiple import files.",
        )
        return False

    if any(
        sample[1] == ""
        for input_table in input_tables.values()
        for sample in input_table
    ):
        messagebox.showerror(
            "Error",
            "At least one row does not have a condition assigned.\n"
            "Please check for multiple imported files.",
        )
        return False

    if any(
        sample[2] == ""
        for input_table in input_tables.values()
        for sample in input_table
    ):
        messagebox.showerror(
            "Error",
            "At least one Replicate is missing.\n"
            "Please check for multiple imported files.",
        )
        return False

    if any(
        sample[3] == ""
        for input_table in input_tables.values()
        for sample in input_table
    ):
        messagebox.showerror(
            "Error",
            "At least one Fraction is missing.\n"
            "Please check for multiple imported files.",
        )
        return False

    if any(
        len(set(sample[2] for sample in input_table)) - 1 < min_replicates
        for input_table in input_tables.values()
    ):
        messagebox.showerror(
            "Error",
            "Not enough replicates! "
            "Load more replicates or reduce threshold in Parameters.",
        )
        return False

    return True


def FDP_exec(
    fract_input: dict[str, FractDataset],
    preparams: dict[str, dict],
    data_ways: dict[str, dict[str, pd.DataFrame]],
    std_ways: dict[str, dict[str, pd.DataFrame]],
    protein_info: dict[str, pd.DataFrame],
    conditions_trans: list[str],
):
    """Execute the Fractionation Data Processing."""
    window = create_fract_processing_window()

    while True:
        event, values = window.read()
        if event == "--cancel--" or event == sg.WIN_CLOSED:
            window.close()
            break

        if event == "--start--":
            window["--start--"].update(disabled=True)
            window["--cancel--"].update(disabled=True)

            if not sample_tables_are_valid(
                {k: dset.table for k, dset in fract_input.items()},
                min_replicates=int(preparams["global"]["minrep"][0]),
            ):
                break

            (
                data_ways,
                std_ways,
                protein_info,
                conditions_trans,
            ) = start_fract_data_processing(
                fract_input,
                preparams,
                window,
            )
            break

    window.close()

    return (
        data_ways,
        std_ways,
        protein_info,
        conditions_trans,
    )
