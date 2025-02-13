"""Total proteome data processing."""

import logging
import math
from tkinter import messagebox
from typing import Any

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from ._utils import unique_preserve_order

logger = logging.getLogger(__package__)


def create_dataset(
    tp_indata, tp_tables, tp_identifiers, tp_conditions, window
):
    tp_conditions.remove("[IDENTIFIER]")

    idents = []
    for path in tp_tables:
        idents = list(
            set(idents + list(tp_indata[path][tp_identifiers[path]]))
        )

    dataset = {}

    for condition in tp_conditions:
        data_new = pd.DataFrame(index=idents)
        for path in tp_tables:
            window["--status2--"].update(condition)
            window.read(timeout=50)
            replicate = 1
            for sample in tp_tables[path]:
                data = pd.DataFrame()
                if sample[1] == condition:
                    samplename = sample[0]
                    data[samplename] = tp_indata[path][sample[0]]
                    data.set_index(
                        tp_indata[path][tp_identifiers[path]], inplace=True
                    )
                    data_new = pd.merge(
                        data_new,
                        data,
                        right_index=True,
                        left_index=True,
                        how="outer",
                    )
                    if condition == "[KEEP]":
                        if samplename + "_x" in data_new.columns:
                            for element in list(data_new.index):
                                if pd.isnull(
                                    data_new[samplename + "_x"][element]
                                ):
                                    data_new[samplename + "_x"][element] = (
                                        data_new[samplename + "_y"][element]
                                    )
                                if pd.isnull(
                                    data_new[samplename + "_y"][element]
                                ):
                                    data_new[samplename + "_y"][element] = (
                                        data_new[samplename + "_y"][element]
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
                                samplename: f"{sample[1]}_Rep.{replicate}"
                            }
                        )
                        replicate += 1

        if data_new.map(lambda x: "," in str(x)).any().any():
            data_new = data_new.map(
                lambda x: str(x).replace(",", ".") if isinstance(x, str) else x
            )
            data_new = data_new.apply(pd.to_numeric, errors="coerce")

        dataset[condition] = data_new

    if "[KEEP]" in dataset:
        data_keep = dataset["[KEEP]"]
        del dataset["[KEEP]"]
    else:
        data_keep = pd.DataFrame()

    return dataset, data_keep, tp_conditions


def filter_missing(data, mincount, window):
    for condition in data:
        window["--status2--"].update(condition)
        window.read(timeout=50)
        data[condition].dropna(thresh=mincount, inplace=True)
    return data


def calculate_correlations(data):
    tp_icorr = {}
    for condition in data:
        corrs = []
        data[condition].dropna(
            thresh=len(data[condition].columns), inplace=True
        )
        for rep_own in data[condition].columns.tolist():
            for rep_other in data[condition].columns.tolist():
                if not rep_own == rep_other:
                    corrs.append(
                        pearsonr(
                            data[condition][rep_own].tolist(),
                            data[condition][rep_other].tolist(),
                        )[0]
                    )
        tp_icorr[condition] = np.mean(corrs)

    return tp_icorr


def transform_data(data, window):
    for condition in data:
        window["--status2--"].update(condition)
        window.read(timeout=50)
        # data[condition] = pd.to_numeric(data[condition], errors='coerce')
        data[condition] = np.log2(data[condition])
    return data


def impute_data(data, window, mode):
    s = 1.8
    w = 0.3
    for condition in data:
        window["--status2--"].update(condition)
        window.read(timeout=50)
        if mode == "normal":
            for sample in data[condition]:
                mean = np.mean(data[condition][sample])
                std = np.std(data[condition][sample])
                mean_imp = mean - s * std
                sigma = std * w
                data[condition][sample] = data[condition][sample].apply(
                    lambda x: np.random.normal(mean_imp, sigma, 1)[0]
                    if math.isnan(x)
                    else x
                )
        elif mode == "constant":
            for sample in data[condition]:
                data[condition][sample] = data[condition][sample].apply(
                    lambda x: 0 if math.isnan(x) else x
                )
    return data


def normalize_data(data):
    for condition in data:
        for replicate in data[condition]:
            q1 = np.percentile(data[condition][replicate], 25)
            q2 = np.percentile(data[condition][replicate], 50)
            q3 = np.percentile(data[condition][replicate], 75)

            data[condition][replicate] = data[condition][replicate].apply(
                lambda x: (x - q2) / (q3 - q2)
                if x - q2 >= 0
                else (x - q2) / (q2 - q1)
            )


def create_window() -> sg.Window:
    """Create the total proteome processing dialog window."""
    layout = [
        [
            sg.Column(
                [
                    [
                        sg.ProgressBar(
                            60,
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
                            pad=(1, 1),
                            key="--status1--",
                        )
                    ],
                    [
                        sg.Text(
                            "for run",
                            font=("Arial", 9),
                            size=(60, 2),
                            pad=(1, 1),
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
            ),
        ]
    ]
    return sg.Window("Processing...", layout, size=(600, 110), modal=True)


def start_total_proteome_processing(
    window: sg.Window,
    tp_data: dict[str, pd.DataFrame],
    tp_tables: dict[str, list[tuple[str, str]]],
    tp_preparams: dict[str, Any],
    tp_identifiers: dict[str, str],
    tp_info: pd.DataFrame,
    tp_icorr: dict,
    tp_indata: dict[str, pd.DataFrame],
    tp_conditions: list,
):
    # validate input
    if not all(
        any("[IDENTIFIER]" == sample[1] for sample in table)
        for table in tp_tables.values()
    ):
        messagebox.showerror("Error", "At least one Identifier is missing.")
        return tp_data, tp_info, tp_conditions, tp_icorr

    if any(
        sample[1] == "" for table in tp_tables.values() for sample in table
    ):
        messagebox.showerror("Error", "At least one Condition is missing.")
        return tp_data, tp_info, tp_conditions, tp_icorr

    # deactivate buttons
    window["--start--"].update(disabled=True)
    window["--cancel--"].update(disabled=True)

    # ---------------------------------------------------------------------
    logger.info("creating dataset...")
    window["--status1--"].update(value="creating dataset...")
    window.read(timeout=50)

    conditions = unique_preserve_order(
        sample[1] for table in tp_tables.values() for sample in table
    )
    tp_data, tp_info, tp_conditions = create_dataset(
        tp_indata,
        tp_tables,
        tp_identifiers,
        conditions,
        window,
    )

    # ---------------------------------------------------------------------
    logger.info("filtering by missing values...")
    progress = 10
    window["--status1--"].update(value="filtering by missing values...")
    window["--progress--"].update(progress)
    window.read(timeout=50)

    tp_data = filter_missing(tp_data, tp_preparams["minrep"], window)

    # ---------------------------------------------------------------------
    logger.info("transforming data...")
    progress = 20
    window["--status1--"].update(value="transforming data...")
    window["--progress--"].update(progress)
    window.read(timeout=50)

    tp_data = transform_data(tp_data, window)

    # ---------------------------------------------------------------------
    logger.info("imputing MissingValues...")
    progress = 30
    window["--status1--"].update(value="imputing MissingValues...")
    window["--progress--"].update(progress)
    window.read(timeout=50)

    tp_data = impute_data(tp_data, window, tp_preparams["imputation"])

    # ---------------------------------------------------------------------
    logger.info("calculating correlations...")
    progress = 40
    window["--status1--"].update(value="calculating correlations...")
    window["--progress--"].update(progress)
    window.read(timeout=50)

    tp_icorr = calculate_correlations(tp_data)

    # ---------------------------------------------------------------------
    logger.info("normalizing data...")
    progress = 50
    window["--status1--"].update(value="normalizing data...")
    window["--progress--"].update(progress)
    window.read(timeout=50)

    normalize_data(tp_data)

    logger.info("done!")
    progress = 60
    window["--status1--"].update(value="normalizing data...")
    window["--progress--"].update(progress)
    window.read(timeout=50)

    return tp_data, tp_info, tp_conditions, tp_icorr


def total_proteome_processing_dialog(
    tp_data: dict[str, pd.DataFrame],
    tp_tables: dict[str, list[tuple[str, str]]],
    tp_preparams: dict[str, Any],
    tp_identifiers: dict[str, str],
    tp_info: pd.DataFrame,
    tp_icorr: dict,
    tp_indata: dict[str, pd.DataFrame],
    tp_conditions: list,
):
    """Show the total proteome processing dialog."""
    window = create_window()

    while True:
        event, values = window.read()

        if event == "--cancel--" or event == sg.WIN_CLOSED:
            break

        if event == "--start--":
            tp_data, tp_info, tp_conditions, tp_icorr = (
                start_total_proteome_processing(
                    window,
                    tp_data,
                    tp_tables,
                    tp_preparams,
                    tp_identifiers,
                    tp_info,
                    tp_icorr,
                    tp_indata,
                    tp_conditions,
                )
            )
            break

    window.close()

    return tp_data, tp_info, tp_conditions, tp_icorr
