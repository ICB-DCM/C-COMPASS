"""Total proteome data processing."""

import logging
import math
from itertools import chain
from tkinter import messagebox
from typing import Any, Literal

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from ._utils import unique_preserve_order
from .core import IDENTIFIER, KEEP, TotalProtDataset

logger = logging.getLogger(__package__)


def create_dataset(
    tp_input: dict[str, TotalProtDataset], window: sg.Window | None
):
    tp_conditions = unique_preserve_order(
        sample[1]
        for dset in tp_input.values()
        for sample in dset.table
        if sample[1] != IDENTIFIER
    )

    identifiers = list(
        set(
            chain.from_iterable(
                dset.df[dset.id_col] for dset in tp_input.values()
            )
        )
    )
    combined = {}

    for condition in tp_conditions:
        data_new = pd.DataFrame(index=identifiers)
        for dset in tp_input.values():
            if window:
                window["--status2--"].update(condition)
                window.read(timeout=50)
            replicate = 1
            for sample in dset.table:
                data = pd.DataFrame()
                if sample[1] == condition:
                    samplename = sample[0]
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

        combined[condition] = data_new

    if KEEP in combined:
        data_keep = combined[KEEP]
        del combined[KEEP]
    else:
        data_keep = pd.DataFrame()

    return combined, data_keep


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
        if window:
            window["--status2--"].update(condition)
            window.read(timeout=50)
        # data[condition] = pd.to_numeric(data[condition], errors='coerce')
        data[condition] = np.log2(data[condition])
    return data


def impute_data(
    data, window: sg.Window | None, mode: Literal["normal", "constant"]
):
    s = 1.8
    w = 0.3
    for condition in data:
        if window:
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
    tp_input: dict[str, TotalProtDataset],
    tp_preparams: dict[str, Any],
    tp_data: dict[str, pd.DataFrame],
    tp_info: pd.DataFrame,
    tp_icorr: dict,
    window: sg.Window | None = None,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict]:
    # validate input
    if not all(
        any(IDENTIFIER == sample[1] for sample in dset.table)
        for dset in tp_input.values()
    ):
        messagebox.showerror("Error", "At least one Identifier is missing.")
        return tp_data, tp_info, tp_icorr

    if any(
        sample[1] == "" for dset in tp_input.values() for sample in dset.table
    ):
        messagebox.showerror(
            "Error", "At least one row does not have a condition assigned."
        )
        return tp_data, tp_info, tp_icorr

    if window:
        # deactivate buttons
        window["--start--"].update(disabled=True)
        window["--cancel--"].update(disabled=True)

    # ---------------------------------------------------------------------
    logger.info("creating dataset...")
    if window:
        window["--status1--"].update(value="creating dataset...")
        window.read(timeout=50)

    tp_data, tp_info = create_dataset(
        tp_input,
        window,
    )

    # ---------------------------------------------------------------------
    logger.info("filtering by missing values...")
    progress = 10
    if window:
        window["--status1--"].update(value="filtering by missing values...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    for df in tp_data.values():
        df.dropna(thresh=tp_preparams["minrep"], inplace=True)

    # ---------------------------------------------------------------------
    logger.info("transforming data...")
    progress = 20
    if window:
        window["--status1--"].update(value="transforming data...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    tp_data = transform_data(tp_data, window)

    # ---------------------------------------------------------------------
    logger.info("imputing MissingValues...")
    progress = 30
    if window:
        window["--status1--"].update(value="imputing MissingValues...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    tp_data = impute_data(tp_data, window, tp_preparams["imputation"])

    # ---------------------------------------------------------------------
    logger.info("calculating correlations...")
    progress = 40
    if window:
        window["--status1--"].update(value="calculating correlations...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    tp_icorr = calculate_correlations(tp_data)

    # ---------------------------------------------------------------------
    logger.info("normalizing data...")
    progress = 50
    if window:
        window["--status1--"].update(value="normalizing data...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    normalize_data(tp_data)

    logger.info("done!")
    progress = 60
    if window:
        window["--status1--"].update(value="normalizing data...")
        window["--progress--"].update(progress)
        window.read(timeout=50)

    return tp_data, tp_info, tp_icorr


def total_proteome_processing_dialog(
    tp_input: dict[str, TotalProtDataset],
    tp_preparams: dict[str, Any],
    tp_data: dict[str, pd.DataFrame],
    tp_info: pd.DataFrame,
    tp_icorr: dict,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict]:
    """Show the total proteome processing dialog."""
    window = create_window()

    while True:
        event, values = window.read()

        if event == "--cancel--" or event == sg.WIN_CLOSED:
            break

        if event == "--start--":
            tp_data, tp_info, tp_icorr = start_total_proteome_processing(
                tp_input,
                tp_preparams,
                tp_data,
                tp_info,
                tp_icorr,
                window,
            )
            break

    window.close()

    return tp_data, tp_info, tp_icorr
