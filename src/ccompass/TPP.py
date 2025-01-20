"""Total proteome data processing."""

import copy
import math
from tkinter import messagebox
from typing import Any

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


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
            window["--status2--"].Update(condition)
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
        window["--status2--"].Update(condition)
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
        window["--status2--"].Update(condition)
        window.read(timeout=50)
        # data[condition] = pd.to_numeric(data[condition], errors='coerce')
        data[condition] = np.log2(data[condition])
    return data


def impute_data(data, window, mode):
    s = 1.8
    w = 0.3
    for condition in data:
        window["--status2--"].Update(condition)
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


def normalize_data(data, window):
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
    return data


def create_window() -> sg.Window:
    """Create the total proteome processing dialog window."""
    layout_TPP = [
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
    return sg.Window("Processing...", layout_TPP, size=(600, 110), modal=True)


def start_total_proteome_processing(
    window_TPP: sg.Window,
    tp_data: dict[str, pd.DataFrame],
    tp_tables: dict[str, list[tuple[str, str]]],
    tp_preparams: dict[str, Any],
    tp_identifiers: dict[str, str],
    tp_intermediate: dict[str, dict[str, pd.DataFrame]],
    tp_info: pd.DataFrame,
    tp_icorr: dict,
    tp_indata: dict[str, pd.DataFrame],
    tp_conditions: list,
):
    is_ident = True
    is_con = True

    conditions = []
    for path in tp_tables:
        for sample in tp_tables[path]:
            if sample[1] not in conditions:
                conditions.append(sample[1])
        if "[IDENTIFIER]" not in conditions:
            is_ident = False
        if "" in conditions:
            is_con = False

    if not is_ident:
        messagebox.showerror("Error", "At least one Identifier is missing.")
    if not is_con:
        messagebox.showerror("Error", "At least one Condition is missing.")

    else:
        window_TPP["--start--"].Update(disabled=True)
        window_TPP["--cancel--"].Update(disabled=True)

        tp_intermediate = {}

        # ---------------------------------------------------------------------
        print("creating dataset...")
        progress = 0
        window_TPP["--status1--"].Update(value="creating dataset...")
        window_TPP.read(timeout=50)

        tp_data, tp_info, tp_conditions = create_dataset(
            tp_indata,
            tp_tables,
            tp_identifiers,
            conditions,
            window_TPP,
        )
        tp_intermediate["[0] abs"] = copy.deepcopy(tp_data)

        # ---------------------------------------------------------------------
        print("filtering by missing values...")
        progress = 10
        window_TPP["--status1--"].Update(
            value="filtering by missing values..."
        )
        window_TPP["--progress--"].Update(progress)
        window_TPP.read(timeout=50)

        minrep = tp_preparams["minrep"]
        tp_data = filter_missing(tp_data, minrep, window_TPP)
        tp_intermediate["[1] f_missing"] = copy.deepcopy(tp_data)

        # ---------------------------------------------------------------------
        print("transforming data...")
        progress = 20
        window_TPP["--status1--"].Update(value="transforming data...")
        window_TPP["--progress--"].Update(progress)
        window_TPP.read(timeout=50)

        tp_data = transform_data(tp_data, window_TPP)
        tp_intermediate["[2] transformed"] = copy.deepcopy(tp_data)

        # ---------------------------------------------------------------------
        print("imputing MissingValues...")
        progress = 30
        window_TPP["--status1--"].Update(value="imputing MissingValues...")
        window_TPP["--progress--"].Update(progress)
        window_TPP.read(timeout=50)

        tp_data = impute_data(tp_data, window_TPP, tp_preparams["imputation"])
        tp_intermediate["[3] imputed"] = copy.deepcopy(tp_data)

        # ---------------------------------------------------------------------
        print("calculating correlations...")
        progress = 40
        window_TPP["--status1--"].Update(value="calculating correlations...")
        window_TPP["--progress--"].Update(progress)
        window_TPP.read(timeout=50)

        tp_icorr = calculate_correlations(tp_data)

        # ---------------------------------------------------------------------
        print("normalizing data...")
        progress = 50
        window_TPP["--status1--"].Update(value="normalizing data...")
        window_TPP["--progress--"].Update(progress)
        window_TPP.read(timeout=50)

        tp_data = normalize_data(tp_data, window_TPP)
        tp_intermediate["[4] normalized"] = copy.deepcopy(tp_data)

        print("done!")
        progress = 60
        window_TPP["--status1--"].Update(value="normalizing data...")
        window_TPP["--progress--"].Update(progress)
        window_TPP.read(timeout=50)

    return tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr


def total_proteome_processing_dialog(
    tp_data: dict[str, pd.DataFrame],
    tp_tables: dict[str, list[tuple[str, str]]],
    tp_preparams: dict[str, Any],
    tp_identifiers: dict[str, str],
    tp_intermediate: dict[str, dict[str, pd.DataFrame]],
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
            tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr = (
                start_total_proteome_processing(
                    window,
                    tp_data,
                    tp_tables,
                    tp_preparams,
                    tp_identifiers,
                    tp_intermediate,
                    tp_info,
                    tp_icorr,
                    tp_indata,
                    tp_conditions,
                )
            )
            break

    window.close()

    return tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr
