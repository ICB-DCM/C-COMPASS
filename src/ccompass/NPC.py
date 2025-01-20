"""Normalized profile creation"""
# FIXME: this module is currently unused

import copy
from datetime import datetime
from tkinter import messagebox, simpledialog

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd

from .NPC_func import (
    calculate_icorr,
    calculate_ocorr,
    create_concat,
    create_dataset,
    create_median,
    create_separate,
    filter_count,
    filter_empty,
    filter_missing,
    implement_icorr,
    list_samples,
    pre_scaling,
    remove_worst,
    remove_zeros,
    scale_area,
)


def create_data(path):  # creates initial dataset with corresponding data table
    data = pd.read_csv(path, sep="\t", header=0)
    data = data.replace("NaN", np.nan)
    data = data.replace("Filtered", np.nan)
    colnames = data.columns.values.tolist()
    table = []
    for name in colnames:
        namelist = [name, "", "", ""]
        table.append(namelist)
    return (data, table)


def modify_table(
    title, prompt, values, tables_all, pos, q, ask
):  # modifies data table after making changes in GUI
    if values["--definition_table--"]:
        path = values["--definition_path--"]
        table = tables_all[path]
        if ask == "integer":
            value = simpledialog.askinteger(title, prompt)
            p = 0
            if value:
                for i in values["--definition_table--"]:
                    table[i][pos] = value + p
                    p = p + q
                tables_all[path] = table
        elif ask == "string":
            value = simpledialog.askstring(title, prompt)
            if value:
                for i in values["--definition_table--"]:
                    table[i][pos] = value
                tables_all[path] = table
    else:
        messagebox.showerror("Error", "No sample selected.")
    return (values, tables_all)


def create_stats(data_con, params):  # summarizes stats
    stats = {}
    index_names = []
    for condition in data_con:
        for replicate in data_con[condition]:
            naming = condition + "_" + replicate
            index_names.append(naming)
        index_names.append(condition + "_" + params["mode"])
    stats["filtered"] = pd.DataFrame(
        columns=[
            "by ValidValues",
            "by count",
            "by InnerCorrelation",
            "by baseline profiles",
        ],
        index=index_names,
    )
    return stats


def refresh_status(window, state, fullstate, text):
    window["--progress--"].Update(state, fullstate)
    window["--status--"].Update(text)
    event_NPC, values_NPC = window.read(timeout=5)


# -----------------------------------------------------------------------------
### RUN SCRIPT:
# -----------------------------------------------------------------------------
def NPC(data_con_in, data_keep_in, params_in, stats_in, data_con_std_in):
    # -----------------------------------------------------------------------------
    ### CREATE GUI:
    data_con = data_con_in
    params = params_in
    stats = stats_in
    data_keep = data_keep_in
    data_con_std = data_con_std_in

    data_all = {}
    tables_all = {}

    input_paths = []
    input_list = []
    ident_pos = {}
    identifier = {}

    layout_NPC = [
        [
            sg.Column(
                layout=[
                    [
                        sg.Frame(
                            layout=[  ## Data Import Frame
                                [
                                    sg.Listbox(
                                        size=(65, 8),
                                        select_mode="single",
                                        horizontal_scroll=True,
                                        values=[],
                                        key="--input_box--",
                                    ),
                                    sg.Column(
                                        layout=[
                                            [
                                                sg.Button(
                                                    "Open...",
                                                    enable_events=True,
                                                    key="--input_open--",
                                                    button_color="grey",
                                                    size=(6, 1),
                                                    disabled=False,
                                                )
                                            ],
                                            [
                                                sg.Button(
                                                    "Save...",
                                                    enable_events=True,
                                                    key="--input_save--",
                                                    button_color="grey",
                                                    size=(6, 1),
                                                    disabled=False,
                                                )
                                            ],
                                            [sg.Text("-" * 13)],
                                            [
                                                sg.Button(
                                                    "Load",
                                                    enable_events=True,
                                                    key="--input_load--",
                                                    button_color="dark green",
                                                    size=(6, 1),
                                                )
                                            ],
                                            [
                                                sg.Button(
                                                    "Reset",
                                                    enable_events=True,
                                                    key="--input_reset--",
                                                    disabled=True,
                                                    button_color="dark red",
                                                    size=(6, 1),
                                                )
                                            ],
                                        ]
                                    ),
                                ],
                                [
                                    sg.Button(
                                        "Add file",
                                        disabled=False,
                                        key="--input_add--",
                                        enable_events=True,
                                        size=(6, 1),
                                        button_color="grey",
                                    ),
                                    sg.Button(
                                        "Remove file",
                                        enable_events=True,
                                        key="--input_remove--",
                                        disabled=True,
                                    ),
                                ],
                            ],
                            title="Data Import",
                            size=(580, 220),
                        )
                    ],
                    [
                        sg.Frame(
                            layout=[  ## Operations Frame
                                [
                                    sg.Text("purpose: "),
                                    sg.Combo(
                                        ["for prediction", "for plots"],
                                        default_value="for prediction",
                                        size=(20, 1),
                                        enable_events=True,
                                        key="--preset--",
                                        readonly=True,
                                        disabled=True,
                                    ),
                                ],
                                [
                                    sg.Checkbox(
                                        "pre-scaling - over condition",
                                        default=False,
                                        enable_events=True,
                                        disabled=True,
                                        key="--scale1_condition--",
                                    ),
                                    sg.Checkbox(
                                        "pre-scaling - over replicate",
                                        default=True,
                                        enable_events=True,
                                        disabled=True,
                                        key="--scale1_replicate--",
                                    ),
                                ],
                                [
                                    sg.Checkbox(
                                        "minimum count of Valid Values (per gradient/replicate):  ",
                                        default=True,
                                        disabled=True,
                                        enable_events=True,
                                        key="--missing_active--",
                                    ),
                                    sg.Spin(
                                        values=(list(range(1, 100))),
                                        initial_value="1",
                                        size=(4, 2),
                                        readonly=False,
                                        disabled=True,
                                        key="--missing_number--",
                                    ),
                                ],
                                [
                                    sg.Checkbox(
                                        "found in at least: ",
                                        default=True,
                                        enable_events=True,
                                        disabled=True,
                                        key="--min_active--",
                                    ),
                                    sg.Spin(
                                        values=(list(range(1, 10))),
                                        initial_value="3",
                                        size=(4, 2),
                                        readonly=True,
                                        disabled=True,
                                        key="--min_count--",
                                    ),
                                    sg.Text(
                                        " replicates", text_color="light grey"
                                    ),
                                ],
                                [
                                    sg.Checkbox(
                                        "exclude proteins from worst correlated replicate",
                                        default=False,
                                        enable_events=False,
                                        disabled=True,
                                        key="--correl_rep_active--",
                                    )
                                ],
                                [
                                    sg.Checkbox(
                                        "median profile",
                                        default=False,
                                        enable_events=True,
                                        disabled=True,
                                        key="--median--",
                                        tooltip=" for consistent/reproducible replicates ",
                                    ),
                                    sg.Checkbox(
                                        "concatenated profiles",
                                        default=True,
                                        enable_events=True,
                                        disabled=True,
                                        key="--concat--",
                                        tooltip=" for variations over replicate ",
                                    ),
                                    sg.Checkbox(
                                        "process separately",
                                        default=False,
                                        enable_events=True,
                                        disabled=True,
                                        key="--separate--",
                                        tooltip=" for other purposes ",
                                    ),
                                ],
                                [
                                    sg.Checkbox(
                                        "post-scaling - MinMax (0-1)",
                                        default=False,
                                        enable_events=True,
                                        disabled=True,
                                        key="--scale2_0--",
                                    ),
                                    sg.Checkbox(
                                        "post-scaling - A=1",
                                        default=False,
                                        enable_events=True,
                                        disabled=True,
                                        key="--scale2_1--",
                                    ),
                                ],
                                [
                                    sg.Checkbox(
                                        "remove baseline profiles (zeroes)",
                                        default=False,
                                        disabled=True,
                                        key="--zeros_active--",
                                    )
                                ],
                                [
                                    sg.Checkbox(
                                        "calculate correlations between conditions",
                                        default=False,
                                        disabled=True,
                                        enable_events=False,
                                        key="--correl_con_active--",
                                    )
                                ],
                            ],
                            title="Operations",
                            size=(580, 290),
                        ),
                    ],
                ]
            ),
            sg.Column(
                layout=[
                    [
                        sg.Frame(
                            layout=[  ## Data Definition Frame
                                [
                                    sg.Combo(
                                        (input_paths),
                                        size=(75, 1),
                                        enable_events=True,
                                        key="--definition_path--",
                                        readonly=True,
                                        disabled=True,
                                    )
                                ],
                                [
                                    sg.Table(
                                        values=[],
                                        num_rows=24,
                                        headings=[
                                            "Sample",
                                            "Condition",
                                            "Replicate",
                                            "Fraction",
                                        ],
                                        auto_size_columns=False,
                                        max_col_width=20,
                                        vertical_scroll_only=False,
                                        col_widths=[34, 10, 7, 7],
                                        key="--definition_table--",
                                    )
                                ],
                                [
                                    sg.Button(
                                        "Remove",
                                        enable_events=True,
                                        key="--definition_remove--",
                                        disabled=True,
                                        button_color="dark red",
                                    ),
                                    sg.Button(
                                        "Keep",
                                        enable_events=True,
                                        key="--definition_keep--",
                                        disabled=True,
                                        button_color="dark grey",
                                        tooltip=" Try to keep gene names! ",
                                    ),
                                    sg.Button(
                                        "Set Condition",
                                        enable_events=True,
                                        key="--definition_edit_condition--",
                                        disabled=True,
                                    ),
                                    sg.Button(
                                        "Set Replicate",
                                        enable_events=True,
                                        key="--definition_edit_replicate--",
                                        disabled=True,
                                    ),
                                    sg.Button(
                                        "Set Fractions",
                                        enable_events=True,
                                        key="--definition_edit_fractions--",
                                        disabled=True,
                                    ),
                                    sg.Button(
                                        "Set Identifier",
                                        enable_events=True,
                                        key="--definition_identifier--",
                                        disabled=True,
                                        button_color="grey",
                                        tooltip=" If possible, use protein groups! ",
                                    ),
                                ],
                            ],
                            title="Data Definition",
                            size=(580, 516),
                        )
                    ]
                ]
            ),
        ],
        [
            sg.Button(
                "Process!",
                enable_events=True,
                key="--start--",
                disabled=True,
                button_color="black",
                size=(8, 1),
            ),
            sg.Button(
                "Accept",
                disabled=True,
                enable_events=True,
                key="--exit--",
                button_color="dark green",
                size=(8, 1),
            ),
            sg.Button(
                "Export",
                enable_events=True,
                key="--export--",
                button_color="grey",
                disabled=True,
                size=(8, 1),
            ),
            sg.ProgressBar(
                100,
                orientation="h",
                size=(60, 2),
                key="--progress--",
                border_width=2,
                visible=False,
            ),
            sg.Text(
                "ready...", key="--status--", visible=False, font=("Arial", 8)
            ),
        ],
        [
            sg.InputText(
                "myExperiment",
                key="--experiment_name--",
                disabled=False,
                readonly=False,
            )
        ],
    ]

    window_NPC = sg.Window(
        "Normalized Profile Creator", layout_NPC, size=(1200, 600)
    )

    # -----------------------------------------------------------------------------
    ### RUN GUI:
    while True:
        event_NPC, values_NPC = window_NPC.read()

        if event_NPC == "--input_add--":
            filename = sg.popup_get_file(
                "Chose dataset",
                no_window=True,
                file_types=(
                    ("Tab Separated Values", "*.tsv"),
                    ("Text (tab delimited)", "*.txt"),
                ),
            )
            if filename:
                input_paths.append(filename)
                window_NPC["--input_box--"].Update(values=input_paths)
                for i in ["--input_load--", "--input_remove--"]:
                    window_NPC[i].Update(disabled=False)

        if event_NPC == "--input_remove--":  # Remove Button
            if window_NPC["--input_box--"].get_indexes():
                selected = window_NPC["--input_box--"].get_indexes()[0]
                rempath = window_NPC["--input_box--"].get_list_values()[
                    selected
                ]
                input_paths.remove(rempath)
                window_NPC["--input_box--"].Update(values=input_paths)

            else:
                messagebox.showerror("Error", "No file selected.")

        if event_NPC == "--input_load--":  # Load Button
            if input_paths:
                input_list.extend(input_paths)
                window_NPC["--definition_path--"].Update(
                    values=input_paths, value=input_paths[0]
                )
                for path in input_paths:
                    data_all[path], tables_all[path] = create_data(path)
                    ident_pos[path] = []
                window_NPC["--definition_table--"].Update(
                    values=tables_all[input_list[0]]
                )
                for i in ["--input_load--", "--input_remove--"]:
                    window_NPC[i].Update(disabled=True)
                for i in [
                    "--input_reset--",
                    "--definition_path--",
                    "--definition_remove--",
                    "--definition_keep--",
                    "--definition_edit_condition--",
                    "--definition_edit_replicate--",
                    "--definition_edit_fractions--",
                    "--definition_identifier--",
                    "--preset--",
                    "--missing_number--",
                    "--min_count--",
                    "--min_active--",
                    "--missing_active--",
                    #'--preset--', '--scale1_condition--', '--scale1_replicate--', '--missing_active--', '--missing_number--', '--min_active--', '--min_count--', '--correl_rep_active--', '--median--', '--concat--', '--separate--', '--scale2_0--', '--scale2_1--', '--zeros_active--', '--correl_con_active--',
                    "--start--",
                ]:
                    window_NPC[i].Update(disabled=False)
                input_paths = []
                window_NPC["--input_box--"].Update(values=input_paths)
            else:
                messagebox.showerror("Error", "No files selected.")

        if event_NPC == "--input_save--":  # save input information
            save_data = {}
            save_data["tables_all"] = tables_all
            save_data["input_list"] = input_list
            save_data["identifier"] = identifier
            save_data["data_all"] = data_all
            filename = sg.popup_get_file(
                "Save Settings",
                no_window=True,
                file_types=(("Numpy", "*.npy"),),
                save_as=True,
            )
            if filename:
                np.save(filename, save_data)

        if event_NPC == "--input_open--":  # open input information
            filename = sg.popup_get_file(
                "Open Settings",
                no_window=True,
                file_types=(("Numpy", "*.npy"),),
            )
            if filename:
                file = np.load(filename, allow_pickle="TRUE").item()
                tables_all = file["tables_all"]
                input_list = file["input_list"]
                identifier = file["identifier"]
                data_all = file["data_all"]
                if input_list:
                    window_NPC["--definition_table--"].Update(
                        values=tables_all[input_list[0]]
                    )
                    window_NPC["--definition_path--"].Update(
                        values=input_list, value=input_list[0]
                    )
                    input_paths = []
                    window_NPC["--input_box--"].Update(values=input_paths)
                    for i in ["--input_load--", "--input_remove--"]:
                        window_NPC[i].Update(disabled=True)
                    for i in [
                        "--input_reset--",
                        "--definition_path--",
                        "--definition_remove--",
                        "--definition_keep--",
                        "--definition_edit_condition--",
                        "--definition_edit_replicate--",
                        "--definition_edit_fractions--",
                        "--definition_identifier--",
                        "--preset--",
                        "--missing_number--",
                        "--min_count--",
                        "--min_active--",
                        "--missing_active--",
                        #'--preset--', '--scale1_condition--', '--scale1_replicate--', '--missing_active--', '--min_active--', '--missing_number--', '--min_count--', '--correl_rep_active--', '--median--', '--concat--', '--separate--', '--scale2_0--', '--scale2_1--', '--zeros_active--', '--correl_con_active--',
                        "--start--",
                    ]:
                        window_NPC[i].Update(disabled=False)
                else:
                    messagebox.showerror("Error", "This file is empty.")

        if event_NPC == "--input_reset--":  # Reset Button
            input_paths = []
            input_list = []
            ident_pos = {}
            identifier = {}
            window_NPC["--definition_table--"].Update(values=[])
            window_NPC["--definition_path--"].Update(values=input_list)
            for i in [
                "--input_load--",
                "--input_reset--",
                "--input_remove--",
                "--definition_path--",
                "--definition_remove--",
                "--definition_keep--",
                "--definition_edit_condition--",
                "--definition_edit_replicate--",
                "--definition_edit_fractions--",
                "--definition_identifier--",
                "--preset--",
                "--scale1_condition--",
                "--scale1_replicate--",
                "--missing_active--",
                "--min_active--",
                "--missing_number--",
                "--min_count--",
                "--correl_rep_active--",
                "--median--",
                "--concat--",
                "--separate--",
                "--scale2_0--",
                "--scale2_1--",
                "--zeros_active--",
                "--correl_con_active--",
                "--start--",
            ]:
                window_NPC[i].Update(disabled=True)
            window_NPC["--input_box--"].Update(values=[])

        if (
            event_NPC == "--definition_path--"
        ):  # Buttons for condition definition
            window_NPC["--definition_table--"].Update(
                values=tables_all[values_NPC["--definition_path--"]]
            )
        if event_NPC == "--definition_edit_replicate--":
            values_NPC, tables_all = modify_table(
                "Set Replicate",
                "Replicate Number:",
                values_NPC,
                tables_all,
                2,
                0,
                "integer",
            )
            window_NPC["--definition_table--"].Update(
                values=tables_all[values_NPC["--definition_path--"]]
            )
        if event_NPC == "--definition_edit_fractions--":
            values_NPC, tables_all = modify_table(
                "Set Fractions",
                "FIRST Fraction Number:",
                values_NPC,
                tables_all,
                3,
                1,
                "integer",
            )
            window_NPC["--definition_table--"].Update(
                values=tables_all[values_NPC["--definition_path--"]]
            )
        if event_NPC == "--definition_edit_condition--":
            values_NPC, tables_all = modify_table(
                "Set Condition",
                "Condition Name:",
                values_NPC,
                tables_all,
                1,
                0,
                "string",
            )
            window_NPC["--definition_table--"].Update(
                values=tables_all[values_NPC["--definition_path--"]]
            )

        if event_NPC == "--definition_identifier--":
            pos = values_NPC["--definition_table--"]
            if pos:
                if len(pos) > 1:
                    messagebox.showerror(
                        "Error", "Please set only one Identifier!"
                    )
                elif len(pos) == 1:
                    path = values_NPC["--definition_path--"]
                    table = tables_all[path]
                    if ident_pos[path]:
                        table[ident_pos[path][0]][1] = ""
                        table[ident_pos[path][0]][2] = ""
                        table[ident_pos[path][0]][3] = ""
                    identifier[path] = table[pos[0]][0]
                    ident_pos[path] = pos
                    table[pos[0]][1] = "[IDENTIFIER]"
                    table[pos[0]][2] = "-"
                    table[pos[0]][3] = "-"
                    tables_all[path] = table
                    window_NPC["--definition_table--"].Update(
                        values=tables_all[values_NPC["--definition_path--"]]
                    )
            else:
                messagebox.showerror("Error", "No sample selected.")

        if event_NPC == "--definition_keep--":
            path = values_NPC["--definition_path--"]
            table = tables_all[path]
            for pos in values_NPC["--definition_table--"]:
                table[pos][1] = "[KEEP]"
                table[pos][2] = "-"
                table[pos][3] = "-"
            tables_all[path] = table
            window_NPC["--definition_table--"].Update(
                values=tables_all[values_NPC["--definition_path--"]]
            )

        if event_NPC == "--definition_remove--":
            path = values_NPC["--definition_path--"]
            selected = values_NPC["--definition_table--"]
            table = tables_all[path]
            for index in sorted(selected, reverse=True):
                del table[index]
            tables_all[path] = table
            window_NPC["--definition_table--"].Update(values=tables_all[path])

        if event_NPC == "--scale1_condition--":
            if values_NPC["--scale1_condition--"] == True:
                window_NPC["--scale1_replicate--"].Update(value=False)
        if event_NPC == "--scale1_replicate--":
            if values_NPC["--scale1_replicate--"] == True:
                window_NPC["--scale1_condition--"].Update(value=False)
        if event_NPC == "--missing_active--":
            window_NPC["--missing_number--"].Update(
                visible=values_NPC["--missing_active--"]
            )

        if event_NPC == "--min_active--":
            window_NPC["--min_count--"].Update(
                disabled=not values_NPC["--min_active--"]
            )
        if event_NPC == "--median--":
            window_NPC["--median--"].Update(value=True)
            window_NPC["--concat--"].Update(value=False)
            window_NPC["--separate--"].Update(value=False)
            window_NPC["--scale2_0--"].Update(disabled=False, value=True)
            window_NPC["--scale2_1--"].Update(disabled=False, value=False)
            window_NPC["--correl_con_active--"].Update(
                disabled=False, value=True
            )
        if event_NPC == "--concat--":
            window_NPC["--concat--"].Update(value=True)
            window_NPC["--median--"].Update(value=False)
            window_NPC["--separate--"].Update(value=False)
            window_NPC["--scale2_0--"].Update(disabled=False, value=False)
            window_NPC["--scale2_1--"].Update(disabled=True, value=False)
            window_NPC["--correl_con_active--"].Update(
                disabled=True, value=False
            )
        if event_NPC == "--separate--":
            window_NPC["--separate--"].Update(value=True)
            window_NPC["--median--"].Update(value=False)
            window_NPC["--concat--"].Update(value=False)
            window_NPC["--scale2_0--"].Update(disabled=False, value=False)
            window_NPC["--scale2_1--"].Update(disabled=False, value=False)
            window_NPC["--correl_con_active--"].Update(
                disabled=True, value=False
            )

        if event_NPC == "--scale2_0--":
            if values_NPC["--scale2_0--"] == True:
                window_NPC["--scale2_1--"].Update(value=False)
        if event_NPC == "--scale2_1--":
            if values_NPC["--scale2_1--"] == True:
                window_NPC["--scale2_0--"].Update(value=False)

        if event_NPC == "--preset--":
            if values_NPC["--preset--"] == "for plots":
                window_NPC["--scale1_condition--"].Update(
                    value=False, disabled=False
                )
                window_NPC["--scale1_replicate--"].Update(
                    value=True, disabled=False
                )
                window_NPC["--min_active--"].Update(value=True)
                window_NPC["--min_count--"].Update(value=3)
                window_NPC["--missing_active--"].Update(value=True)
                window_NPC["--missing_number--"].Update(value=1)
                window_NPC["--correl_rep_active--"].Update(
                    value=True, disabled=False
                )
                window_NPC["--median--"].Update(value=True, disabled=False)
                window_NPC["--concat--"].Update(value=False, disabled=False)
                window_NPC["--separate--"].Update(value=False, disabled=False)
                window_NPC["--scale2_0--"].Update(value=True, disabled=False)
                window_NPC["--scale2_1--"].Update(value=False, disabled=False)
                window_NPC["--zeros_active--"].Update(
                    value=True, disabled=False
                )
                window_NPC["--correl_con_active--"].Update(
                    value=True, disabled=False
                )
            if values_NPC["--preset--"] == "for prediction":
                window_NPC["--scale1_condition--"].Update(
                    value=False, disabled=True
                )
                window_NPC["--scale1_replicate--"].Update(
                    value=True, disabled=True
                )
                window_NPC["--min_active--"].Update(value=True)
                window_NPC["--min_count--"].Update(value=3)
                window_NPC["--missing_active--"].Update(value=True)
                window_NPC["--missing_number--"].Update(value=1)
                window_NPC["--correl_rep_active--"].Update(
                    value=False, disabled=True
                )
                window_NPC["--median--"].Update(value=False, disabled=True)
                window_NPC["--concat--"].Update(value=True, disabled=True)
                window_NPC["--separate--"].Update(value=False, disabled=True)
                window_NPC["--scale2_0--"].Update(value=False, disabled=True)
                window_NPC["--scale2_1--"].Update(value=False, disabled=True)
                window_NPC["--zeros_active--"].Update(
                    value=False, disabled=True
                )
                window_NPC["--correl_con_active--"].Update(
                    value=False, disabled=True
                )

        if event_NPC == sg.WIN_CLOSED:
            data_con_out = data_con_in
            params_out = params_in
            stats_out = stats_in
            data_con_std_out = data_con_std_in
            break
        if event_NPC == "--exit--":
            try:
                data_con_out = data_con
            except Exception:
                data_con_out = data_con_in
            try:
                params_out = params
            except Exception:
                params_out = params_in
            try:
                stats_out = stats
            except Exception:
                stats_out = stats_in
            try:
                data_con_std_out = data_con_std
            except Exception:
                data_con_std_out = data_con_std_in
            break

        if event_NPC == "--export--":
            export_folder = path = sg.popup_get_folder("Select folder")
            if export_folder:
                experiment = values_NPC["--experiment_name--"]
                now = datetime.now()
                time = "\\" + now.strftime("%Y%m%d%H%M%S") + "_"
                data_export = data_keep
                for condition in data_con:
                    data_export = pd.merge(
                        data_export,
                        data_con[condition],
                        left_index=True,
                        right_index=True,
                        how="outer",
                    )
                export_name = "AllConditions"
                export_path = (
                    export_folder
                    + time
                    + experiment
                    + " "
                    + export_name
                    + ".txt"
                )
                data_export.to_csv(
                    export_path,
                    header=True,
                    index=True,
                    index_label="Identifier",
                    sep="\t",
                    mode="a",
                )
                stats["filtered"].to_csv(
                    export_folder + time + experiment + "_filtered" + ".txt",
                    header=True,
                    index=True,
                    index_label="Identifier",
                    sep="\t",
                    mode="a",
                )

        # -----------------------------------------------------------------------------
        ## RUN ANALYSIS:

        if event_NPC == "--start--":  ## START Button
            is_ident = True  # check if identifier was set
            is_con = True
            is_rep = True
            is_fract = True
            fract_ok = True
            conditions = []

            for path in tables_all:
                conlist = []
                replist = []
                fractlist = []
                for sample in tables_all[path]:
                    conlist.append(sample[1])
                    replist.append(sample[2])
                    fractlist.append(sample[3])
                    if sample[1] not in conditions:
                        conditions.append(sample[1])
                if "[IDENTIFIER]" not in conlist:
                    is_ident = False
                if "" in conlist:
                    is_con = False
                if "" in replist:
                    is_rep = False
                if "" in fractlist:
                    is_fract = False
                if len(list(set(replist))) - 1 < values_NPC["--min_count--"]:
                    fract_ok = False

            if is_ident and is_con and is_rep and is_fract and fract_ok:
                # ------------------------------------------------
                window_NPC["--progress--"].Update(visible=True)
                window_NPC["--status--"].Update(visible=True)
                state = 0
                fullstate = 14
                # ------------------------------------------------

                # create dataset:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC, state, fullstate, "creating dataset..."
                )
                # ------------------------------------------------
                data_con_abs, data_keep, params = create_dataset(
                    data_all, tables_all, identifier, conditions, values_NPC
                )
                data_con = copy.deepcopy(data_con_abs)

                # scale data (pre-scaling):
                # ------------------------------------------------
                state += 1
                refresh_status(window_NPC, state, fullstate, "pre-scaling...")
                # ------------------------------------------------
                if values_NPC["--scale1_condition--"] == True:
                    mode = "condition"
                elif values_NPC["--scale1_replicate--"] == True:
                    mode = "replicate"
                else:
                    mode = "none"
                data_con_scaled = pre_scaling(data_con, data_con_abs, mode)
                data_con = copy.deepcopy(data_con_scaled)

                stats = create_stats(data_con, params)

                # filter by empty profiles:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "filtering by empty profiles...",
                )
                # ------------------------------------------------
                data_con_fempty = filter_empty(data_con)
                data_con = copy.deepcopy(data_con_fempty)

                # filter by missing values:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "filtering by missing values...",
                )
                # ------------------------------------------------
                if values_NPC["--missing_active--"] == True:
                    data_con_fmissing = filter_missing(
                        data_con, stats, values_NPC
                    )
                    data_con = copy.deepcopy(data_con_fmissing)

                # filter by count over replicates:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "filtering by replicate coverage...",
                )
                # ------------------------------------------------
                if values_NPC["--min_active--"] == True:
                    mincount = values_NPC["--min_count--"]
                    data_con_fcount, protlist_con, stats = filter_count(
                        data_con, mincount, stats
                    )
                    data_con = copy.deepcopy(data_con_fcount)
                else:
                    mincount = 1

                # find and list all fractions in each condition:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "finding and listing fractions...",
                )
                # ------------------------------------------------
                fracts_con, fracts_count, fracts_corr = list_samples(data_con)

                # calculate inner correlations:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "calculating inner correlations...",
                )
                # ------------------------------------------------
                icorr = calculate_icorr(data_con, fracts_corr, protlist_con)

                # remove worst correlated profiles:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "removing worst correlations...",
                )
                # ------------------------------------------------
                if values_NPC["--correl_rep_active--"] == True:
                    data_con_cleaned = remove_worst(
                        data_con, protlist_con, mincount, stats, icorr
                    )
                    data_con = copy.deepcopy(data_con_cleaned)

                    # re-calculate inner correlations:
                    # ------------------------------------------------
                    state += 1
                    refresh_status(
                        window_NPC,
                        state,
                        fullstate,
                        "re-calculating inner correlations...",
                    )
                    # ------------------------------------------------
                    icorr = calculate_icorr(
                        data_con, fracts_corr, protlist_con
                    )
                else:
                    state += 1

                # copy inner correlations to dataset:
                data_keep = implement_icorr(data_keep, icorr)

                # create median profiles & scale:
                if values_NPC["--median--"] == True:
                    # ------------------------------------------------
                    state += 1
                    refresh_status(
                        window_NPC,
                        state,
                        fullstate,
                        "creating median profiles...",
                    )
                    # ------------------------------------------------
                    data_con_median, data_con_std = create_median(
                        data_con, fracts_con, values_NPC["--scale2_0--"]
                    )
                    data_con = copy.deepcopy(data_con_median)

                # create concatenated profiles & scale:
                elif values_NPC["--concat--"] == True:
                    # ------------------------------------------------
                    state += 1
                    refresh_status(
                        window_NPC,
                        state,
                        fullstate,
                        "creating concatenated profiles...",
                    )
                    # ------------------------------------------------
                    data_con_concat = create_concat(
                        data_con, values_NPC["--scale2_0--"]
                    )
                    data_con = copy.deepcopy(data_con_concat)

                # create separate profiles:
                elif values_NPC["--separate--"] == True:
                    # ------------------------------------------------
                    state += 1
                    refresh_status(
                        window_NPC,
                        state,
                        fullstate,
                        "creating separated profiles...",
                    )
                    # ------------------------------------------------
                    data_con_sep = create_separate(data_con)
                    data_con = copy.deepcopy(data_con_sep)

                # scale data by area:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC, state, fullstate, "scaling profiles..."
                )
                # ------------------------------------------------
                if values_NPC["--scale2_1--"] == True:
                    data_con_area = scale_area(data_con)
                    data_con = copy.deepcopy(data_con_area)

                # remove baseline profiles:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "removing baseline profiles...",
                )
                # ------------------------------------------------
                if values_NPC["--zeros_active--"] == True:
                    data_con_nozeros = remove_zeros(data_con, stats, params)
                    data_con = copy.deepcopy(data_con_nozeros)

                # caluclate outer correlations:
                # ------------------------------------------------
                state += 1
                refresh_status(
                    window_NPC,
                    state,
                    fullstate,
                    "calculating outer correlations...",
                )
                # ------------------------------------------------
                if values_NPC["--correl_con_active--"] == True:
                    data_keep = calculate_ocorr(data_con, data_keep)

                # ------------------------------------------------
                state += 1
                refresh_status(window_NPC, state, fullstate, "done!")
                window_NPC["--progress--"].Update(bar_color=("green", "grey"))
                window_NPC["--export--"].Update(disabled=False)
                # ------------------------------------------------

                print("preprocessing finished!")

            else:
                if not is_ident:
                    messagebox.showerror(
                        "Error", "At least one Identifier is missing."
                    )
                elif not is_con:
                    messagebox.showerror(
                        "Error", "At least one Condition is missing."
                    )
                elif not is_rep:
                    messagebox.showerror(
                        "Error", "At least one Replicate is missing."
                    )
                elif not is_fract:
                    messagebox.showerror(
                        "Error", "At least one Fraction is missing."
                    )

            stats["additional info"] = data_keep
            stats["InnerCorrelations"] = icorr
            window_NPC["--exit--"].Update(disabled=False)

    window_NPC.close()
    return data_con_out, data_keep, params_out, stats_out, data_con_std_out
