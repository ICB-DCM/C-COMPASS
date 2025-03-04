"""C-COMPASS main window."""

import copy
import logging
import pickle
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from tkinter import messagebox, simpledialog
from typing import Any, Literal

import FreeSimpleGUI as sg
import pandas as pd

from . import MOA, RP
from ._gui_utils import wait_cursor
from .core import (
    IDENTIFIER,
    KEEP,
    NA,
    AppSettings,
    FractDataset,
    MarkerSet,
    SessionModel,
    SessionStatusModel,
    TotalProtDataset,
    app_name,
    create_fullprofiles,
    create_identity_conversion,
    create_marker_profiles,
    create_markerlist,
    read_fract_table,
    read_marker_table,
    read_tp_table,
    readthedocs_url,
    repository_url,
    write_class_changes_reports,
    write_comparison_reports,
    write_global_changes_reports,
    write_statistics_reports,
)

logger = logging.getLogger(__package__)


def create_fractionation_tab(fract_paths: Sequence[Any]) -> sg.Tab:
    """Create the "Fractionation" tab."""
    layout_fractionation = [
        [
            sg.Button(
                "Add file...",
                size=(8, 1),
                key="-fractionation_add-",
                disabled=False,
                enable_events=True,
                button_color="grey",
            ),
            sg.Combo(
                fract_paths,
                size=(58, 1),
                key="-fractionation_path-",
                disabled=False,
                enable_events=True,
                readonly=True,
                expand_x=True,
            ),
            sg.Button(
                "Remove",
                size=(8, 1),
                key="-fractionation_remove-",
                disabled=False,
                enable_events=True,
                button_color="darkred",
            ),
        ],
        [
            sg.Table(
                values=[],
                num_rows=15,
                headings=["Sample", "Condition", "Replicate", "Fraction"],
                col_widths=[34, 10, 10, 9],
                max_col_width=20,
                key="-fractionation_table-",
                auto_size_columns=False,
                vertical_scroll_only=False,
                expand_x=True,
                expand_y=True,
            )
        ],
        [
            sg.Button(
                "Remove",
                size=(8, 1),
                key="-fractionation_edit_remove-",
                disabled=False,
                enable_events=True,
                button_color="dark red",
            ),
            sg.Button(
                "Keep",
                size=(8, 1),
                key="-fractionation_edit_keep-",
                disabled=False,
                enable_events=True,
                button_color="grey",
                tooltip=" Try to keep gene names! ",
            ),
            sg.Button(
                "Set Condition",
                size=(11, 1),
                key="-fractionation_edit_condition-",
                disabled=False,
                enable_events=True,
            ),
            sg.Button(
                "Set Replicate",
                size=(11, 1),
                key="-fractionation_edit_replicate-",
                disabled=False,
                enable_events=True,
            ),
            sg.Button(
                "Set Fractions",
                size=(11, 1),
                key="-fractionation_edit_fractions-",
                disabled=False,
                enable_events=True,
            ),
            sg.Button(
                "Set Identifier",
                size=(11, 1),
                key="-fractionation_edit_identifier-",
                disabled=False,
                enable_events=True,
                button_color="grey",
                tooltip=" If possible, use protein groups! ",
            ),
        ],
        [sg.HSep()],
        [
            sg.Column(
                [
                    [
                        sg.Button(
                            "Parameters...",
                            size=(15, 1),
                            key="-fractionation_parameters-",
                            disabled=False,
                            enable_events=True,
                            button_color="black",
                        )
                    ],
                    [
                        sg.Button(
                            "Reset Fract.",
                            size=(15, 1),
                            key="-fractionation_reset-",
                            disabled=True,
                            enable_events=True,
                            button_color="dark red",
                        )
                    ],
                ],
                size=(140, 70),
            ),
            sg.Column(
                [
                    [
                        sg.Button(
                            "Process Fract.!",
                            size=(30, 1),
                            key="-fractionation_start-",
                            disabled=False,
                            enable_events=True,
                            button_color="darkgreen",
                        )
                    ],
                ],
                size=(260, 70),
            ),
            sg.Column(
                [
                    [
                        sg.Button(
                            "Plot/Export...",
                            size=(15, 1),
                            key="-fractionation_summary-",
                            disabled=True,
                            enable_events=True,
                            button_color="grey",
                        )
                    ],
                ],
                size=(140, 70),
            ),
        ],
    ]

    return sg.Tab(
        " - Fractionation - ",
        layout_fractionation,
        expand_x=True,
        expand_y=True,
    )


def create_total_proteome_tab(tp_paths: Sequence[Any]) -> sg.Tab:
    """Create the "Total Proteome" tab."""
    layout = [
        [
            sg.Button(
                "Add file...",
                size=(8, 1),
                key="-tp_add-",
                disabled=False,
                enable_events=True,
                button_color="grey",
            ),
            sg.Combo(
                tp_paths,
                size=(58, 1),
                key="-tp_path-",
                disabled=False,
                enable_events=True,
                readonly=True,
            ),
            sg.Button(
                "Remove",
                size=(8, 1),
                key="-tp_remove-",
                disabled=False,
                enable_events=True,
                button_color="darkred",
            ),
        ],
        [
            sg.Table(
                values=[],
                num_rows=15,
                headings=["Sample", "Condition"],
                col_widths=[43, 20],
                max_col_width=20,
                key="-tp_table-",
                auto_size_columns=False,
                vertical_scroll_only=False,
                expand_x=True,
                expand_y=True,
            )
        ],
        [
            sg.Button(
                "Remove",
                size=(8, 1),
                key="-tp_edit_remove-",
                disabled=False,
                enable_events=True,
                button_color="dark red",
            ),
            sg.Button(
                "Keep",
                size=(8, 1),
                key="-tp_edit_keep-",
                disabled=False,
                enable_events=True,
                button_color="grey",
                tooltip=" Try to keep gene names! ",
            ),
            sg.Button(
                "Set Condition",
                size=(11, 1),
                key="-tp_edit_condition-",
                disabled=False,
                enable_events=True,
            ),
            sg.Button(
                "Set Identifier",
                size=(11, 1),
                key="-tp_edit_identifier-",
                disabled=False,
                enable_events=True,
                button_color="grey",
                tooltip=" If possible, use protein groups! ",
            ),
        ],
        [sg.HSep()],
        [
            sg.Column(
                [
                    [
                        sg.Button(
                            "Parameters...",
                            size=(15, 1),
                            key="-tp_parameters-",
                            disabled=False,
                            enable_events=True,
                            button_color="black",
                        )
                    ],
                    [
                        sg.Button(
                            "Reset TP",
                            size=(15, 1),
                            key="-tp_reset-",
                            disabled=True,
                            enable_events=True,
                            button_color="dark red",
                        )
                    ],
                ],
                size=(140, 70),
            ),
            sg.Column(
                [
                    [
                        sg.Button(
                            "Process TP!",
                            size=(30, 1),
                            key="-tp_start-",
                            disabled=False,
                            enable_events=True,
                            button_color="darkgreen",
                        )
                    ],
                ],
                size=(260, 70),
            ),
            sg.Column(
                [
                    [
                        sg.Button(
                            "TP Summary",
                            size=(15, 1),
                            key="-tp_summary-",
                            disabled=True,
                            enable_events=True,
                            button_color="grey",
                            visible=False,
                        )
                    ],
                    [
                        sg.Button(
                            "TP Export...",
                            size=(15, 1),
                            key="-tp_export-",
                            disabled=True,
                            enable_events=True,
                            button_color="grey",
                        )
                    ],
                ],
                size=(140, 70),
            ),
        ],
    ]

    return sg.Tab(
        " - TotalProteomes - ",
        layout,
        key="-total_tab-",
        expand_x=True,
        expand_y=True,
    )


def create_data_import_frame(
    fract_paths: Sequence[Any], tp_paths: Sequence[Any]
) -> sg.Frame:
    """Create the "Data Import" frame."""
    return sg.Frame(
        layout=[
            [
                sg.TabGroup(
                    [
                        [
                            create_fractionation_tab(fract_paths),
                            create_total_proteome_tab(tp_paths),
                        ]
                    ],
                    tab_location="topleft",
                    tab_background_color="grey",
                    size=(600, 450),
                    expand_x=True,
                    expand_y=True,
                )
            ]
        ],
        title="Data Import",
        size=(620, 480),
        expand_x=True,
        expand_y=True,
    )


def create_spatial_prediction_frame() -> sg.Frame:
    """Create the "Spatial Prediction" frame."""
    static_statistics_frame = sg.Frame(
        layout=[
            [
                sg.Frame(
                    layout=[
                        [
                            sg.Button(
                                "Predict Proteome!",
                                key="-statistic_predict-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark blue",
                                expand_x=True,
                            )
                        ],
                        [sg.HSep()],
                        [
                            sg.Button(
                                "Export Prediction",
                                size=(15, 1),
                                key="-statistic_export-",
                                disabled=True,
                                enable_events=True,
                                button_color="black",
                            ),
                            sg.Button(
                                "Import Prediction",
                                size=(15, 1),
                                key="-statistic_import-",
                                disabled=False,
                                enable_events=True,
                                button_color="black",
                            ),
                        ],
                        [sg.HSep()],
                        [sg.VPush()],
                        [
                            sg.Button(
                                "Report...",
                                size=(15, 1),
                                key="-statistic_report-",
                                disabled=True,
                                enable_events=True,
                                button_color="grey",
                            ),
                            sg.Button(
                                "Reset",
                                size=(15, 1),
                                key="-statistic_reset-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark red",
                            ),
                        ],
                        [
                            sg.Button(
                                "Static Heatmap",
                                size=(15, 1),
                                key="-statistic_heatmap-",
                                disabled=True,
                                enable_events=True,
                            ),
                            sg.Button(
                                "Distribution Plots",
                                size=(15, 1),
                                key="-statistic_distribution-",
                                disabled=True,
                                enable_events=True,
                            ),
                        ],
                    ],
                    title="Proteome",
                    size=(280, 195),
                ),
                sg.Frame(
                    layout=[
                        [
                            sg.Button(
                                "Predict Lipidome!",
                                size=(32, 1),
                                key="-lipidome_predict-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark blue",
                            )
                        ],
                        [sg.VPush()],
                        [sg.HSep()],
                        [
                            sg.Button(
                                "Report...",
                                size=(15, 1),
                                key="-lipidome_report-",
                                disabled=True,
                                enable_events=True,
                                button_color="grey",
                            ),
                            sg.Button(
                                "Reset",
                                size=(15, 1),
                                key="-lipidome_reset-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark red",
                            ),
                        ],
                        [
                            sg.Button(
                                "Heatmap",
                                size=(15, 1),
                                key="-lipidome_heatmap-",
                                disabled=True,
                                enable_events=True,
                            ),
                            sg.Button(
                                "Reorganization Plot",
                                size=(15, 1),
                                key="-lipidome_reorganization-",
                                disabled=True,
                                enable_events=True,
                            ),
                        ],
                        [
                            sg.Button(
                                "Density Plot",
                                size=(15, 1),
                                key="-lipidome_density-",
                                disabled=True,
                                enable_events=True,
                            ),
                            sg.Button(
                                "Class Compositions",
                                size=(15, 1),
                                key="-lipidome_composition-",
                                disabled=True,
                                enable_events=True,
                            ),
                        ],
                    ],
                    title="Lipidome",
                    size=(290, 195),
                    visible=False,
                ),
            ],
        ],
        title="Static Statistics",
        size=(580, 230),
    )
    status_column = sg.Column(
        [
            [
                sg.Text("Prot. Fractionation:"),
                sg.Push(),
                sg.Text(
                    "none",
                    key="-status_fract-",
                    text_color="black",
                ),
            ],  # or: 'ready'
            [
                sg.Text("Total Proteome:"),
                sg.Push(),
                sg.Text(
                    "none",
                    key="-status_tp-",
                    text_color="black",
                ),
            ],  # or: 'ready'
            [
                sg.Text("Marker Proteins:"),
                sg.Push(),
                sg.Text(
                    "none",
                    key="-status_marker-",
                    text_color="black",
                ),
            ],  # or: 'ready'
            [
                sg.Text("Lipidome Fractionation:", visible=False),
                sg.Push(),
                sg.Text(
                    "none",
                    key="-status_fract_lipid-",
                    text_color="black",
                    visible=False,
                ),
            ],
            [
                sg.Text("Total Lipidome:", visible=False),
                sg.Push(),
                sg.Text(
                    "none",
                    key="-status_total_lipid-",
                    text_color="black",
                    visible=False,
                ),
            ],
        ],
        size=(240, 180),
    )
    additional_import_frame = sg.Frame(
        layout=[
            [
                sg.Text("Data Type:"),
                sg.Push(),
                sg.Combo(
                    [
                        "Lipidomics",
                    ],
                    default_value="Lipidomics",
                    size=(15, 1),
                    readonly=True,
                ),
            ],
            [sg.HSep()],
            [sg.VPush()],
            [
                sg.Button(
                    "Import/Edit Fractionation...",
                    size=(39, 1),
                )
            ],
            [
                sg.Button(
                    "Import/Edit Total Lysate...",
                    size=(39, 1),
                )
            ],
        ],
        title="Additional Import",
        size=(340, 130),
        vertical_alignment="top",
        visible=False,
    )

    return sg.Frame(
        layout=[
            [
                sg.Button(
                    "Parameters...",
                    size=(15, 1),
                    key="-classification_parameters-",
                    disabled=False,
                    enable_events=True,
                    button_color="black",
                ),
                sg.Button(
                    f"Train {app_name}!",
                    size=(25, 1),
                    key="-classification_MOP-",
                    disabled=True,
                    enable_events=True,
                    button_color="dark blue",
                ),
                sg.Button(
                    "ML Validation",
                    size=(15, 1),
                    key="-classification_validation-",
                    disabled=True,
                    enable_events=True,
                    visible=False,
                ),
                sg.Button(
                    "Reset",
                    size=(10, 1),
                    key="-classification_reset-",
                    disabled=True,
                    enable_events=True,
                    button_color="dark red",
                ),
            ],
            [sg.HSep()],
            [status_column, additional_import_frame],
            [sg.VPush()],
            [static_statistics_frame],
        ],
        title="Spatial Prediction",
        size=(600, 480),
        expand_x=True,
        expand_y=True,
    )


def create_marker_selection_frame() -> sg.Frame:
    """Create the "Marker Selection" frame."""
    tt_fract_key = (
        f"The column in the fractionation data (marked {IDENTIFIER} or {KEEP})"
        " to match the key column of the marker table."
    )
    tt_marker_key = (
        "The column in the marker table to match the key column "
        "of the fractionation data, e.g. a gene name."
    )
    tt_marker_class = (
        "The column in the marker table that contains the "
        "class (compartment) associated with the marker."
    )

    import_frame = sg.Frame(
        layout=[
            [
                sg.Listbox(
                    values=[],
                    size=(40, 4),
                    key="-marker_list-",
                    disabled=False,
                    enable_events=True,
                    horizontal_scroll=True,
                )
            ],
            [
                sg.Column(
                    layout=[
                        [
                            sg.Button(
                                "Add...",
                                size=(7, 1),
                                key="-marker_add-",
                                disabled=False,
                                enable_events=True,
                            )
                        ],
                        [
                            sg.Button(
                                "Remove",
                                size=(7, 1),
                                key="-marker_remove-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark red",
                            )
                        ],
                    ],
                ),
                sg.Column(
                    layout=[
                        [
                            sg.Text("Key column:", tooltip=tt_marker_key),
                            sg.Push(),
                            sg.Combo(
                                [],
                                size=(10, 1),
                                key="-marker_key-",
                                enable_events=True,
                                readonly=True,
                                tooltip=tt_marker_key,
                            ),
                        ],
                        [
                            sg.Text("Class column:", tooltip=tt_marker_class),
                            sg.Push(),
                            sg.Combo(
                                [],
                                size=(10, 1),
                                key="-marker_class-",
                                enable_events=True,
                                readonly=True,
                                tooltip=tt_marker_class,
                            ),
                        ],
                    ],
                    size=(220, 70),
                ),
            ],
        ],
        title="Import",
        size=(320, 190),
    )

    return sg.Frame(
        layout=[
            [
                import_frame,
                sg.Column(
                    layout=[
                        [
                            sg.Button(
                                "Load preset-list",
                                size=(32, 1),
                                key="-marker_preset-",
                                disabled=False,
                                enable_events=True,
                                button_color="dark blue",
                                visible=False,
                            )
                        ],
                        [sg.HSep()],
                        [
                            sg.Text("Fract. Key:", tooltip=tt_fract_key),
                            sg.Push(),
                            sg.Combo(
                                [IDENTIFIER],
                                key="-marker_fractkey-",
                                size=(18, 1),
                                readonly=True,
                                tooltip=tt_fract_key,
                                enable_events=True,
                            ),
                        ],
                        [sg.HSep()],
                        [
                            sg.Button(
                                "Parameters...",
                                size=(15, 1),
                                key="-marker_parameters-",
                                disabled=False,
                                enable_events=True,
                                button_color="black",
                            ),
                            sg.Button(
                                "Manage...",
                                size=(15, 1),
                                key="-marker_manage-",
                                disabled=True,
                                enable_events=True,
                                button_color="black",
                            ),
                        ],
                        [
                            sg.Button(
                                "Correlations...",
                                size=(15, 1),
                                key="-marker_test-",
                                disabled=True,
                                enable_events=True,
                                button_color="grey",
                            ),
                            sg.Button(
                                "Profiles...",
                                size=(15, 1),
                                key="-marker_profiles-",
                                disabled=True,
                                enable_events=True,
                                button_color="grey",
                            ),
                        ],
                        [sg.HSep()],
                        [
                            sg.Button(
                                "Reset",
                                size=(15, 1),
                                key="-marker_reset-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark red",
                            ),
                            sg.Button(
                                "Match!",
                                size=(15, 1),
                                key="-marker_accept-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark green",
                            ),
                        ],
                    ],
                ),
            ],
        ],
        title="Marker Selection",
        size=(620, 220),
        expand_x=True,
        expand_y=True,
    )


def create_conditional_comparison_frame() -> sg.Frame:
    """Create the "Conditional Comparison" frame."""
    global_changes_frame = sg.Frame(
        layout=[
            # [sg.Text('Statistics required!')],
            [
                sg.Button(
                    "Calculate global Changes!",
                    key="-global_run-",
                    disabled=True,
                    enable_events=True,
                    button_color="dark blue",
                    expand_x=True,
                )
            ],
            [sg.HSep()],
            [
                sg.Button(
                    "Change Heatmap",
                    size=(15, 1),
                    key="-global_heatmap-",
                    disabled=True,
                    enable_events=True,
                ),
                sg.Button(
                    "Distance Plot",
                    size=(15, 1),
                    key="-global_distance-",
                    disabled=True,
                    enable_events=True,
                ),
            ],
            [
                sg.Button(
                    "Report...",
                    size=(15, 1),
                    key="-global_report-",
                    disabled=True,
                    enable_events=True,
                    button_color="grey",
                ),
                sg.Button(
                    "Reset",
                    size=(15, 1),
                    key="-global_reset-",
                    disabled=True,
                    enable_events=True,
                    button_color="dark red",
                ),
            ],
        ],
        title="Global Changes",
        size=(290, 190),
    )
    class_centric_changes_frame = sg.Frame(
        layout=[
            [
                sg.Button(
                    "Calculate class-centric Changes!",
                    key="-class_run-",
                    disabled=True,
                    enable_events=True,
                    button_color="dark blue",
                    expand_x=True,
                )
            ],
            [sg.HSep()],
            [
                sg.Button(
                    "Class Heatmaps",
                    size=(15, 1),
                    key="-class_heatmap-",
                    disabled=True,
                    enable_events=True,
                    visible=False,
                ),
                sg.Button(
                    "Class Reorg. Plots",
                    size=(15, 1),
                    key="-class_reorganization-",
                    disabled=True,
                    enable_events=True,
                    visible=False,
                ),
            ],
            [
                sg.Button(
                    "Report...",
                    size=(15, 1),
                    key="-class_report-",
                    disabled=True,
                    enable_events=True,
                    button_color="grey",
                ),
                sg.Button(
                    "Reset",
                    size=(15, 1),
                    key="-class_reset-",
                    disabled=True,
                    enable_events=True,
                    button_color="dark red",
                ),
            ],
        ],
        title="Class-centric Changes",
        size=(290, 190),
    )
    return sg.Frame(
        layout=[[global_changes_frame, class_centric_changes_frame]],
        title="Conditional Comparison",
        size=(600, 220),
        expand_x=True,
        expand_y=True,
    )


def create_main_window(
    model: SessionModel, app_settings: AppSettings
) -> sg.Window:
    """Create the C-COMPASS main window."""

    # The main menu
    recent_files = [
        f"&{i}: {f}::recent_file_{i}"
        for i, f in enumerate(app_settings.recent_files)
    ]
    menu_def = [
        [
            "&File",
            [
                "&New",
                "&Open...",
                "&Save As...",
                "&Recent Files",
                recent_files,
                "---",
                "E&xit",
            ],
        ],
        ["&Settings", ["&Settings..."]],
        ["&Help", ["&About...", "Open &Website", "&Manual"]],
    ]

    layout = [
        [sg.Menu(menu_def, key="-menu-", tearoff=False)],
        [
            create_data_import_frame(
                fract_paths=model.fract_paths,
                tp_paths=model.tp_paths,
            ),
            create_spatial_prediction_frame(),
        ],
        [
            create_marker_selection_frame(),
            create_conditional_comparison_frame(),
        ],
    ]

    main_window = sg.Window(
        app_name,
        layout,
        size=(1260, 720),
        resizable=True,
    )
    return main_window


class MainController:
    """The main controller for the C-COMPASS application."""

    def __init__(self, model: SessionModel):
        self.model = model
        self.app_settings = AppSettings.load()
        self.main_window = create_main_window(
            model=model, app_settings=self.app_settings
        )

    def run(self):
        """Run the C-COMPASS application."""

        # The event loop
        while True:
            event, values = self.main_window.read()

            if event == sg.WIN_CLOSED or event == "Exit":
                break

            if event == "-fractionation_add-":
                fract_add(self.main_window, model=self.model)
            elif event == "-fractionation_remove-":
                if values["-fractionation_path-"]:
                    fract_remove_file(
                        values,
                        self.main_window,
                        model=self.model,
                    )
            elif event == "-fractionation_path-":
                fract_refreshtable(
                    self.main_window,
                    self.model.fract_input[
                        values["-fractionation_path-"]
                    ].table,
                )
            elif event == "-fractionation_edit_remove-":
                if values["-fractionation_table-"]:
                    fract_remove_row(
                        values, self.main_window, self.model.fract_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_keep-":
                if values["-fractionation_table-"]:
                    fract_set_keep(
                        values, self.main_window, self.model.fract_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_condition-":
                if values["-fractionation_table-"]:
                    fract_define_condition(
                        values, self.main_window, self.model.fract_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_replicate-":
                if values["-fractionation_table-"]:
                    fract_define_replicate(
                        values, self.main_window, self.model.fract_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_fractions-":
                if values["-fractionation_table-"]:
                    fract_handle_set_fraction(
                        values, self.main_window, self.model.fract_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_identifier-":
                fract_handle_set_identifier(
                    values,
                    self.main_window,
                    self.model.fract_input,
                )
            elif event == "-fractionation_parameters-":
                from .fractionation_parameters_dialog import show_dialog

                self.model.fract_preparams = show_dialog(
                    self.model.fract_preparams
                )
            elif event == "-fractionation_reset-":
                self._handle_reset_fract_data()
            elif event == "-fractionation_start-":
                self._handle_process_fract_data()
            elif event == "-fractionation_summary-":
                RP.RP_gradient_heatmap(self.model.fract_data)

            elif event == "-tp_add-":
                tp_add_dataset(
                    self.main_window,
                    self.model.tp_input,
                )
            elif event == "-tp_remove-":
                tp_remove_dataset(
                    values,
                    self.main_window,
                    self.model.tp_input,
                )
            elif event == "-tp_path-":
                tp_refreshtable(
                    self.main_window,
                    self.model.tp_input[values["-tp_path-"]].table,
                )
            elif event == "-tp_edit_remove-":
                if values["-tp_table-"]:
                    tp_remove_row(
                        values, self.main_window, self.model.tp_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-tp_edit_keep-":
                if values["-tp_table-"]:
                    tp_set_keep(values, self.main_window, self.model.tp_input)
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-tp_edit_condition-":
                if values["-tp_table-"]:
                    tp_set_condition(
                        values, self.main_window, self.model.tp_input
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-tp_edit_identifier-":
                tp_set_identifier(
                    values,
                    self.main_window,
                    self.model.tp_input,
                )
            elif event == "-tp_parameters-":
                from .total_proteome_parameters_dialog import show_dialog

                self.model.tp_preparams = show_dialog(self.model.tp_preparams)
            elif event == "-tp_reset-":
                self._handle_reset_total_proteome()
            elif event == "-tp_start-":
                self._handle_process_total_proteome()
            elif event == "-tp_export-":
                export_folder = sg.popup_get_folder("Export Folder")
                if export_folder:
                    experiment = simpledialog.askstring(
                        "Export", "Experiment Name: "
                    )

                    tp_export(
                        export_folder,
                        experiment,
                        self.model.tp_data,
                        self.model.tp_info,
                    )

            elif event == "-marker_add-":
                marker_add(self.main_window, self.model.marker_sets)
            elif event == "-marker_remove-":
                marker_remove(self.main_window, values, self.model.marker_sets)
            elif event == "-marker_list-":
                refresh_markercols(
                    self.main_window, values, self.model.marker_sets
                )
            elif event == "-marker_key-":
                marker_setkey(values, self.model.marker_sets)
            elif event == "-marker_class-":
                self.model.marker_conv = marker_setclass(
                    values, self.model.marker_sets
                )
            elif event == "-marker_parameters-":
                from .marker_parameters_dialog import show_dialog

                self.model.marker_params = show_dialog(
                    self.model.marker_params
                )
            elif event == "-marker_manage-":
                from .class_manager_dialog import show_class_manager_dialog

                if check_markers(self.model.marker_sets):
                    self.model.marker_conv = show_class_manager_dialog(
                        self.model.marker_conv
                    )
                else:
                    messagebox.showerror(
                        "Error", "Please define key and class column."
                    )
            elif event == "-marker_fractkey-":
                self.model.marker_fractkey = values["-marker_fractkey-"]
            elif event == "-marker_test-":
                self._open_marker_correlations()
            elif event == "-marker_profiles-":
                self._open_marker_profiles()
            elif event == "-marker_accept-":
                self._handle_match_markers(values)
            elif event == "-marker_reset-":
                self.model.reset_marker()
            elif event == "-classification_parameters-":
                from .training_parameters_dialog import show_dialog

                self.model.NN_params = show_dialog(self.model.NN_params)

            elif event == "-classification_MOP-":
                self._handle_training()
            elif event == "-classification_reset-":
                self.model.reset_classification()

            elif event == "-statistic_predict-":
                with wait_cursor(self.main_window):
                    self.model.results = MOA.stats_proteome(
                        self.model.learning_xyz,
                        self.model.fract_data,
                        self.model.fract_conditions,
                        self.model.NN_params.reliability,
                    )

            elif event == "-statistic_export-":
                filename = sg.popup_get_file(
                    "Export Statistics",
                    no_window=True,
                    file_types=(("Pickle", "*.pkl"),),
                    save_as=True,
                )
                if filename:
                    with open(filename, "wb") as file:
                        pickle.dump(self.model.results, file)

            elif event == "-statistic_import-":
                self._handle_import_prediction()

            elif event == "-statistic_report-":
                self._handle_statistics_report()
            elif event == "-global_report-":
                self._handle_global_changes_report()

            elif event == "-class_report-":
                self._handle_class_changes_report()
            elif event == "-statistic_reset-":
                self.model.reset_static_statistics()

            elif event == "-statistic_heatmap-":
                RP.RP_stats_heatmap(self.model.results)

            elif event == "-statistic_distribution-":
                RP.RP_stats_distribution(self.model.results)

            elif event == "-global_heatmap-":
                RP.RP_global_heatmap(self.model.comparison)

            elif event == "-global_distance-":
                RP.RP_global_distance(self.model.comparison)

            elif event == "-class_heatmap-":
                RP.RP_class_heatmap(self.model.results)

            elif event == "-global_run-":
                with wait_cursor(self.main_window):
                    self.model.comparison = MOA.global_comparisons(
                        results=self.model.results,
                        max_processes=self.app_settings.max_processes,
                    )
            elif event == "-global_reset-":
                self.model.reset_global_changes()

            elif event == "-class_run-":
                with wait_cursor(self.main_window):
                    MOA.class_comparisons(
                        self.model.tp_data,
                        self.model.results,
                        self.model.comparison,
                    )
            elif event == "-class_reset-":
                self.model.reset_class_centric_changes()
            elif event == "-export_statistics-":
                self._handle_export_statistics()
            elif event == "-export_comparison-":
                self._handle_export_comparison()
            elif event == "Save As...":
                self._handle_session_save()
            elif event == "Open...":
                self._handle_session_open()
            elif event == "New":
                self._handle_session_new()
            elif "::recent_file_" in event:
                self._handle_session_open_recent(event)
            elif event == "About...":
                from .about_dialog import show_about_dialog

                show_about_dialog()
            elif event == "Open Website":
                import webbrowser

                webbrowser.open(repository_url)
            elif event == "Manual":
                import webbrowser

                webbrowser.open(readthedocs_url)
            elif event == "Settings...":
                from .settings_dialog import show_settings_dialog

                show_settings_dialog(self.app_settings)
            else:
                logger.debug(f"Unhandled event: {event}")

            refresh_window(self.main_window, self.model.status)

        self.main_window.close()

    def _open_marker_correlations(self):
        """Open the marker correlation window."""
        if check_markers(self.model.marker_sets):
            from .marker_correlation_dialog import (
                show_marker_correlation_dialog,
            )

            try:
                self.model.marker_list = create_markerlist(
                    self.model.marker_sets,
                    self.model.marker_conv,
                    **self.model.marker_params,
                )
                show_marker_correlation_dialog(
                    self.model.fract_data,
                    self.model.fract_info,
                    self.model.marker_list,
                    self.model.marker_fractkey,
                )
            except Exception:
                logger.exception("Error")
                messagebox.showerror(
                    "Error",
                    "Something is wrong with your marker list.",
                )
        else:
            messagebox.showerror(
                "Error", "Please define key and class column."
            )

    def _open_marker_profiles(self):
        """Open the marker profiles window."""
        if not check_markers(self.model.marker_sets):
            messagebox.showerror(
                "Error", "Please define key and class column."
            )
            return

        from .marker_profiles_dialog import show_marker_profiles_dialog

        try:
            self.model.marker_list = create_markerlist(
                self.model.marker_sets,
                self.model.marker_conv,
                **self.model.marker_params,
            )
            show_marker_profiles_dialog(
                self.model.fract_data,
                self.model.fract_info,
                self.model.marker_list,
                self.model.marker_fractkey,
            )
        except Exception:
            logger.exception("Error")
            messagebox.showerror(
                "Error",
                "Something is wrong with your marker list.",
            )

    def _handle_match_markers(self, values: dict):
        if not values["-marker_list-"]:
            messagebox.showerror(
                "Error", "Please import at least one Marker List!"
            )
            return

        if not self.model.fract_data["class"]:
            messagebox.showerror(
                "Error", "Please import Fractionation Data first!"
            )
            return

        if (
            values["-marker_fractkey-"] == ""
            or values["-marker_class-"] == ""
            or values["-marker_key-"] == ""
        ):
            messagebox.showerror(
                "Error", "Please select key and class columns!"
            )
            return

        try:
            with wait_cursor(self.main_window):
                self.model.marker_list = create_markerlist(
                    self.model.marker_sets,
                    self.model.marker_conv,
                    **self.model.marker_params,
                )
                logger.info("Marker list created")
                (
                    self.model.fract_marker,
                    self.model.fract_marker_vis,
                    self.model.fract_test,
                ) = create_marker_profiles(
                    self.model.fract_data,
                    self.model.marker_fractkey,
                    self.model.fract_info,
                    self.model.marker_list,
                )
                logger.info("Marker profiles created")
                self.model.fract_full = create_fullprofiles(
                    self.model.fract_marker, self.model.fract_test
                )
                logger.info("Full profiles created")
        except Exception:
            logger.exception("Error matching markers")
            messagebox.showerror("Error", "Incompatible Fractionation Key!")
            self.model.reset_marker()

    def _handle_training(self):
        """Handle click on "Train C-COMPASS!" button."""
        with wait_cursor(self.main_window):
            from .MOP import MOP_exec

            self.model.learning_xyz = MOP_exec(
                self.model.fract_full,
                self.model.fract_marker,
                self.model.fract_test,
                self.model.NN_params,
                max_processes=self.app_settings.max_processes,
            )

    def _handle_import_prediction(self):
        filename = sg.popup_get_file(
            "Import Statistics",
            no_window=True,
            file_types=(("Pickle", "*.pkl"),),
        )
        if not filename:
            return

        with open(filename, "rb") as file:
            results_new = pickle.load(file)
        try:
            for condition in results_new:
                if condition in self.model.results:
                    messagebox.showerror(
                        "Error",
                        "There are already statistics for "
                        f"{condition} in your current session.",
                    )
                else:
                    self.model.results[condition] = copy.deepcopy(
                        results_new[condition]
                    )
        except Exception:
            logger.exception("Error importing prediction")
            messagebox.showerror("Error", "Incompatible file type!")

    def _handle_statistics_report(self):
        export_folder = sg.popup_get_folder("Statistics Report")
        if not export_folder:
            return

        Path(export_folder).mkdir(parents=True, exist_ok=True)

        for condition, result in self.model.results.items():
            fname = Path(export_folder, f"CCMPS_statistics_{condition}.xlsx")
            selected_columns = [
                col for col in result.metrics.columns if col.startswith("fCC_")
            ] + ["SVM_winner", "fNN_winner", "marker"]
            df_out = result.metrics[selected_columns]
            df_out.columns = [
                col.replace("fCC_", "CC_ClassContribution_")
                if col.startswith("fCC_")
                else "C-CMPS_MainClass"
                if col == "fNN_winner"
                else col
                for col in df_out.columns
            ]
            df_out.to_excel(fname, index=True)

    def _handle_global_changes_report(self):
        export_folder = sg.popup_get_folder("Global Changes Report")
        if not export_folder:
            return

        write_global_changes_reports(self.model.comparison, export_folder)

    def _handle_class_changes_report(self):
        export_folder = sg.popup_get_folder("Class-centric Changes Report")
        if not export_folder:
            return

        write_class_changes_reports(self.model, export_folder)

    def _handle_export_statistics(self):
        export_folder = sg.popup_get_folder("Export Statistics")
        if not export_folder:
            return

        write_statistics_reports(self.model, export_folder)

    def _handle_export_comparison(self):
        export_folder = sg.popup_get_folder("Export Comparison")
        if not export_folder:
            return

        write_comparison_reports(self.model, export_folder)

    def _handle_session_new(self):
        sure = sg.popup_yes_no(
            "Are you sure to close the session and start a new one?"
        )
        if sure != "Yes":
            return

        self.model.reset()

        fract_clearinput(self.main_window)
        tp_clearinput(self.main_window)

        self.main_window["-marker_list-"].update(values=[])
        self.main_window["-marker_key-"].update(
            values=[], size=self.main_window["-marker_key-"].Size
        )
        self.main_window["-marker_class-"].update(
            values=[], size=self.main_window["-marker_class-"].Size
        )
        self.main_window["-marker_fractkey-"].update(
            values=[IDENTIFIER] + list(self.model.fract_info)
        )

    def _handle_session_open(self):
        """'Open session' was clicked."""
        filename = sg.popup_get_file(
            "Open Session",
            initial_folder=str(self.app_settings.last_session_dir),
            no_window=True,
            file_types=(
                ("Numpy", "*.npy"),
                ("C-COMPASS zip", "*.ccompass"),
            ),
        )
        if not filename:
            return

        self.open_session(filename)

    def open_session(self, filename: str):
        try:
            with wait_cursor(self.main_window):
                session_open(
                    self.main_window,
                    filename,
                    model=self.model,
                )
        except Exception as e:
            logger.exception("Error opening session")
            messagebox.showerror(
                "Error",
                "An error occurred while opening the session:\n\n" + str(e),
            )

        self.main_window["-marker_fractkey-"].update(
            values=[IDENTIFIER] + list(self.model.fract_info),
            value=self.model.marker_fractkey,
        )

        self.app_settings.add_recent_file(filename)
        self.app_settings.last_session_dir = Path(filename).parent
        self.app_settings.save()

        self._update_recent_files()

    def _handle_session_save(self):
        """'Save session' was clicked."""
        filename = sg.popup_get_file(
            "Save Session",
            no_window=True,
            file_types=(
                ("Numpy", "*.npy"),
                ("C-COMPASS zip", "*.ccompass"),
            ),
            save_as=True,
            initial_folder=str(self.app_settings.last_session_dir),
        )
        if not filename:
            return

        self.app_settings.add_recent_file(filename)
        self.app_settings.last_session_dir = Path(filename).parent
        self.app_settings.save()

        with wait_cursor(self.main_window):
            if str(filename).endswith(".ccompass"):
                self.model.to_zip(filename)
            else:
                self.model.to_numpy(filename)

        self._update_recent_files()

    def _handle_session_open_recent(self, event: str):
        """Open a session from the recent files list."""
        file_idx = int(event.split("::")[-1].removeprefix("recent_file_"))
        filename = self.app_settings.recent_files[file_idx]
        self.open_session(filename)

    def _update_recent_files(self):
        """Update the recent files menu."""
        menu = self.main_window["-menu-"]
        menu_def = menu.MenuDefinition
        recent_files = [
            f"&{i}: {f}::recent_file_{i}"
            for i, f in enumerate(self.app_settings.recent_files)
        ]

        # find "Recent Files" entry
        assert menu_def[0][1][3] == "Recent Files"
        menu_def[0][1][4] = recent_files
        menu.update(menu_definition=menu_def)

    def _handle_process_total_proteome(self):
        """Button-click "process total proteome data"."""
        if not self.model.tp_input:
            messagebox.showerror("No dataset!", "Please import a TP dataset.")
            return

        from .TPP import total_proteome_processing_dialog

        (
            self.model.tp_data,
            self.model.tp_info,
            self.model.tp_icorr,
        ) = total_proteome_processing_dialog(
            self.model.tp_input,
            self.model.tp_preparams,
            self.model.tp_data,
            self.model.tp_info,
            self.model.tp_icorr,
        )

    def _handle_process_fract_data(self):
        """Button-click "process fractionation data"."""
        if not self.model.fract_input:
            messagebox.showerror(
                "No dataset!", "Please import a fractionation dataset."
            )
            return

        from .FDP import FDP_exec

        (
            self.model.fract_data,
            self.model.fract_std,
            self.model.fract_info,
            self.model.fract_conditions,
        ) = FDP_exec(
            self.model.fract_input,
            self.model.fract_preparams,
            self.model.fract_data,
            self.model.fract_std,
            self.model.fract_info,
            self.model.fract_conditions,
        )

        self.main_window["-marker_fractkey-"].update(
            values=[IDENTIFIER] + list(self.model.fract_info)
        )

    def _handle_reset_fract_data(self):
        """Button-click "reset fractionation data"."""
        sure = sg.popup_yes_no(
            "Reset Fractionation Pre-Processing? "
            "You have to run it again to use your data."
        )
        if sure != "Yes":
            return

        self.model.reset_fractionation()
        fract_buttons(self.main_window, False)

        self.main_window["-marker_fractkey-"].update(
            values=[IDENTIFIER], value=""
        )

    def _handle_reset_total_proteome(self):
        """Button-click "reset total proteome data"."""
        sure = sg.popup_yes_no(
            "Reset TotalProteome Pre-Processing? "
            "You have to run it again to use your data."
        )
        if sure != "Yes":
            return

        self.model.reset_tp()
        if self.model.status.comparison_class:
            self.model.reset_class_centric_changes()


def fract_refreshtable(window: sg.Window, table: list):
    window["-fractionation_table-"].update(values=table)


def tp_refreshtable(window: sg.Window, table: list):
    window["-tp_table-"].update(values=table)


def fract_modifytable(
    title: str,
    prompt: str,
    values: dict,
    fract_input: dict[str, FractDataset],
    pos: int,
    q: int,
    ask: Literal["integer", "string"],
):
    """Show dialog for updating values in the fractionation table."""
    if not values["-fractionation_table-"]:
        messagebox.showerror("Error", "Select (a) sample(s).")
        return

    path = values["-fractionation_path-"]
    table = fract_input[path].table
    if ask == "integer":
        value = simpledialog.askinteger(title, prompt)
        if value:
            p = 0
            for i in values["-fractionation_table-"]:
                table[i][pos] = value + p
                p = p + q
    elif ask == "string":
        value = simpledialog.askstring(title, prompt)
        if value:
            for i in values["-fractionation_table-"]:
                table[i][pos] = value
    else:
        raise ValueError(f"Invalid ask value: {ask}")


def fract_buttons(window: sg.Window, status: bool) -> None:
    active = [
        "-fractionation_add-",
        "-fractionation_remove-",
        "-fractionation_edit_remove-",
        "-fractionation_edit_keep-",
        "-fractionation_edit_condition-",
        "-fractionation_edit_replicate-",
        "-fractionation_edit_fractions-",
        "-fractionation_edit_identifier-",
        "-fractionation_parameters-",
        "-fractionation_start-",
    ]
    inactive = [
        "-fractionation_reset-",
        "-fractionation_summary-",
    ]
    for button in active:
        window[button].update(disabled=status)
    for button in inactive:
        window[button].update(disabled=not status)


def tp_buttons(window: sg.Window, status: bool) -> None:
    active = [
        "-tp_add-",
        "-tp_remove-",
        "-tp_edit_remove-",
        "-tp_edit_keep-",
        "-tp_edit_condition-",
        "-tp_edit_identifier-",
        "-tp_start-",
        "-tp_parameters-",
    ]
    inactive = ["-tp_reset-", "-tp_summary-", "-tp_export-"]
    for button in active:
        window[button].update(disabled=status)
    for button in inactive:
        window[button].update(disabled=not status)


def fract_clearinput(window):
    window["-fractionation_path-"].update(
        values=[], size=window["-fractionation_path-"].Size
    )
    window["-fractionation_table-"].update(values=[])


def tp_clearinput(window):
    window["-tp_path-"].update(values=[], size=window["-tp_path-"].Size)
    window["-tp_table-"].update(values=[])


def fract_add(
    window,
    model: SessionModel,
):
    """Add a fractionation dataset."""
    filename = sg.popup_get_file(
        "Chose dataset",
        no_window=True,
        file_types=(
            ("Tab Separated Values", "*.tsv"),
            ("Text (tab delimited)", "*.txt"),
        ),
    )
    if not filename:
        return

    try:
        df = read_fract_table(filename)
    except Exception:
        logger.exception(
            f"Error reading fractionation dataset from {filename}"
        )
        messagebox.showerror("Error", "Invalid file format.")
        return

    table = [[name, "", "", ""] for name in df.columns]
    model.fract_input[filename] = FractDataset(df=df, table=table)

    window["-fractionation_path-"].update(
        values=model.fract_paths, value=filename
    )

    fract_refreshtable(window, table)


def fract_remove_file(values, window, model: SessionModel):
    """Remove a fractionation dataset."""
    sure = sg.popup_yes_no("Remove data from list?")
    if sure != "Yes":
        return

    filepath = values["-fractionation_path-"]
    del model.fract_input[filepath]

    filepath = next(iter(model.fract_input)) if model.fract_input else []
    fract_refreshtable(
        window, model.fract_input[filepath].table if filepath else []
    )
    window["-fractionation_path-"].update(
        values=model.fract_paths, value=filepath
    )


def fract_remove_row(values, window, fract_input: dict[str, FractDataset]):
    """Remove a row from the fractionation table."""
    path = values["-fractionation_path-"]
    selected = values["-fractionation_table-"]
    table = fract_input[path].table
    for index in sorted(selected, reverse=True):
        del table[index]
    window["-fractionation_table-"].update(values=table)


def fract_set_keep(values, window, fract_input: dict[str, FractDataset]):
    """Set the selected rows to "keep"."""
    path = values["-fractionation_path-"]
    table = fract_input[path].table
    for pos in values["-fractionation_table-"]:
        table[pos][1] = KEEP
        table[pos][2] = NA
        table[pos][3] = NA
    window["-fractionation_table-"].update(values=table)


def fract_define_condition(
    values, window, fract_input: dict[str, FractDataset]
):
    """Define the condition for the selected rows."""
    fract_modifytable(
        "Set Condition",
        "Condition Name:",
        values,
        fract_input,
        1,
        0,
        "string",
    )
    window["-fractionation_table-"].update(
        values=fract_input[values["-fractionation_path-"]].table
    )


def fract_define_replicate(
    values, window, fract_input: dict[str, FractDataset]
):
    """Define the replicate for the selected rows."""
    fract_modifytable(
        "Set Replicate",
        "Replicate Number:",
        values,
        fract_input,
        2,
        0,
        "integer",
    )
    window["-fractionation_table-"].update(
        values=fract_input[values["-fractionation_path-"]].table
    )


def fract_handle_set_fraction(
    values, window, fract_input: dict[str, FractDataset]
):
    """Set the fraction for the selected rows."""
    fract_modifytable(
        "Set Fractions",
        "FIRST Fraction Number:",
        values,
        fract_input,
        3,
        1,
        "integer",
    )
    window["-fractionation_table-"].update(
        values=fract_input[values["-fractionation_path-"]].table
    )


def fract_handle_set_identifier(
    values,
    window: sg.Window,
    fract_input: dict[str, FractDataset],
) -> None:
    """Set the identifier for the selected rows."""
    if not (pos := values["-fractionation_table-"]):
        messagebox.showerror("Error", "No sample selected.")
        return

    if len(pos) > 1:
        messagebox.showerror("Error", "Please set only one Identifier!")
        return

    path = values["-fractionation_path-"]
    dataset = fract_input[path]
    dataset.id_col = pos[0]

    window["-fractionation_table-"].update(values=dataset.table)


def tp_add_dataset(
    window: sg.Window,
    tp_input: dict[str, TotalProtDataset],
):
    """Add a total proteome dataset."""
    filename = sg.popup_get_file(
        "Chose dataset",
        no_window=True,
        file_types=(
            ("Tab Separated Values", "*.tsv"),
            ("Text (tab delimited)", "*.txt"),
        ),
    )
    if not filename:
        return

    # read file
    try:
        df = read_tp_table(filename)
    except Exception:
        logger.exception(
            f"Error reading total proteome dataset from {filename}"
        )
        messagebox.showerror("Error", "Invalid file format.")
        return

    table = [[name, ""] for name in df.columns]
    tp_input[filename] = TotalProtDataset(df=df, table=table)

    tp_refreshtable(window, table)
    window["-tp_path-"].update(values=list(tp_input), value=filename)


def tp_remove_dataset(
    values,
    window: sg.Window,
    tp_input: dict[str, TotalProtDataset],
):
    """Remove a total proteome dataset."""
    if not (selected_file := values["-tp_path-"]):
        return

    sure = sg.popup_yes_no("Remove data from list?")
    if sure != "Yes":
        return

    del tp_input[selected_file]

    files = list(tp_input)
    selected_file_new = files[0] if files else []
    window["-tp_path-"].update(
        values=files, value=selected_file_new, size=window["-tp_path-"].Size
    )

    cur_table = tp_input[selected_file_new].table if selected_file_new else []
    tp_refreshtable(window, cur_table)


def tp_remove_row(values, window, tp_input):
    """Remove a row from the total proteome data table."""
    path = values["-tp_path-"]
    selected = values["-tp_table-"]
    table = tp_input[path].table
    for index in sorted(selected, reverse=True):
        del table[index]
    window["-tp_table-"].update(values=table)


def tp_set_keep(values, window, tp_input):
    """Set a total proteome row to [KEEP]."""
    path = values["-tp_path-"]
    table = tp_input[path].table
    for pos in values["-tp_table-"]:
        table[pos][1] = KEEP
    window["-tp_table-"].update(values=table)


def tp_set_condition(values, window, tp_input):
    """Set the condition of the selected total proteome rows."""
    if not values["-tp_table-"]:
        messagebox.showerror("Error", "Select (a) sample(s).")
        return

    path = values["-tp_path-"]
    table = tp_input[path].table
    value = simpledialog.askstring("Set Condition", "Condition Name")
    if value:
        for i in values["-tp_table-"]:
            table[i][1] = value

    window["-tp_table-"].update(values=table)


def tp_set_identifier(
    values,
    window: sg.Window,
    tp_input: dict[str, TotalProtDataset],
) -> None:
    """Set the identifier of the selected total proteome row."""
    if not (pos := values["-tp_table-"]):
        messagebox.showerror("Error", "No sample selected.")
        return

    if len(pos) > 1:
        messagebox.showerror("Error", "Please set only one Identifier!")
        return

    path = values["-tp_path-"]
    tp_input[path].id_col = pos[0]
    window["-tp_table-"].update(values=tp_input[path].table)


def tp_export(export_folder: str | Path, experiment: str, tp_data, tp_info):
    """Export total proteome data."""
    export_folder = Path(export_folder)
    export_folder.mkdir(parents=True, exist_ok=True)
    now = datetime.now()
    time = now.strftime("%Y%m%d%H%M%S")

    export_full = pd.DataFrame()
    for condition in tp_data:
        path = export_folder / f"{time}_{experiment}_{condition}.txt"
        tp_data[condition].to_csv(
            path,
            header=True,
            index=True,
            index_label="Identifier",
            sep="\t",
            mode="a",
        )

        export_full = pd.merge(
            export_full,
            tp_data[condition],
            left_index=True,
            right_index=True,
            how="outer",
        )

    # add [KEEP] columns
    for info in tp_info:
        export_full = pd.merge(
            export_full,
            tp_info[info],
            left_index=True,
            right_index=True,
            how="left",
        )

    path = export_folder / f"{time}_{experiment}_combined.txt"
    export_full.to_csv(
        path,
        header=True,
        index=True,
        index_label="Identifier",
        sep="\t",
        mode="a",
    )


def check_markers(marker_sets: dict[str, MarkerSet]) -> bool:
    """Check if identifier and class columns are set for all marker sets.

    :return: True if all marker sets have identifier and class columns set.
    """
    if not marker_sets:
        return False

    for marker_set in marker_sets.values():
        if marker_set.identifier_col == NA or marker_set.class_col == NA:
            return False

    return True


def refresh_markertable(window: sg.Window, marker_sets: dict[str, MarkerSet]):
    """Update the marker table according to the marker sets."""
    file_list = list(marker_sets.keys())
    window["-marker_list-"].update(values=file_list)

    if file_list:
        window["-marker_list-"].update(set_to_index=0)
        cur_marker_set = marker_sets[file_list[0]]
        column_ids = cur_marker_set.df.columns.tolist()
        window["-marker_key-"].update(
            values=column_ids,
            value=cur_marker_set.identifier_col,
            size=window["-marker_key-"].Size,
        )
        window["-marker_class-"].update(
            values=column_ids,
            value=cur_marker_set.class_col,
            size=window["-marker_class-"].Size,
        )


def refresh_markercols(window, values, marker_sets: dict[str, MarkerSet]):
    try:
        marker_filename = values["-marker_list-"][0]
        marker_set = marker_sets[marker_filename]
        marker_set_col_ids = marker_set.df.columns.tolist()
        window["-marker_key-"].update(
            values=marker_set_col_ids,
            value=marker_set.identifier_col,
            size=window["-marker_key-"].Size,
        )
        window["-marker_class-"].update(
            values=marker_set_col_ids,
            value=marker_set.class_col,
            size=window["-marker_class-"].Size,
        )
    except Exception:
        logger.exception("Error")

        window["-marker_key-"].update(
            values=[], value=NA, size=window["-marker_key-"].Size
        )
        window["-marker_class-"].update(
            values=[], value=NA, size=window["-marker_class-"].Size
        )


def marker_add(window, marker_sets: dict[str, MarkerSet]):
    """Show marker list selection dialog
    and add the selected list to the marker sets."""
    filename = sg.popup_get_file(
        "Select a new Marker List!",
        no_window=True,
        file_types=(
            ("Tab delimited Text", "*.txt"),
            ("Tab Separated Values", "*.tsv"),
        ),
    )
    if not filename:
        return

    try:
        df = read_marker_table(filename)
    except Exception:
        logger.exception(f"Error reading marker list from {filename}")
        messagebox.showerror("Error", "Invalid file format.")
        return

    marker_sets[filename] = MarkerSet(df=df)
    refresh_markertable(window, marker_sets)


def marker_remove(window, values, marker_sets):
    marker_filename = values["-marker_list-"][0]
    del marker_sets[marker_filename]
    refresh_markertable(window, marker_sets)
    if not len(marker_sets) > 0:
        window["-marker_test-"].update(disabled=True)
        window["-marker_profiles-"].update(disabled=True)
        window["-marker_remove-"].update(disabled=True)


def marker_setkey(values, marker_sets: dict[str, MarkerSet]):
    """Set the identifier column for the selected marker list."""
    marker_filename = values["-marker_list-"][0]
    marker_sets[marker_filename].identifier_col = values["-marker_key-"]


def marker_setclass(values, marker_sets):
    """Set the class column for the selected marker list."""
    marker_filename = values["-marker_list-"][0]
    marker_sets[marker_filename].class_col = values["-marker_class-"]
    return create_identity_conversion(marker_sets.values())


def session_open(window: sg.Window, filename: str, model: SessionModel):
    """Read session data from file and update the window."""
    # Update session data
    if str(filename).endswith(".ccompass"):
        tmp_session = SessionModel.from_zip(filename)
    else:
        tmp_session = SessionModel.from_numpy(filename)
    model.reset(tmp_session)

    # update GUI
    if fract_filepaths := model.fract_paths:
        fract_refreshtable(window, model.fract_input[fract_filepaths[0]].table)
        window["-fractionation_path-"].update(
            values=fract_filepaths, value=fract_filepaths[0]
        )
    else:
        fract_refreshtable(window, [])
        window["-fractionation_path-"].update(values=[], value="")

    fract_buttons(window, bool(model.fract_data["class"]))

    if tp_filepaths := model.tp_paths:
        tp_refreshtable(window, model.tp_input[tp_filepaths[0]].table)
        window["-tp_path-"].update(values=tp_filepaths, value=tp_filepaths[0])
    else:
        tp_refreshtable(window, [])
        window["-tp_path-"].update(values=[], value="")

    tp_buttons(window, bool(model.tp_data))

    if model.marker_sets:
        refresh_markertable(window, model.marker_sets)

        event, values = window.read(timeout=50)
        refresh_markercols(window, values, model.marker_sets)


def enable_markersettings(window: sg.Window, status: SessionStatusModel):
    """Enable / disable marker settings based on the current session status."""
    for element in ["-marker_remove-", "-marker_manage-", "-marker_accept-"]:
        if status.marker_matched:
            window[element].update(disabled=True)
        else:
            window[element].update(disabled=not status.marker_file)

    for element in ["-marker_test-", "-marker_profiles-"]:
        window[element].update(disabled=not status.marker_file)

    for element in ["-marker_reset-"]:
        window[element].update(disabled=not status.marker_matched)

    for element in ["-marker_add-", "-marker_parameters-", "-marker_preset-"]:
        window[element].update(disabled=status.marker_matched)

    if status.marker_matched:
        window["-status_marker-"].update("ready", text_color="dark green")
    else:
        window["-status_marker-"].update("none", text_color="black")


def refresh_window(window: sg.Window, status: SessionStatusModel):
    """Update the window based on the current session status."""
    fract_buttons(window, status.fractionation_data)
    tp_buttons(window, status.tp_data)
    enable_markersettings(window, status)

    for element in ["-statistic_import-"]:
        window[element].update(disabled=status.comparison_global)

    if status.fractionation_data:
        window["-status_fract-"].update("ready", text_color="dark green")
    else:
        window["-status_fract-"].update("none", text_color="black")

    if status.tp_data:
        window["-status_tp-"].update("ready", text_color="dark green")
    else:
        window["-status_tp-"].update("none", text_color="black")

    if status.lipidome_data:
        window["-status_fract_lipid-"].update("ready", text_color="dark green")
    else:
        window["-status_fract_lipid-"].update("none", text_color="black")

    if status.lipidome_total:
        window["-status_total_lipid-"].update("ready", text_color="dark green")
    else:
        window["-status_total_lipid-"].update("none", text_color="black")

    for element in ["-classification_MOP-"]:
        if status.fractionation_data and status.marker_matched:
            window[element].update(disabled=status.training)
        else:
            window[element].update(disabled=True)

    for element in ["-classification_validation-", "-classification_reset-"]:
        window[element].update(disabled=not status.training)

    if status.training:
        for element in ["-statistic_predict-"]:
            window[element].update(disabled=status.proteome_prediction)

        for element in [
            "-statistic_export-",
            "-statistic_report-",
            "-statistic_reset-",
            "-statistic_heatmap-",
            "-statistic_distribution-",
        ]:
            window[element].update(disabled=not status.proteome_prediction)

        if status.proteome_prediction:
            for element in ["-global_run-"]:
                window[element].update(disabled=status.comparison_global)

            for element in [
                "-global_heatmap-",
                "-global_distance-",
                "-global_report-",
                "-global_reset-",
            ]:
                window[element].update(disabled=not status.comparison_global)

            if status.comparison_global and status.tp_data:
                for element in ["-class_run-"]:
                    window[element].update(disabled=status.comparison_class)

                for element in [
                    "-class_heatmap-",
                    "-class_reorganization-",
                    "-class_report-",
                    "-class_reset-",
                ]:
                    window[element].update(
                        disabled=not status.comparison_class
                    )
            else:
                for element in [
                    "-class_run-",
                    "-class_heatmap-",
                    "-class_reorganization-",
                    "-class_report-",
                    "-class_reset-",
                ]:
                    window[element].update(disabled=True)

                if status.lipidome_data:
                    for element in ["-lipidome_predict-"]:
                        window[element].update(
                            disabled=status.lipidome_prediction
                        )

                    for element in [
                        "-lipidome_report-",
                        "-lipidome_reset-",
                        "-lipidome_heatmap-",
                        "-lipidome_reorganization-",
                        "-lipidome_density-",
                        "-lipidome_composition-",
                    ]:
                        window[element].update(
                            disabled=not status.lipidome_prediction
                        )
                else:
                    for element in [
                        "-lipidome_predict-",
                        "-lipidome_report-",
                        "-lipidome_reset-",
                        "-lipidome_heatmap-",
                        "-lipidome_reorganization-",
                        "-lipidome_density-",
                        "-lipidome_composition-",
                    ]:
                        window[element].update(disabled=True)

        else:
            for element in [
                "-lipidome_predict-",
                "-lipidome_report-",
                "-lipidome_reset-",
                "-lipidome_heatmap-",
                "-lipidome_reorganization-",
                "-lipidome_density-",
                "-lipidome_composition-",
                "-global_run-",
                "-global_heatmap-",
                "-global_distance-",
                "-global_report-",
                "-global_reset-",
                "-class_run-",
                "-class_heatmap-",
                "-class_reorganization-",
                "-class_report-",
                "-class_reset-",
            ]:
                window[element].update(disabled=True)

    else:
        for element in [
            "-statistic_predict-",
            "-statistic_export-",
            "-statistic_report-",
            "-statistic_reset-",
            "-statistic_heatmap-",
            "-statistic_distribution-",
            "-lipidome_predict-",
            "-lipidome_report-",
            "-lipidome_reset-",
            "-lipidome_heatmap-",
            "-lipidome_reorganization-",
            "-lipidome_density-",
            "-lipidome_composition-",
            "-global_run-",
            "-global_heatmap-",
            "-global_distance-",
            "-global_report-",
            "-global_reset-",
            "-class_run-",
            "-class_heatmap-",
            "-class_reorganization-",
            "-class_report-",
            "-class_reset-",
        ]:
            window[element].update(disabled=True)
