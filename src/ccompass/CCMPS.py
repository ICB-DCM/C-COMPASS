"""C-COMPASS main window."""

import copy
import os
import pickle
import random
from datetime import datetime
from pathlib import Path
from tkinter import messagebox, simpledialog
from typing import Any

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict

from . import MOA, PPMS, RP

# tensorflow logging
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"


def default_status():
    status = {
        "fractionation_data": False,
        "tp_data": False,
        "lipidome_data": False,
        "lipidome_total": False,
        "marker_file": False,
        "marker_matched": False,
        "training": False,
        "proteome_prediction": False,
        "lipidome_prediction": False,
        "comparison_global": False,
        "comparison_class": False,
    }
    return status


class SessionModel(BaseModel):
    """Data for a C-COMPASS session."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Filepaths for fractionation data
    fract_paths: list[str] = []
    #: Fractionation column assignments
    #  filepath => [column ID, condition, replicate, fraction]
    fract_tables: dict[
        str, list[tuple[str, int | str, int | str, int | str]]
    ] = {}
    #: ??
    fract_pos: dict[str, list[int]] = {}
    #: Fractionation input files: filepath => DataFrame
    fract_indata: dict[str, pd.DataFrame] = {}
    #: Fractionation data for classification and visualization
    #  One DataFrame for each condition x replicate
    #  ("{condition}_Rep.{replicate}")
    fract_data: dict[str, dict[str, pd.DataFrame]] = {"class": {}, "vis": {}}
    #: ??
    #  for visualization and classification, each containing one DataFrame
    #  per condition with columns "{condition}_std_Fr.{fraction}"
    fract_std: dict[str, dict[str, pd.DataFrame]] = {"class": {}, "vis": {}}
    #: ??
    #  *something* => "{condition}" => "Rep.{replicate}" => DataFrame
    fract_intermediate: dict[str, dict[str, dict[str, pd.DataFrame]]] = {}
    #: Identifier column for the fractionation: filepath => column id
    fract_identifiers: dict[str, str] = {}
    #: Addition ("keep") columns from the fractionation data
    #  column ID => DataFrame
    fract_info: dict[str, pd.DataFrame] = {}
    #: Fractionation preprocessing parameters.
    #  global/classification/visualization
    #  "global"|"class"|"vis" => option => value
    fract_preparams: dict[str, dict[str, Any]] = PPMS.fract_default()
    #: Conditions in the fractionation data, including "[KEEP]"
    fract_conditions: list[str] = []
    #: Fractionation data for the different conditions x replicates
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_full: dict[str, pd.DataFrame] = {}
    #: Fractionation data after upsampling
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_full_up: dict[str, pd.DataFrame] = {}
    #: Marker abundance in the different fractions
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_marker: dict[str, pd.DataFrame] = {}
    #: Marker abundance in the different fractions for visualization
    #  "{condition}_median" => DataFrame
    fract_marker_vis: dict[str, pd.DataFrame] = {}
    #: Marker abundance in the different fractions after upsampling
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_marker_up: dict[str, pd.DataFrame] = {}
    #: ??
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_mixed_up: dict[str, pd.DataFrame] = {}
    #: ??
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_test: dict[str, pd.DataFrame] = {}

    #: Filepaths for total proteome data
    tp_paths: list[str] = []
    #: Total proteome column assignments
    #  filepath => [column ID, condition]
    tp_tables: dict[str, list[tuple[str, str]]] = {}
    #: ??
    tp_pos: dict[str, list[int]] = {}
    #: Total proteome input files: filepath => DataFrame
    tp_indata: dict[str, pd.DataFrame] = {}
    #: Total proteome data for the different conditions
    #  One DataFrame for each condition containing all replicates
    #  (column names are "{condition}_Rep.{replicate}")
    tp_data: dict[str, pd.DataFrame] = {}
    #: ??
    #  *something* => "{condition}" => DataFrame
    tp_intermediate: dict[str, dict[str, pd.DataFrame]] = {}
    #: Identifier column for the total proteome: filepath => column id
    tp_identifiers: dict[str, str] = {}
    #: ??
    tp_icorr: dict = {}
    #: ??
    tp_conditions: list = []
    #: ??
    tp_info: pd.DataFrame = pd.DataFrame()
    #: Total proteome preprocessing parameters
    tp_preparams: dict[str, Any] = PPMS.tp_default()

    #: Marker files, classes and annotations
    #  filepath => {'table'->pd.DataFrame,
    #  'identifier_col'-> column ID ("key column" in GUI),
    #  'class_col': column ID with class names in the marker file,
    #  'classes': list[str] class names
    #  }
    marker_sets: dict[str, dict[str, Any]] = {}
    #: Marker selection parameters
    marker_params: dict[str, Any] = {"how": "exclude", "what": "unite"}
    #: Mapping of compartment names to class names
    #  nan-values indicate that the compartment is not to be used
    marker_conv: dict[str, str | float] = {}
    #: Marker list "name" (gene name) => "class" (class name)
    marker_list: pd.DataFrame = pd.DataFrame()
    #: The column ID of the fractionation DataFrame that contains that is
    #  to be used for matching the markers
    marker_fractkey: str = "[IDENTIFIER]"

    #: SVM marker prediction
    # "{condition}_Rep.{replicate}" => DataFrame
    #  columns are the fractions + ["class", "svm_prediction", "svm_probability"]
    svm_marker: dict[str, pd.DataFrame] = {}
    #: SVM test data
    # "{condition}_Rep.{replicate}" => DataFrame
    #  columns are the fractions + ["class", "svm_prediction", "svm_probability"]
    svm_test: dict[str, pd.DataFrame] = {}
    #: SVM classification metrics for each condition x replicate
    # "{condition}_Rep.{replicate}" => dict(
    #   "accuracy" -> float,
    #   "precision" -> float,
    #   "recall" -> float,
    #   "f1" -> float,
    #   "confusion" -> pd.DataFrame,
    # )
    svm_metrics: dict[str, dict[str, Any]] = {}

    #: Neural network data
    # "{condition}_Rep.{replicate}" => dict(
    #  {w,W,x,X,y,Y,z,Z}_... => ...
    # )
    learning_xyz: dict[str, dict[str, Any]] = {}
    #: Nerural network hyperparameters
    NN_params: dict[str, Any] = PPMS.NN_default()

    #: SVM results (?)
    # "{condition}" => dict(
    #  "metrics" -> DataFrame,
    #  "SVM" -> dict("winner_combined" => DataFrame,
    #                "prob_combined" => DataFrame),
    # )
    results: dict[str, dict[str, Any]] = {}
    #: ??
    comparison: dict = {}
    #: Indicates which of the individual analysis steps
    #  have already been performed or not
    status: dict[str, bool] = default_status()

    def reset_global_changes(self):
        self.comparison = {}
        self.status["comparison_global"] = False
        self.status["comparison_class"] = False

    def reset_static_statistics(self):
        self.reset_global_changes()
        self.results = {}
        self.status["proteome_prediction"] = False
        self.status["lipidome_prediction"] = False

    def reset_input_tp(self):
        self.tp_paths = []
        self.tp_tables = {}
        self.tp_pos = {}
        self.tp_data = {}

    def reset_input_fract(self):
        self.fract_paths = []
        self.fract_tables = {}
        self.fract_pos = {}
        self.fract_data = {}

    def reset_infract(self):
        self.fract_indata = {}
        self.fract_identifiers = {}

    def reset_intp(self):
        self.tp_indata = {}
        self.tp_identifiers = {}

    def reset_fract(self):
        self.fract_data = {"class": [], "vis": []}
        self.fract_std = {"class": [], "vis": []}
        self.fract_intermediate = {}
        self.fract_info = {}
        self.fract_conditions = []

    def reset_tp(self):
        self.tp_data = {}
        self.tp_intermediate = {}
        self.tp_info = pd.DataFrame()
        self.tp_conditions = []
        self.tp_icorr = {}

    def reset_fractionation(self):
        self.reset_fract()
        self.reset_marker()
        self.status["fractionation_data"] = False

    def reset_classification(self):
        self.reset_static_statistics()

        self.reset_global_changes()
        self.svm_marker = {}
        self.svm_test = {}
        self.svm_metrics = {}
        self.learning_xyz = {}

        self.status["training"] = False

    def reset_marker(self):
        self.marker_list = pd.DataFrame()
        self.fract_marker = {}
        self.fract_marker_vis = {}
        self.fract_test = {}
        self.fract_full = {}
        self.reset_classification()
        self.status["marker_matched"] = False

    def to_numpy(self, filepath: Path | str):
        """Serialize using np.save."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, "wb") as f:
            np.save(f, self.model_dump(), allow_pickle=True)

    @classmethod
    def from_numpy(cls, filepath: Path | str):
        """Deserialize using np.load."""
        filepath = Path(filepath)
        with open(filepath, "rb") as f:
            data = np.load(f, allow_pickle=True).item()
            return cls(**data)


def create_fractionation_tab(fract_paths) -> sg.Tab:
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
                (fract_paths),
                size=(58, 1),
                key="-fractionation_path-",
                disabled=False,
                enable_events=True,
                readonly=True,
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
                button_color="dark grey",
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
                    [
                        sg.Text(
                            "...ready!",
                            font=("Arial", 8),
                            size=(40, 1),
                            key="-fractionation_status-",
                            visible=True,
                            justification="center",
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
                    # [sg.Button('Fract. Export...', size = (15,1), key = '-fractionation_export-', disabled = True, enable_events = True, button_color = 'grey')],
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


def create_total_proteome_tab(tp_paths) -> sg.Tab:
    """Create the "Total Proteome" tab."""

    layout_TP = [
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
                (tp_paths),
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
                button_color="dark grey",
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
                    [
                        sg.Text(
                            "...ready!",
                            font=("Arial", 8),
                            size=(40, 1),
                            key="-tp_status-",
                            visible=True,
                            justification="center",
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
        layout_TP,
        key="-total_tab-",
        expand_x=True,
        expand_y=True,
    )


def get_data_import_frame(fract_paths, tp_paths) -> sg.Frame:
    """Create the "Data Import" frame."""

    # The fractionation tab

    # layout_additional =     [[sg.Button('Add file...', size = (8,1), key = '-additional_add-', disabled = False, enable_events = True, button_color = 'grey'),
    #                         sg.Combo((fract_paths), size=(58, 1), key = '-additional_path-', disabled = False, enable_events = True, readonly = True),
    #                         sg.Button('Remove', size = (8,1), key = '-additional_remove-', disabled = False, enable_events = True, button_color = 'darkred')],
    #                         [sg.Table(values = [], num_rows = 15, headings = ['Sample', 'Condition', 'Replicate', 'Fraction'], col_widths = [34,10,10,9], max_col_width = 20, key = '-additional_table-', auto_size_columns = False, vertical_scroll_only = False)],
    #                         [sg.Button('Remove', size = (8,1), key = '-additional_edit_remove-', disabled = False, enable_events = True, button_color = 'dark red'),
    #                          sg.Button('Keep', size = (8,1), key = '-additional_edit_keep-', disabled = False, enable_events = True, button_color = 'dark grey', tooltip = ' Try to keep gene names! '),
    #                          sg.Button('Set Condition', size = (11,1), key = '-additional_edit_condition-', disabled = False, enable_events = True),
    #                          sg.Button('Set Replicate', size = (11,1), key = '-additional_edit_replicate-', disabled = False, enable_events = True),
    #                          sg.Button('Set Fractions', size = (11,1), key = '-additional_edit_fractions-', disabled = False, enable_events = True),
    #                          sg.Button('Set Identifier', size = (11,1), key = '-additional_edit_identifier-', disabled = False, enable_events = True, button_color = 'grey', tooltip = ' If possible, use protein groups! ')],
    #                         [sg.HSep()],
    #                         [sg.Column([
    #                              [sg.Button('Parameters...', size = (15,1), key = '-additional_parameters-', disabled = False, enable_events = True, button_color = 'black')],
    #                              [sg.Button('Reset Add.', size = (15,1), key = '-additional_reset-', disabled = True, enable_events = True, button_color = 'dark red')],
    #                              ], size = (140,70)),
    #                          sg.Column([
    #                              [sg.Button('Process Add.!', size = (30,1), key = '-additional_start-', disabled = False, enable_events = True, button_color = 'darkgreen')],
    #                              [sg.Text('...ready!', font = ('Arial', 8), size = (40,1), key = '-additional_status-', visible = True, justification = 'center')],
    #                              ], size = (260,70)),
    #                          sg.Column([
    #                              [sg.Button('Add. Summary', size = (15,1), key = '-additional_summary-', disabled = True, enable_events = True, button_color = 'grey')],
    #                              [sg.Button('Add. Export...', size = (15,1), key = '-additional_export-', disabled = True, enable_events = True, button_color = 'grey')],
    #                              ], size = (140,70))],
    #                         ]

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


def get_spatial_prediction_frame() -> sg.Frame:
    """Create the "Spatial Prediction" frame."""
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
                    "Train C-COMPASS!",
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
            [
                sg.Column(
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
                            sg.Text("Lipidome Fractionation:"),
                            sg.Push(),
                            sg.Text(
                                "none",
                                key="-status_fract_lipid-",
                                text_color="black",
                            ),
                        ],
                        [
                            sg.Text("Total Lipidome:"),
                            sg.Push(),
                            sg.Text(
                                "none",
                                key="-status_total_lipid-",
                                text_color="black",
                            ),
                        ],
                    ],
                    size=(240, 180),
                ),
                sg.Frame(
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
                ),
            ],
            # [sg.Frame(layout = [
            #     [sg.Button('calculate full static Statistics', size = (25,1), key = '-classification_statistics-', disabled = True, enable_events = True, button_color = 'dark blue'),
            #      sg.Text('(TP required!)', font = ('Helvetica', 8)),
            #      sg.Button('rouch statistics', size =(14,1), key = '-classification_statistics_rough-', disabled = False, enable_events = True, button_color = 'grey'),
            #      sg.Text('(w/o TP but biased!)', font = ('Helvetica', 8))],
            #     [sg.Text('missing', key = '-status_statistics-')],
            #     [sg.HSep()],
            #     [sg.Button('calculate conditional Changes', size = (25,1), key = '-classification_comparison-', disabled = False, enable_events = True, button_color = 'dark blue'),
            #      sg.Text('(TP required!)', font = ('Helvetica', 8)),
            #      sg.Button('rough Comparison', size = (14,1), key = '-classification_comparison_rough-', disabled = False, enable_events = True, button_color = 'grey'),
            #      sg.Text('(w/o TP but biased!)', font = ('Helvetica', 8))],
            #     [sg.Text('missing', key = '-status_comparison-')],
            #     ], title = 'Statistics', size = (590,150))],
            [sg.VPush()],
            [
                sg.Frame(
                    layout=[
                        [
                            sg.Frame(
                                layout=[
                                    [
                                        sg.Button(
                                            "Predict Proteome!",
                                            size=(32, 1),
                                            key="-statistic_predict-",
                                            disabled=True,
                                            enable_events=True,
                                            button_color="dark blue",
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
            ],
        ],
        title="Spatial Prediction",
        size=(600, 480),
        expand_x=True,
        expand_y=True,
    )


def get_marker_selection_frame() -> sg.Frame:
    """Create the "Marker Selection" frame."""
    return sg.Frame(
        layout=[
            [
                sg.Frame(
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
                                size=(72, 70),
                            ),
                            sg.Column(
                                layout=[
                                    [
                                        sg.Text("key column:\t"),
                                        sg.Combo(
                                            [],
                                            size=(10, 1),
                                            key="-marker_key-",
                                            enable_events=True,
                                            readonly=True,
                                        ),
                                    ],
                                    [
                                        sg.Text("class column:\t"),
                                        sg.Combo(
                                            [],
                                            size=(10, 1),
                                            key="-marker_class-",
                                            enable_events=True,
                                            readonly=True,
                                        ),
                                    ],
                                ],
                                size=(220, 70),
                            ),
                        ],
                    ],
                    title="Import",
                    size=(320, 190),
                ),
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
                            sg.Text("Fract. Key:\t"),
                            sg.Combo(
                                ["[IDENTIFIER]"],
                                key="-marker_fractkey-",
                                size=(18, 1),
                                readonly=True,
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
                    size=(280, 180),
                ),
            ],
        ],
        title="Marker Selection",
        size=(620, 220),
        expand_x=True,
        expand_y=True,
    )


def get_conditional_comparison_frame() -> sg.Frame:
    """Create the "Conditional Comparison" frame."""
    return sg.Frame(
        layout=[
            [
                sg.Frame(
                    layout=[
                        # [sg.Text('Statistics required!')],
                        [
                            sg.Button(
                                "Calculate global Changes!",
                                size=(32, 1),
                                key="-global_run-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark blue",
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
                ),
                sg.Frame(
                    layout=[
                        [
                            sg.Button(
                                "Calculate class-centric Changes!",
                                size=(32, 1),
                                key="-class_run-",
                                disabled=True,
                                enable_events=True,
                                button_color="dark blue",
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
                ),
            ]
        ],
        title="Conditional Comparison",
        size=(600, 220),
        expand_x=True,
        expand_y=True,
    )


def create_main_window(model: SessionModel) -> sg.Window:
    """Create the C-COMPASS main window."""

    # The main menu
    menu_def = [
        ["Session", ["Save...", "Open...", "New", "Exit"]],
        ["Help", ["About...", "Open Website", "Manual"]],
    ]

    layout_CCMPS = [
        [sg.Menu(menu_def, tearoff=False)],
        [
            get_data_import_frame(
                fract_paths=model.fract_paths, tp_paths=model.tp_paths
            ),
            get_spatial_prediction_frame(),
        ],
        [
            get_marker_selection_frame(),
            get_conditional_comparison_frame(),
        ],
    ]

    main_window = sg.Window(
        "C-COMPASS", layout_CCMPS, size=(1260, 720), resizable=True
    )
    return main_window


def main():
    """The entry point for the C-COMPASS application."""
    sg.theme("Dark Blue 3")

    model = SessionModel()
    controller = MainController(model=model)
    controller.run()

    # import dill
    # filepath = 'session.pkl'
    # dill.dump_session(filepath) # Save the session
    # dill.load_session(filepath) # Load the session


class MainController:
    """The main controller for the C-COMPASS application."""

    def __init__(self, model: SessionModel):
        self.model = model
        self.main_window = create_main_window(model=model)

    def run(self):
        """Run the C-COMPASS application."""

        # The event loop
        while True:
            event, values_CCMPS = self.main_window.read()
            # refresh_window(window_CCMPS, status)

            # if status['fractionation_data']:
            #     window_CCMPS['-status_fract-'].Update('ready')
            # else:
            #     window_CCMPS['-status_fract-'].Update('missing')

            if event == "-fractionation_add-":
                fract_add(values_CCMPS, self.main_window, model=self.model)
            elif event == "-fractionation_remove-":
                if values_CCMPS["-fractionation_path-"]:
                    fract_rem(
                        values_CCMPS,
                        self.main_window,
                        model=self.model,
                    )
            elif event == "-fractionation_path-":
                fract_refreshtable(
                    self.main_window,
                    self.model.fract_tables[
                        values_CCMPS["-fractionation_path-"]
                    ],
                )
            elif event == "-fractionation_edit_remove-":
                if values_CCMPS["-fractionation_table-"]:
                    fract_defrem(
                        values_CCMPS, self.main_window, self.model.fract_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_keep-":
                if values_CCMPS["-fractionation_table-"]:
                    fract_defkeep(
                        values_CCMPS, self.main_window, self.model.fract_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_condition-":
                if values_CCMPS["-fractionation_table-"]:
                    fract_defcon(
                        values_CCMPS, self.main_window, self.model.fract_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_replicate-":
                if values_CCMPS["-fractionation_table-"]:
                    fract_defrep(
                        values_CCMPS, self.main_window, self.model.fract_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_fractions-":
                if values_CCMPS["-fractionation_table-"]:
                    fract_deffract(
                        values_CCMPS, self.main_window, self.model.fract_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_edit_identifier-":
                if values_CCMPS["-fractionation_table-"]:
                    self.model.fract_identifiers = fract_defident(
                        values_CCMPS,
                        self.main_window,
                        self.model.fract_tables,
                        self.model.fract_pos,
                        self.model.fract_identifiers,
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-fractionation_parameters-":
                self.model.fract_preparams = PPMS.PPMS_exec_fract(
                    self.model.fract_preparams
                )
            elif event == "-fractionation_reset-":
                sure = sg.popup_yes_no(
                    "Reset Fractionation Pre-Processing? You have to run it again to use your data."
                )
                if sure == "Yes":
                    self.model.reset_fractionation()
                    fract_buttons(self.main_window, False)
                    # window_CCMPS['-marker_fractkey-'].Update(values = ['[IDENTIFIER]'] + list(fract_info))
                    # enable_markersettings(window_CCMPS, True)

                    self.main_window["-marker_fractkey-"].Update(
                        values=["[IDENTIFIER]"], value=""
                    )
                    # window_CCMPS['-classification_SVM-'].Update(disabled = True)

                else:
                    pass
            elif event == "-fractionation_start-":
                if self.model.fract_paths:
                    from .FDP import FDP_exec

                    (
                        self.model.fract_data,
                        self.model.fract_std,
                        self.model.fract_intermediate,
                        self.model.fract_info,
                        self.model.fract_conditions,
                    ) = FDP_exec(
                        self.main_window,
                        self.model.fract_tables,
                        self.model.fract_preparams,
                        self.model.fract_identifiers,
                        self.model.fract_data,
                        self.model.fract_std,
                        self.model.fract_intermediate,
                        self.model.fract_info,
                        self.model.fract_conditions,
                        self.model.fract_indata,
                    )
                    self.main_window["-marker_fractkey-"].Update(
                        values=["[IDENTIFIER]"] + list(self.model.fract_info)
                    )
                    if self.model.fract_data["class"]:
                        self.model.status["fractionation_data"] = True
                    #     fract_buttons(window_CCMPS, True)
                else:
                    messagebox.showerror(
                        "No dataset!", "Please import a fractionation dataset."
                    )
            elif event == "-fractionation_summary-":
                RP.RP_gradient_heatmap(self.model.fract_data)
                # FSD.FSD_exec(fract_preparams, fract_data)
            # if event_CCMPS == '-fractionation_export-':
            #     fract_export(values_CCMPS, fract_data, fract_info)

            elif event == "-tp_add-":
                tp_add(
                    values_CCMPS,
                    self.main_window,
                    self.model.tp_paths,
                    self.model.tp_tables,
                    self.model.tp_indata,
                    self.model.tp_pos,
                    self.model.tp_identifiers,
                )
            elif event == "-tp_remove-":
                if values_CCMPS["-tp_path-"]:
                    tp_rem(
                        values_CCMPS,
                        self.main_window,
                        self.model.tp_paths,
                        self.model.tp_tables,
                    )
            elif event == "-tp_path-":
                tp_refreshtable(
                    self.main_window,
                    self.model.tp_tables[values_CCMPS["-tp_path-"]],
                )
            elif event == "-tp_edit_remove-":
                if values_CCMPS["-tp_table-"]:
                    tp_defrem(
                        values_CCMPS, self.main_window, self.model.tp_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-tp_edit_keep-":
                if values_CCMPS["-tp_table-"]:
                    tp_defkeep(
                        values_CCMPS, self.main_window, self.model.tp_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-tp_edit_condition-":
                if values_CCMPS["-tp_table-"]:
                    tp_defcon(
                        values_CCMPS, self.main_window, self.model.tp_tables
                    )
                else:
                    messagebox.showerror("Error", "Select (a) row(s).")
            elif event == "-tp_edit_identifier-":
                if values_CCMPS["-tp_table-"]:
                    self.model.tp_identifiers = tp_defident(
                        values_CCMPS,
                        self.main_window,
                        self.model.tp_tables,
                        self.model.tp_pos,
                        self.model.tp_identifiers,
                    )
            elif event == "-tp_parameters-":
                self.model.tp_preparams = PPMS.PPMS_exec_TP(
                    self.model.tp_preparams
                )
            elif event == "-tp_reset-":
                sure = sg.popup_yes_no(
                    "Reset TotalProteome Pre-Processing? You have to run it again to use your data."
                )
                if sure == "Yes":
                    self.model.reset_tp()

                    self.model.status["tp_data"] = False
                    if self.model.status["comparison_class"]:
                        self.model.results, self.model.comparison = (
                            MOA.class_reset(
                                self.model.results, self.model.comparison
                            )
                        )
                        self.model.status["comparison_class"] = False
                else:
                    pass
            elif event == "-tp_start-":
                if self.self.model.tp_paths:
                    from .TPP import TPP_exec

                    (
                        self.model.tp_data,
                        self.model.tp_intermediate,
                        self.model.tp_info,
                        self.model.tp_conditions,
                        self.model.tp_icorr,
                    ) = TPP_exec(
                        self.main_window,
                        self.model.tp_data,
                        self.model.tp_tables,
                        self.model.tp_preparams,
                        self.model.tp_identifiers,
                        self.model.tp_intermediate,
                        self.model.tp_info,
                        self.model.tp_icorr,
                        self.model.tp_indata,
                        self.model.tp_conditions,
                    )

                    if self.model.tp_data:
                        self.model.status["tp_data"] = True
                        # tp_buttons(window_CCMPS, True)
                else:
                    messagebox.showerror(
                        "No dataset!", "Please import a TP dataset."
                    )
            elif event == "-tp_export-":
                tp_export(self.model.tp_data, self.model.tp_info)

            elif event == "-marker_add-":
                marker_add(
                    self.main_window, values_CCMPS, self.model.marker_sets
                )
                event, values_CCMPS = self.main_window.read(timeout=50)
                if self.model.marker_sets:
                    self.model.status["marker_file"] = True
                else:
                    self.model.status["marker_file"] = False
            elif event == "-marker_remove-":
                try:
                    marker_remove(
                        self.main_window, values_CCMPS, self.model.marker_sets
                    )
                except Exception:
                    pass
                if self.model.marker_sets:
                    self.model.status["marker_file"] = True
                else:
                    self.model.status["marker_file"] = False
            elif event == "-marker_list-":
                refresh_markercols(
                    self.main_window, values_CCMPS, self.model.marker_sets
                )
            elif event == "-marker_key-":
                marker_setkey(values_CCMPS, self.model.marker_sets)
            elif event == "-marker_class-":
                self.model.marker_conv = marker_setclass(
                    values_CCMPS, self.model.marker_sets
                )

            elif event == "-marker_parameters-":
                self.model.marker_params = PPMS.PPMS_exec_marker(
                    self.model.marker_params
                )

            elif event == "-marker_manage-":
                from .CM import CM_exec

                if check_markers(self.model.marker_sets):
                    self.model.marker_conv = CM_exec(
                        self.model.marker_sets, self.model.marker_conv
                    )
                else:
                    messagebox.showerror(
                        "Error", "Please define key and class column."
                    )

            elif event == "-marker_test-":
                if check_markers(self.model.marker_sets):
                    from .TM import TM_exec

                    try:
                        self.model.marker_list = create_markerlist(
                            self.model.marker_sets,
                            self.model.marker_conv,
                            self.model.marker_params,
                        )
                        TM_exec(
                            self.model.fract_data,
                            self.model.fract_info,
                            self.model.marker_list,
                            values_CCMPS["-marker_fractkey-"],
                        )
                    except Exception:
                        messagebox.showerror(
                            "Error",
                            "Something is wrong with your marker list.",
                        )
                else:
                    messagebox.showerror(
                        "Error", "Please define key and class column."
                    )
            elif event == "-marker_profiles-":
                if check_markers(self.model.marker_sets):
                    from .SM import SM_exec

                    try:
                        self.model.marker_list = create_markerlist(
                            self.model.marker_sets,
                            self.model.marker_conv,
                            self.model.marker_params,
                        )
                        SM_exec(
                            self.model.fract_data,
                            self.model.fract_info,
                            self.model.marker_list,
                            values_CCMPS["-marker_fractkey-"],
                        )
                    except Exception:
                        messagebox.showerror(
                            "Error",
                            "Something is wrong with your marker list.",
                        )
                else:
                    messagebox.showerror(
                        "Error", "Please define key and class column."
                    )

            elif event == "-marker_accept-":
                if values_CCMPS["-marker_list-"] == []:
                    messagebox.showerror(
                        "Error", "Please import at least one Marker List!"
                    )
                elif self.model.fract_data["class"] == []:
                    messagebox.showerror(
                        "Error", "Please import Fractionation Data first!"
                    )
                elif (
                    values_CCMPS["-marker_fractkey-"] == ""
                    or values_CCMPS["-marker_class-"] == ""
                    or values_CCMPS["-marker_key-"] == ""
                ):
                    messagebox.showerror(
                        "Error", "Please select key and class columns!"
                    )
                else:
                    from .MOP import create_fullprofiles

                    try:
                        # print('Starting try block')
                        self.model.marker_list = create_markerlist(
                            self.model.marker_sets,
                            self.model.marker_conv,
                            self.model.marker_params,
                        )
                        # print('check1: marker_list created')
                        (
                            self.model.fract_marker,
                            self.model.fract_marker_vis,
                            self.model.fract_test,
                            self.model.classnames,
                        ) = create_markerprofiles(
                            self.model.fract_data,
                            values_CCMPS["-marker_fractkey-"],
                            self.model.fract_info,
                            self.model.marker_list,
                        )
                        # print('check2: marker profiles created')
                        self.model.fract_full = create_fullprofiles(
                            self.model.fract_marker, self.model.fract_test
                        )
                        self.model.status["marker_matched"] = True
                        # print('check3: full profiles created')
                        # enable_markersettings(window_CCMPS, False)
                        # print('check4: marker settings enabled')
                        # window_CCMPS['-classification_MOP-'].Update(disabled = False)
                        # print('check5: classification MOP updated')
                        # window_CCMPS['-classification_SVM-'].Update(disabled = False)
                        # print('check6: classification SVM updated')
                    except Exception as e:
                        print(f"An error occurred: {e}")
                        import traceback

                        traceback.print_exc()
                        messagebox.showerror(
                            "Error", "Incompatible Fractionation Key!"
                        )

            elif event == "-marker_reset-":
                self.model.reset_marker()
                # enable_markersettings(window_CCMPS, True)
                # window_CCMPS['-classification_MOP-'].Update(disabled = True)
                # window_CCMPS['-classification_SVM-'].Update(disabled = True)

            elif event == "-classification_parameters-":
                self.model.NN_params = PPMS.PPMS_exec_NN(self.model.NN_params)

            elif event == "-classification_MOP-":
                from .MOP import MOP_exec

                (
                    self.model.learning_xyz,
                    self.model.NN_params,
                    self.model.fract_full_up,
                    self.model.fract_marker_up,
                    self.model.fract_mixed_up,
                    self.model.fract_unmixed_up,
                    self.model.svm_marker,
                    self.model.svm_test,
                    self.model.svm_metrics,
                ) = MOP_exec(
                    self.model.fract_conditions,
                    self.model.fract_full,
                    self.model.fract_marker,
                    self.model.fract_test,
                    self.model.fract_std,
                    self.model.fract_info,
                    values_CCMPS["-marker_fractkey-"],
                    self.model.NN_params,
                )
                # window_CCMPS['-classification_statistics-'].Update(disabled = False)
                # window_CCMPS['-status_comparison-'].Update('done!')
                self.model.status["training"] = True

            elif event == "-classification_reset-":
                self.model.reset_classification()

            elif event == "-statistic_predict-":
                self.model.results = MOA.stats_proteome(
                    self.model.learning_xyz,
                    self.model.NN_params,
                    self.model.fract_data,
                    self.model.fract_conditions,
                )
                self.model.status["proteome_prediction"] = True

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
                filename = sg.popup_get_file(
                    "Import Statistics",
                    no_window=True,
                    file_types=(("Pickle", "*.pkl"),),
                )
                if filename:
                    with open(filename, "rb") as file:
                        results_new = pickle.load(file)
                    try:
                        for condition in results_new:
                            if condition in self.model.results:
                                messagebox.showerror(
                                    "Error",
                                    "There are already statistics for "
                                    + condition
                                    + " in your current session.",
                                )
                            else:
                                self.model.results[condition] = copy.deepcopy(
                                    results_new[condition]
                                )
                        self.model.status[
                            "proteome_prediction"
                        ] = self.model.status["training"] = True
                    except Exception:
                        messagebox.showerror(
                            "Error", "Incompatible file type!"
                        )

            elif event == "-statistic_report-":
                export_folder = sg.popup_get_folder("Statistics Report")
                if export_folder:
                    for condition in self.model.results:
                        fname = Path(
                            export_folder, f"CCMPS_statistics_{condition}.xlsx"
                        )
                        selected_columns = [
                            col
                            for col in self.model.results[condition][
                                "metrics"
                            ].columns
                            if col.startswith("fCC_")
                        ] + ["SVM_winner", "fNN_winner", "marker"]
                        df_out = self.model.results[condition]["metrics"][
                            selected_columns
                        ]
                        df_out.columns = [
                            col.replace("fCC_", "CC_ClassContribution_")
                            if col.startswith("fCC_")
                            else "C-CMPS_MainClass"
                            if col == "fNN_winner"
                            else col
                            for col in df_out.columns
                        ]
                        df_out.to_excel(fname, index=True)

            elif event == "-global_report-":
                export_folder = sg.popup_get_folder("Global Changes Report")
                if export_folder:
                    for comb in self.model.comparison:
                        fname = Path(
                            export_folder,
                            f"CCMPS_comparison_{comb[0]}_{comb[1]}.xlsx",
                        )
                        selected_columns = [
                            col
                            for col in self.model.comparison[comb][
                                "metrics"
                            ].columns
                            if col.startswith("fRL_")
                        ] + ["fRLS", "DS", "P(t)_RLS"]
                        df_out = self.model.comparison[comb]["metrics"][
                            selected_columns
                        ]
                        df_out.columns = [
                            col.replace("fRL_", "RL_Relocalization_")
                            if col.startswith("fRL_")
                            else "RLS_ReLocalizationScore"
                            if col == "fRLS"
                            else "DS_DistanceScore"
                            if col == "DS"
                            else "P-Value"
                            if col == "P(t)_RLS"
                            else col
                            for col in df_out.columns
                        ]
                        df_out.to_excel(fname, index=True)

            elif event == "-class_report-":
                export_folder = sg.popup_get_folder(
                    "Class-centric Changes Report"
                )
                if export_folder:
                    for condition in self.model.results:
                        fname = Path(
                            export_folder,
                            f"CCMPS_ClassComposition_{condition}.xlsx",
                        )
                        selected_columns = [
                            col
                            for col in self.model.results[condition][
                                "metrics"
                            ].columns
                            if col.startswith("nCPA")
                        ] + ["TPA"]
                        df_out = self.model.results[condition]["metrics"][
                            selected_columns
                        ]
                        df_out.columns = [
                            col.replace(
                                "nCPA_imp_",
                                "nCPA_normalizedClasscentrigProteinAmount_",
                            )
                            if col.startswith("nCPA_")
                            else "TPA_TotalProteinAmount"
                            if col == "TPA"
                            else col
                            for col in df_out.columns
                        ]
                        df_out.to_excel(fname, index=True)
                    for comb in self.model.comparison:
                        fname = Path(
                            export_folder,
                            f"CCMPS_ClassComparison_{comb[0]}_{comb[1]}.xlsx",
                        )
                        selected_columns = [
                            col
                            for col in self.model.comparison[comb][
                                "metrics"
                            ].columns
                            if col.startswith("nCFC_")
                        ]
                        df_out = self.model.comparison[comb]["metrics"][
                            selected_columns
                        ]
                        df_out.columns = [
                            col.replace(
                                "nCFC_",
                                "nCFC_normalizedClasscentricFoldChange_",
                            )
                            if col.startswith("nCFC_")
                            else col
                            for col in df_out.columns
                        ]
                        df_out.to_excel(fname, index=True)

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

            elif event == "-class_reorganization-":
                RP.RP_class_reorganization(self.model.comparison)

            elif event == "-global_run-":
                self.model.comparison = MOA.global_comparison(
                    self.model.results
                )
                self.model.status["comparison_global"] = True

            elif event == "-global_reset-":
                self.model.reset_global_changes()

            elif event == "-class_run-":
                self.model.comparison = MOA.class_comparison(
                    self.model.tp_data,
                    self.model.fract_conditions,
                    self.model.results,
                    self.model.comparison,
                )
                self.model.status["comparison_class"] = True

            elif event == "-class_reset-":
                self.model.results, self.model.comparison = MOA.class_reset(
                    self.model.results, self.model.comparison
                )
                self.model.status["comparison_class"] = False

            # if event_CCMPS == '-classification_comparison-':
            #     # results = MOP_stats.comp_exec(learning_xyz, results)
            #     comparison = MOP_stats.comp_exec3('deep', results, learning_xyz)

            # if event_CCMPS == '-classification_comparison_rough-':
            #     comparison = MOP_stats.comp_exec3('rough', results, learning_xyz)

            elif event == "-export_statistics-":
                export_folder = sg.popup_get_folder("Export Statistics")
                if export_folder:
                    for condition in self.model.results:
                        fname = Path(
                            export_folder, f"CCMPS_statistics_{condition}.tsv"
                        )
                        df_out = pd.merge(
                            self.model.fract_data["vis"][
                                condition + "_median"
                            ],
                            self.model.results[condition]["metrics"],
                            left_index=True,
                            right_index=True,
                            how="outer",
                        )
                        for colname in self.model.fract_info:
                            df_out = pd.merge(
                                df_out,
                                self.model.fract_info[colname],
                                left_index=True,
                                right_index=True,
                                how="left",
                            )
                        df_out.to_csv(
                            fname,
                            sep="\t",
                            index=True,
                            index_label="Identifier",
                        )
                        # results[condition]['metrics'].to_csv(fname, sep='\t', index=True, index_label='Identifier')

            elif event == "-export_comparison-":
                export_folder = sg.popup_get_folder("Export Statistics")
                if export_folder:
                    for comb in self.model.comparison:
                        fname = Path(
                            export_folder,
                            f"CCMPS_comparison_{comb[0]}_{comb[1]}.tsv",
                        )
                        df_out = pd.DataFrame(
                            index=self.model.comparison[comb][
                                "intersection_data"
                            ].index
                        )
                        df_out = pd.merge(
                            df_out,
                            self.model.comparison[comb]["metrics"],
                            left_index=True,
                            right_index=True,
                            how="left",
                        )
                        for colname in self.model.fract_info:
                            df_out = pd.merge(
                                df_out,
                                self.model.fract_info[colname],
                                left_index=True,
                                right_index=True,
                                how="left",
                            )
                        df_out.to_csv(
                            fname,
                            sep="\t",
                            index=True,
                            index_label="Identifier",
                        )
                        # comparison[comb]['metrics'].to_csv(fname, sep='\t', index=True, index_label='Identifier')

            # if event_CCMPS == '-export_statistics-':
            # path = sg.popup_get_folder('Select a folder')
            # print(path)
            # if fract_data['class']:
            #     path = sg.FolderBrowse()
            #     for condition in fract_data['class']:
            #         df_export_stats = pd.merge(fract_data['class'][condition], fract_data['vis'][condition], left_index = True, right_index = True, how = 'outer')
            #         df_export_stats = pd.merge(df_export_stats, results[condition]['metrics'], left_index = True, right_index = True, how = 'outer')
            #         full_path = path+'/CCMPS_export_CombinedData.tsv'
            #         df_export_stats.to_csv(full_path, sep = '\t', index = False)
            #     for comp in results:
            #         df_export_part = results[comp]
            #         full_path = path+'/CCMPS_export_'+comp+'.tsv'
            #         df_export_part.to_csv(full_path, sep = '\t', index = False)

            elif event == "Save...":
                session_save(self.model)
            elif event == "Open...":
                filename = sg.popup_get_file(
                    "Open Session",
                    no_window=True,
                    file_types=(("Numpy", "*.npy"),),
                )
                if filename:
                    session_open(
                        self.main_window,
                        values_CCMPS,
                        filename,
                        model=self.model,
                    )
                    # window_CCMPS['-marker_tpkey-'].Update(values = ['[IDENTIFIER]'] + tp_info.columns.tolist())
                    self.main_window["-marker_fractkey-"].Update(
                        values=["[IDENTIFIER]"] + list(self.model.fract_info),
                        value=self.model.marker_fractkey,
                    )
            elif event == "New":
                sure = sg.popup_yes_no(
                    "Are you sure to close the session and start a new one?"
                )
                if sure == "Yes":
                    self.model.reset_input_fract()
                    fract_clearinput(self.main_window)
                    tp_clearinput(self.main_window)
                    self.model.reset_infract()
                    self.model.reset_fract()
                    self.model.reset_input_tp()
                    self.model.reset_intp()
                    self.model.reset_tp()

                    self.model.fract_preparams = PPMS.fract_default()
                    self.model.tp_preparams = PPMS.tp_default()

                    self.model.marker_sets = {}
                    self.model.marker_conv = {}
                    self.model.fract_full_up = {}
                    self.model.fract_marker_up = {}
                    self.model.fract_mixed_up = {}
                    self.model.reset_marker()
                    self.main_window["-marker_list-"].Update(values=[])
                    self.main_window["-marker_key-"].Update(values=[])
                    self.main_window["-marker_class-"].Update(values=[])
                    self.model.marker_params = {
                        "how": "exclude",
                        "what": "unite",
                    }
                    self.model.NN_params = PPMS.NN_default()

                    self.model.status = default_status()

                    self.main_window["-marker_fractkey-"].Update(
                        values=["[IDENTIFIER]"] + list(self.model.fract_info)
                    )
                else:
                    pass

            # -----------------------------------------------------------------------------------------------------------------------------

            elif event == sg.WIN_CLOSED or event == "Exit":
                break

            refresh_window(self.main_window, self.model.status)

        self.main_window.close()


def fract_refreshtable(window, table):
    window["-fractionation_table-"].Update(values=table)
    return


def tp_refreshtable(window, table):
    window["-tp_table-"].Update(values=table)
    return


def fract_modifytable(title, prompt, values, fract_tables, pos, q, ask):
    if values["-fractionation_table-"]:
        path = values["-fractionation_path-"]
        table = fract_tables[path]
        if ask == "integer":
            value = simpledialog.askinteger(title, prompt)
            p = 0
            if value:
                for i in values["-fractionation_table-"]:
                    table[i][pos] = value + p
                    p = p + q
                fract_tables[path] = table
        elif ask == "string":
            value = simpledialog.askstring(title, prompt)
            if value:
                for i in values["-fractionation_table-"]:
                    table[i][pos] = value
                fract_tables[path] = table
    else:
        messagebox.showerror("Error", "Select (a) sample(s).")
    return (values, fract_tables)


def fract_buttons(window, status):
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
        window[button].Update(disabled=status)
    for button in inactive:
        window[button].Update(disabled=not status)
    if status:
        window["-fractionation_status-"].Update(
            value="done!", text_color="dark green"
        )
    else:
        window["-fractionation_status-"].Update(
            value="...ready!", text_color="white"
        )
    return


def tp_buttons(window, status):
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
        window[button].Update(disabled=status)
    for button in inactive:
        window[button].Update(disabled=not status)
    if status:
        window["-tp_status-"].Update(value="done!", text_color="dark green")
    else:
        window["-tp_status-"].Update(value="...ready!", text_color="white")


def fract_clearinput(window):
    window["-fractionation_path-"].Update(values=[])
    window["-fractionation_table-"].Update(values=[])
    return


def tp_clearinput(window):
    window["-tp_path-"].Update(values=[])
    window["-tp_table-"].Update(values=[])
    return


def fract_add(
    values,
    window,
    model: SessionModel,
):
    filename = sg.popup_get_file(
        "Chose dataset",
        no_window=True,
        file_types=(
            ("Tab Separated Values", "*.tsv"),
            ("Text (tab delimited)", "*.txt"),
        ),
    )
    if filename:
        model.fract_paths.append(filename)
        window["-fractionation_path-"].Update(
            values=model.fract_paths, value=filename
        )
        data = pd.read_csv(filename, sep="\t", header=0)
        data = data.replace("NaN", np.nan)
        data = data.replace("Filtered", np.nan)
        colnames = data.columns.values.tolist()
        table = []
        for name in colnames:
            namelist = [name, "", "", ""]
            table.append(namelist)
        model.fract_tables[filename] = table
        model.fract_indata[filename] = data
        model.fract_pos[filename] = []
        model.fract_identifiers[filename] = []

        fract_refreshtable(window, table)
    return


def fract_rem(values, window, model: SessionModel):
    sure = sg.popup_yes_no("Remove data from list?")
    if sure == "Yes":
        model.fract_paths.remove(values["-fractionation_path-"])
        del model.fract_tables[values["-fractionation_path-"]]
        if model.fract_paths:
            curr = model.fract_paths[0]
            fract_refreshtable(window, model.fract_tables[curr])
        else:
            curr = []
            fract_refreshtable(window, curr)
        window["-fractionation_path-"].Update(
            values=model.fract_paths, value=curr
        )
    else:
        pass
    return


def fract_defrem(values, window, fract_tables):
    path = values["-fractionation_path-"]
    selected = values["-fractionation_table-"]
    table = fract_tables[path]
    for index in sorted(selected, reverse=True):
        del table[index]
    fract_tables[path] = table
    window["-fractionation_table-"].Update(values=fract_tables[path])
    return


def fract_defkeep(values, window, fract_tables):
    path = values["-fractionation_path-"]
    table = fract_tables[path]
    for pos in values["-fractionation_table-"]:
        table[pos][1] = "[KEEP]"
        table[pos][2] = "-"
        table[pos][3] = "-"
    fract_tables[path] = table
    window["-fractionation_table-"].Update(values=fract_tables[path])
    return


def fract_defcon(values, window, fract_tables):
    values, fract_tables = fract_modifytable(
        "Set Condition",
        "Condition Name:",
        values,
        fract_tables,
        1,
        0,
        "string",
    )
    window["-fractionation_table-"].Update(
        values=fract_tables[values["-fractionation_path-"]]
    )
    return


def fract_defrep(values, window, fract_tables):
    values, fract_tables = fract_modifytable(
        "Set Replicate",
        "Replicate Number:",
        values,
        fract_tables,
        2,
        0,
        "integer",
    )
    window["-fractionation_table-"].Update(
        values=fract_tables[values["-fractionation_path-"]]
    )
    return


def fract_deffract(values, window, fract_tables):
    values, fract_tables = fract_modifytable(
        "Set Fractions",
        "FIRST Fraction Number:",
        values,
        fract_tables,
        3,
        1,
        "integer",
    )
    window["-fractionation_table-"].Update(
        values=fract_tables[values["-fractionation_path-"]]
    )
    return


def fract_defident(values, window, input_tables, ident_pos, identifiers):
    pos = values["-fractionation_table-"]
    if pos:
        if len(pos) > 1:
            messagebox.showerror("Error", "Please set only one Identifier!")
        elif len(pos) == 1:
            path = values["-fractionation_path-"]
            table = input_tables[path]
            if ident_pos[path]:
                table[ident_pos[path][0]][1] = ""
                table[ident_pos[path][0]][2] = ""
                table[ident_pos[path][0]][3] = ""
            identifiers[path] = table[pos[0]][0]
            ident_pos[path] = pos
            table[pos[0]][1] = "[IDENTIFIER]"
            table[pos[0]][2] = "-"
            table[pos[0]][3] = "-"
            input_tables[path] = table
            window["-fractionation_table-"].Update(
                values=input_tables[values["-fractionation_path-"]]
            )
        else:
            messagebox.showerror("Error", "No sample selected.")
    return identifiers


def fract_export(values, data, protein_info):
    export_folder = sg.popup_get_folder("Export Folder")
    if export_folder:
        experiment = simpledialog.askstring("Export", "Experiment Name: ")
        now = datetime.now()
        time = now.strftime("%Y%m%d%H%M%S")
        # export_full = {'class' : pd.DataFrame(),
        #                'vis' : pd.DataFrame()}
        for way in data:
            data_way = pd.DataFrame()
            for condition in data[way]:
                path = (
                    export_folder
                    + "/"
                    + time
                    + "_"
                    + experiment
                    + "_"
                    + way
                    + "_"
                    + condition
                    + ".txt"
                )
                export_data = data[way][condition]
                for info in protein_info:
                    export_data = pd.merge(
                        export_data,
                        protein_info[info],
                        left_index=True,
                        right_index=True,
                        how="left",
                    )
                export_data.to_csv(
                    path,
                    header=True,
                    index=True,
                    index_label="Identifier",
                    sep="\t",
                    mode="a",
                )

                data_way = pd.merge(
                    data_way,
                    data[way][condition],
                    left_index=True,
                    right_index=True,
                    how="outer",
                )

                # for replicate in data[way][condition]:
                #     path = export_folder + '/' + time + '_' + experiment + '_' + way + '_' + condition + '_' + replicate + '.txt'
                #     data[way][condition][replicate].to_csv(path, header = True, index = True, index_label = 'Identifier', sep = '\t', mode = 'a')
                #     export_full[way] = pd.merge(export_full[way], data[way][condition][replicate], left_index = True, right_index = True, how = 'outer')

            # for info in protein_info:
            #     export_full[way] = pd.merge(export_full[way], protein_info[info], left_index = True, right_index = True, how = 'left')

            for info in protein_info:
                data_way = pd.merge(
                    data_way,
                    protein_info[info],
                    left_index=True,
                    right_index=True,
                    how="left",
                )

            path = (
                export_folder
                + "/"
                + time
                + "_"
                + experiment
                + "_"
                + way
                + "_"
                + "COMBINED"
                + ".txt"
            )
            data_way.to_csv(
                path,
                header=True,
                index=True,
                index_label="Identifier",
                sep="\t",
                mode="a",
            )
    return


# ------------------------------------------------------------------------------


def is_float(element):
    try:
        float(element)
        return True
    except (ValueError, TypeError):
        return False


def convert_to_float(x):
    try:
        return float(x)
    except ValueError:
        return x


def tp_add(
    values, window, tp_paths, tp_tables, tp_indata, tp_pos, tp_identifiers
):
    filename = sg.popup_get_file(
        "Chose dataset",
        no_window=True,
        file_types=(
            ("Tab Separated Values", "*.tsv"),
            ("Text (tab delimited)", "*.txt"),
        ),
    )
    if filename:
        tp_paths.append(filename)
        window["-tp_path-"].Update(values=tp_paths, value=filename)
        data = pd.read_csv(filename, sep="\t", header=0)
        data = data.replace("NaN", np.nan)
        data = data.replace("Filtered", np.nan)
        data = data.applymap(convert_to_float)

        rows_with_float = data.applymap(is_float).any(axis=1)
        data = data[rows_with_float]

        colnames = data.columns.values.tolist()
        table = []
        for name in colnames:
            namelist = [name, ""]
            table.append(namelist)
        tp_tables[filename] = table
        tp_indata[filename] = data
        tp_pos[filename] = []
        tp_identifiers[filename] = []
        tp_refreshtable(window, table)
    return


def tp_rem(values, window, tp_paths, tp_tables):
    sure = sg.popup_yes_no("Remove data from list?")
    if sure == "Yes":
        tp_paths.remove(values["-tp_path-"])
        del tp_tables[values["-tp_path-"]]
        # del tp_data[values['-tp_path-']]
        if tp_paths:
            curr = tp_paths[0]
            tp_refreshtable(window, tp_tables[curr])
        else:
            curr = []
            tp_refreshtable(window, curr)
        window["-tp_path-"].Update(values=tp_paths, value=curr)
    else:
        pass
    return


def tp_defrem(values, window, tp_tables):
    path = values["-tp_path-"]
    selected = values["-tp_table-"]
    table = tp_tables[path]
    for index in sorted(selected, reverse=True):
        del table[index]
    tp_tables[path] = table
    window["-tp_table-"].Update(values=tp_tables[path])
    return


def tp_defkeep(values, window, tp_tables):
    path = values["-tp_path-"]
    table = tp_tables[path]
    for pos in values["-tp_table-"]:
        table[pos][1] = "[KEEP]"
    tp_tables[path] = table
    window["-tp_table-"].Update(values=tp_tables[path])
    return


def tp_defcon(values, window, tp_tables):
    if values["-tp_table-"]:
        path = values["-tp_path-"]
        table = tp_tables[path]
        value = simpledialog.askstring("Set Condition", "Condition Name")
        if value:
            for i in values["-tp_table-"]:
                table[i][1] = value
            tp_tables[path] = table
    else:
        messagebox.showerror("Error", "Select (a) sample(s).")
    window["-tp_table-"].Update(values=tp_tables[values["-tp_path-"]])
    return


def tp_defident(values, window, tp_tables, tp_pos, tp_identifiers):
    pos = values["-tp_table-"]
    if pos:
        if len(pos) > 1:
            messagebox.showerror("Error", "Please set only one Identifier!")
        elif len(pos) == 1:
            path = values["-tp_path-"]
            table = tp_tables[path]
            if tp_pos[path]:
                table[tp_pos[path][0]][1] = ""
            tp_identifiers[path] = table[pos[0]][0]
            tp_pos[path] = pos
            table[pos[0]][1] = "[IDENTIFIER]"
            tp_tables[path] = table
            window["-tp_table-"].Update(values=tp_tables[values["-tp_path-"]])
        else:
            messagebox.showerror("Error", "No sample selected.")
    return tp_identifiers


def tp_export(tp_data, tp_info):
    export_folder = sg.popup_get_folder("Export Folder")
    if export_folder:
        experiment = simpledialog.askstring("Export", "Experiment Name: ")
        now = datetime.now()
        time = now.strftime("%Y%m%d%H%M%S")
        export_full = pd.DataFrame()
        for condition in tp_data:
            path = (
                export_folder
                + "/"
                + time
                + "_"
                + experiment
                + "_"
                + condition
                + ".txt"
            )
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
        for info in tp_info:
            export_full = pd.merge(
                export_full,
                tp_info[info],
                left_index=True,
                right_index=True,
                how="left",
            )
        path = (
            export_folder
            + "/"
            + time
            + "_"
            + experiment
            + "_"
            + "combined"
            + ".txt"
        )
        export_full.to_csv(
            path,
            header=True,
            index=True,
            index_label="Identifier",
            sep="\t",
            mode="a",
        )
    return


# ------------------------------------------------------------------------------


def check_markers(marker_sets):
    is_markers = True

    if marker_sets:
        for file in marker_sets:
            if (
                marker_sets[file]["identifier_col"] == "-"
                or marker_sets[file]["class_col"] == "-"
            ):
                is_markers = False
    else:
        is_markers = False

    # print(marker_sets[file])
    return is_markers


def refresh_markertable(window, values, marker_sets):
    file_list = []
    for markerfile in marker_sets:
        file_list.append(markerfile)
    window["-marker_list-"].Update(values=file_list)
    if file_list:
        window["-marker_list-"].Update(set_to_index=0)
        window["-marker_key-"].Update(
            values=marker_sets[file_list[0]]["table"].columns.tolist(),
            value=marker_sets[file_list[0]]["identifier_col"],
        )
        window["-marker_class-"].Update(
            values=marker_sets[file_list[0]]["table"].columns.tolist(),
            value=marker_sets[file_list[0]]["class_col"],
        )
    return


def refresh_markercols(window, values, marker_sets):
    try:
        window["-marker_key-"].Update(
            values=marker_sets[values["-marker_list-"][0]][
                "table"
            ].columns.tolist(),
            value=marker_sets[values["-marker_list-"][0]]["identifier_col"],
        )
        window["-marker_class-"].Update(
            values=marker_sets[values["-marker_list-"][0]][
                "table"
            ].columns.tolist(),
            value=marker_sets[values["-marker_list-"][0]]["class_col"],
        )
    except Exception:
        window["-marker_key-"].Update(values=[], value="-")
        window["-marker_class-"].Update(values=[], value="-")
    return


def marker_add(window, values, marker_sets):
    filename = sg.popup_get_file(
        "Select a new Marker List!",
        no_window=True,
        file_types=(
            ("Tab delimited Text", "*.txt"),
            ("Tab Separated Values", "*.tsv"),
        ),
    )
    if filename:
        marker_sets[filename] = {}
        # marker_sets[filename]['table'] = pd.read_csv(filename, sep = "\t", header = 0).apply(lambda x: x.astype(str))
        marker_sets[filename]["table"] = pd.read_csv(
            filename, sep="\t", header=0
        ).apply(lambda x: x.astype(str).str.upper())
        marker_sets[filename]["identifier_col"] = "-"
        marker_sets[filename]["class_col"] = "-"
        marker_sets[filename]["classes"] = []
        refresh_markertable(window, values, marker_sets)

        # window['-marker_test-'].Update(disabled = False)
        # window['-marker_profiles-'].Update(disabled = False)
        # window['-marker_remove-'].Update(disabled = False)
    return


def marker_remove(window, values, marker_sets):
    del marker_sets[values["-marker_list-"][0]]
    refresh_markertable(window, values, marker_sets)
    if not len(marker_sets) > 0:
        window["-marker_test-"].Update(disabled=True)
        window["-marker_profiles-"].Update(disabled=True)
        window["-marker_remove-"].Update(disabled=True)
    return marker_sets


def marker_setkey(values, marker_sets):
    marker_sets[values["-marker_list-"][0]]["identifier_col"] = values[
        "-marker_key-"
    ]
    return


def marker_setclass(values, marker_sets):
    marker_sets[values["-marker_list-"][0]]["class_col"] = values[
        "-marker_class-"
    ]
    marker_sets[values["-marker_list-"][0]]["classes"] = list(
        set(
            marker_sets[values["-marker_list-"][0]]["table"][
                values["-marker_class-"]
            ]
        )
    )
    marker_conv = create_conversion(marker_sets)
    return marker_conv


def create_conversion(marker_sets):
    marker_conv = {}
    for path in marker_sets:
        for classname in marker_sets[path]["classes"]:
            marker_conv[classname] = classname
    return marker_conv


def create_markerlist(
    marker_sets: dict[str, dict[str, Any]],
    marker_conv: dict[str, str | float],
    marker_params: dict[str, Any],
) -> pd.DataFrame:
    markerset = pd.DataFrame(columns=["name"])
    counter = 1
    for path in marker_sets:
        id_col = marker_sets[path]["identifier_col"]
        class_col = marker_sets[path]["class_col"]
        mset = marker_sets[path]["table"][[id_col, class_col]]
        mset.rename(
            columns={
                id_col: "name",
                class_col: f"class{counter}",
            },
            inplace=True,
        )
        for classname in marker_conv:
            mset.replace(
                {f"class{counter}": {classname: marker_conv[classname]}},
                inplace=True,
            )
            mset.replace(
                {f"class{counter}": {r"^\s*$": np.nan}},
                regex=True,
                inplace=True,
            )
            mset = mset[mset[f"class{counter}"].notna()]

        markerset = pd.merge(markerset, mset, on="name", how="outer")
        counter += 1

    markerset.set_index("name", inplace=True)

    if marker_params["what"] == "unite":
        pass
    elif marker_params["what"] == "intersect":
        markerset.dropna(inplace=True)
    else:
        raise ValueError("Invalid marker parameter")

    if marker_params["how"] == "majority":
        markerset_final = pd.DataFrame(
            markerset.mode(axis=1, dropna=True)[0]
        ).rename(columns={0: "class"})
    elif marker_params["how"] == "exclude":
        markerset_final = markerset.mode(axis=1, dropna=True).fillna(np.nan)
        if 1 in markerset_final.columns:
            markerset_final = pd.DataFrame(
                markerset_final[markerset_final[1].isnull()][0]
            ).rename(columns={0: "class"})
        else:
            markerset_final.rename(columns={0: "class"}, inplace=True)
    else:
        raise ValueError("Invalid marker parameter")
    return markerset_final


def session_save(session: SessionModel):
    filename = sg.popup_get_file(
        "Save Session",
        no_window=True,
        file_types=(("Numpy", "*.npy"),),
        save_as=True,
    )
    if filename:
        session.to_numpy(filename)


def session_open(window, values, filename, model: SessionModel):
    # filename = sg.popup_get_file('Open Session', no_window=True, file_types=(('Numpy', '*.npy'),))
    # if filename:

    # Update session data
    tmp_session = SessionModel.from_numpy(filename)
    for key, value in tmp_session.dict().items():
        setattr(model, key, value)

    if model.fract_paths:
        fract_refreshtable(window, model.fract_tables[model.fract_paths[0]])
        window["-fractionation_path-"].Update(
            values=model.fract_paths, value=model.fract_paths[0]
        )
    else:
        fract_refreshtable(window, [])
        window["-fractionation_path-"].Update(
            values=model.fract_paths, value=""
        )
    if model.fract_data["class"]:
        fract_buttons(window, True)
    else:
        fract_buttons(window, False)
    if model.tp_paths:
        tp_refreshtable(window, model.tp_tables[model.tp_paths[0]])
        window["-tp_path-"].Update(
            values=model.tp_paths, value=model.tp_paths[0]
        )
    else:
        tp_refreshtable(window, [])
        window["-tp_path-"].Update(values=model.tp_paths, value="")
    if model.tp_data:
        tp_buttons(window, True)
    else:
        tp_buttons(window, False)

    if model.marker_sets:
        refresh_markertable(window, values, model.marker_sets)

        event, values = window.read(timeout=50)
        # marker_setkey(values, marker_sets)
        # marker_setclass(values, marker_sets)
        refresh_markercols(window, values, model.marker_sets)

    # if marker_list.empty:
    #     CCMPS_actions.enable_markersettings(window, True)
    #     window['-marker_test-'].Update(disabled = False)
    #     window['-marker_profiles-'].Update(disabled = False)
    #     window['-marker_remove-'].Update(disabled = False)
    # else:
    #     CCMPS_actions.enable_markersettings(window, False)
    #     window['-marker_test-'].Update(disabled = True)
    #     window['-marker_profiles-'].Update(disabled = True)
    #     window['-marker_remove-'].Update(disabled = True)

    # if fract_data['class'] and not marker_list.empty:
    #     # print('positive')
    #     window['-classification_MOP-'].Update(disabled = False)
    #     #window['-classification_SVM-'].Update(disabled = False)
    # else:
    #     # print('negative')
    #     window['-classification_MOP-'].Update(disabled = True)
    # window['-classification_SVM-'].Update(disabled = True)

    # window['-marker_fractkey-'].Update(values = ['[IDENTIFIER]'] + list(fract_info))

    # if fract_data['class']:
    #     window['-classification_MOP-'].Update(disabled = False)
    # else:
    #     window['-classification_MOP-'].Update(disabled = True)

    # if learning_xyz:
    #     window['-classification_statistics-'].Update(disabled = False)
    # else:
    #     window['-classification_statistics-'].Update(disabled = True)

    # if results:
    #     #window['-classification_comparison-'].Update(disabled = False)
    #     window['-status_statistics-'].Update('done!')
    #     window['-export_statistics-'].Update(disabled = False)
    # else:
    #     #window['-classification_comparison-'].Update(disabled = True)
    #     window['-status_statistics-'].Update('missing')
    #     window['-export_statistics-'].Update(disabled = True)

    # if comparison:
    #     window['-status_comparison-'].Update('done!')
    #     window['-export_comparison-'].Update(disabled = False)
    # else:
    #     window['-status_comparison-'].Update('missing')
    #     window['-export_comparison-'].Update(disabled = True)


# ------------------------------------------------------------------------------


# def convert_markers(markers, conversion, mode):
#     markerset = pd.DataFrame(columns = ['name'])
#     counter = 1
#     for path in markers:
#         mset = markers[path]['table'][[ markers[path]['identifier_col'] , markers[path]['class_col']         ]]
#         mset.rename(columns = {markers[path]['identifier_col'] : 'name'  , markers[path]['class_col'] : 'class'+str(counter)}, inplace = True)
#         for classname in conversion:
#             mset['class'+str(counter)].replace({classname : conversion[classname]}, inplace = True)
#             mset['class'+str(counter)].replace(r'^\s*$', np.nan, regex=True, inplace = True)
#             mset = mset[mset['class'+str(counter)].notna()]
#         markerset = pd.merge(markerset, mset, on = 'name', how = 'outer')
#         counter +=1
#     markerset.set_index('name', inplace = True)
#     if mode[0] == 'unite':
#         pass
#     elif mode[0] == 'intersect':
#         markerset.dropna(inplace = True)
#     if mode [1] == 'majority':
#         markerset_final = pd.DataFrame(markerset.mode(axis = 1, dropna = True)[0]).rename(columns = {0 : 'class'})
#     if mode [1] == 'exclude':
#         markerset_final = markerset.mode(axis = 1, dropna = True).fillna(np.nan)
#         if 1 in markerset_final.columns:
#             markerset_final = pd.DataFrame(markerset_final[markerset_final[1].isnull()][0]).rename(columns = {0 : 'class'})
#         else:
#             markerset_final.rename(columns = {0 : 'class'}, inplace = True)
#     return markerset_final


def create_markerprofiles(fract_data, key, fract_info, marker_list):
    profiles = {}
    profiles_vis = {}
    for condition in fract_data["class"]:
        profiles[condition] = copy.deepcopy(fract_data["class"][condition])
    for condition in fract_data["vis"]:
        profiles_vis[condition] = copy.deepcopy(fract_data["vis"][condition])

    fract_marker = {}
    fract_marker_vis = {}
    fract_test = {}

    if key == "[IDENTIFIER]":
        for condition in profiles:
            fract_marker[condition] = pd.merge(
                profiles[condition],
                marker_list,
                left_index=True,
                right_index=True,
                how="left",
            ).dropna(subset=["class"])
            fract_test[condition] = pd.merge(
                profiles[condition],
                marker_list,
                left_index=True,
                right_index=True,
                how="left",
            )
            fract_test[condition] = fract_test[condition][
                fract_test[condition]["class"].isna()
            ]
        for condition in profiles_vis:
            fract_marker_vis[condition] = pd.merge(
                profiles_vis[condition],
                marker_list,
                left_index=True,
                right_index=True,
            )

            # profiles_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_index = True, right_index = True)
            # fract_marker_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_index = True, right_index = True, how = 'left').dropna(subset = ['class'])

    else:
        for condition in profiles:
            profiles[condition] = pd.merge(
                profiles[condition],
                fract_info[key].astype(str).applymap(str.upper),
                left_index=True,
                right_index=True,
            )
            # profiles[condition] = pd.merge(profiles[condition], fract_info[key].applymap(str.upper), left_index = True, right_index = True)

            # fract_info_upper = fract_info[key].applymap(str.upper)
            fract_marker[condition] = (
                pd.merge(
                    profiles[condition],
                    marker_list,
                    left_on=key,
                    right_index=True,
                    how="left",
                )
                .drop(key, axis=1)
                .dropna(subset=["class"])
            )
            fract_test[condition] = pd.merge(
                profiles[condition],
                marker_list,
                left_on=key,
                right_index=True,
                how="left",
            ).drop(key, axis=1)
            fract_test[condition] = fract_test[condition][
                fract_test[condition]["class"].isna()
            ]
        for condition in profiles_vis:
            profiles_vis[condition] = pd.merge(
                profiles_vis[condition],
                fract_info[key],
                left_index=True,
                right_index=True,
            )
            fract_marker_vis[condition] = (
                pd.merge(
                    profiles_vis[condition],
                    marker_list,
                    left_on=key,
                    right_index=True,
                    how="left",
                )
                .drop(key, axis=1)
                .dropna(subset=["class"])
            )

        # for condition in profiles:
        #     print('check1')
        #     profiles[condition] = pd.merge(profiles[condition], fract_info[key], left_index = True, right_index = True)
        #     print('check2')
        #     profiles_vis[condition] = pd.merge(profiles_vis[condition], fract_info[key], left_index = True, right_index = True)
        #     print('check3')
        #     #print(profiles_vis[condition])
        #     fract_marker[condition] = pd.merge(profiles[condition], marker_list, left_on = key, right_index = True, how = 'left').drop(key, axis = 1).dropna(subset = ['class'])
        #     print('check4')
        #     fract_marker_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_on = key, right_index = True, how = 'left').drop(key, axis = 1).dropna(subset = ['class'])
        #     print('check5')
        #     fract_test[condition] = pd.merge(profiles[condition], marker_list, left_on = key, right_index = True, how = 'left'). drop(key, axis = 1)
        #     print('check6')
        #     fract_test[condition] = fract_test[condition][fract_test[condition]['class'].isna()]

    classnames = {}
    for condition in profiles:
        classnames[condition] = []
        for classname in list(set(fract_marker[condition]["class"].tolist())):
            classnames[condition].append(classname)

    return fract_marker, fract_marker_vis, fract_test, classnames


# def create_markerprofiles (fract_data, key, fract_info, marker_list):
#     profiles = {}

#     if key == '[IDENTIFIER]':
#         for condition in fract_data['class']:
#             profiles[condition] = fract_data['class'][condition]
#     else:
#         for condition in fract_data['class']:
#             profiles[condition] = pd.merge(fract_data['class'][condition], fract_info[key], left_index = True, right_index = True, how = 'left').set_index(key)

#     fract_marker = {}
#     fract_test = {}
#     for condition in profiles:
#         fract_marker[condition] = pd.merge(profiles[condition], marker_list, left_index = True, right_index = True, how = 'left').dropna(subset = ['class'])
#         fract_test[condition] = pd.merge(profiles[condition], marker_list, left_index = True, right_index = True, how = 'left')
#         fract_test[condition] = fract_test[condition][fract_test[condition]['class'].isna()]


#     # for condition in fract_data['class']:
#     #     fract_profiles[condition] = fract_data['class'][condition].set_index(key)


#     return fract_marker, fract_test


def upscale(fract_marker, fract_std, key, fract_info, mode):
    stds = {}
    if not key == "[IDENTIFIER]":
        for condition in fract_std["class"]:
            stds[condition] = pd.merge(
                fract_std["class"][condition],
                fract_info[key],
                left_index=True,
                right_index=True,
                how="left",
            ).set_index(key)

    fract_marker_up = {}
    for condition in fract_marker:
        print("condition", condition)

        class_sizes = {}
        for classname in list(set(fract_marker[condition]["class"])):
            class_sizes[classname] = list(
                fract_marker[condition]["class"]
            ).count(classname)
        class_maxsize = max(class_sizes.values())

        fract_marker_up[condition] = fract_marker[condition]
        k = 1
        for classname in list(set(fract_marker[condition]["class"])):
            print("class", classname)

            data_class_temp = fract_marker[condition].loc[
                fract_marker[condition]["class"] == classname
            ]
            data_class = data_class_temp.drop(columns=["class"])

            class_difference = abs(class_maxsize - class_sizes[classname])
            # print('maxsize:', class_maxsize)
            # print('class size:', class_sizes[classname])
            # print('difference:', class_difference)

            if class_sizes[classname] > class_maxsize:
                ID_rnd = random.sample(
                    list(data_class.index), class_difference - 1
                )
                fract_marker_up[condition].drop(ID_rnd, inplace=True)

            if class_sizes[classname] < class_maxsize:
                class_up = pd.DataFrame(columns=data_class.columns)

                class_std = data_class.std(axis=0).to_frame().transpose()
                class_std_flat = class_std.values.flatten()
                # print(class_std)

                if mode == "noised":
                    for i in range(class_difference):
                        ID_rnd = random.choice(list(data_class.index))
                        name_up = "up_" + str(k) + "_" + ID_rnd
                        k += 1

                        std_rnd = stds[condition].loc[[ID_rnd]]
                        std_rnd = std_rnd[
                            ~std_rnd.index.duplicated(keep="first")
                        ]

                        profile_rnd = data_class.loc[[ID_rnd]]
                        profile_rnd = profile_rnd[
                            ~profile_rnd.index.duplicated(keep="first")
                        ]

                        list_noised = []

                        for j in range(len(profile_rnd.columns)):
                            col_val = profile_rnd.columns[j]
                            suffix = profile_rnd.columns[j][
                                profile_rnd.columns[j].rfind("_") + 1 :
                            ]
                            col_std = std_rnd.columns[
                                std_rnd.columns.str.endswith(suffix)
                            ]
                            sigma = 2 * std_rnd[col_std].iloc[0]
                            nv = np.random.normal(
                                profile_rnd[col_val][0], sigma
                            )
                            nv = 0.0 if nv < 0 else 1.0 if nv > 1 else nv[0]
                            list_noised.append(nv)

                        profile_noised = pd.DataFrame(
                            [list_noised], columns=list(profile_rnd.columns)
                        )
                        profile_noised.index = [name_up]
                        profile_noised["class"] = [classname]
                        class_up = class_up.append(profile_noised)
                    fract_marker_up[condition] = fract_marker_up[
                        condition
                    ].append(class_up)

                if mode == "average":
                    for i in range(class_difference):
                        ID_rnd_1 = random.choice(list(data_class.index))
                        ID_rnd_2 = random.choice(list(data_class.index))
                        ID_rnd_3 = random.choice(list(data_class.index))
                        name_up = (
                            "up_"
                            + str(k)
                            + "_"
                            + ID_rnd_1
                            + "_"
                            + ID_rnd_2
                            + "_"
                            + ID_rnd_3
                        )
                        k += 1

                        profile_rnd_1 = data_class.loc[[ID_rnd_1]]
                        profile_rnd_1 = profile_rnd_1[
                            ~profile_rnd_1.index.duplicated(keep="first")
                        ]
                        profile_rnd_2 = data_class.loc[[ID_rnd_2]]
                        profile_rnd_2 = profile_rnd_2[
                            ~profile_rnd_2.index.duplicated(keep="first")
                        ]
                        profile_rnd_3 = data_class.loc[[ID_rnd_3]]
                        profile_rnd_3 = profile_rnd_3[
                            ~profile_rnd_3.index.duplicated(keep="first")
                        ]

                        profile_av = (
                            pd.concat(
                                [profile_rnd_1, profile_rnd_2, profile_rnd_3]
                            )
                            .median(axis=0)
                            .to_frame()
                            .transpose()
                        )

                        profile_av.index = [name_up]
                        profile_av["class"] = [classname]
                        class_up = class_up.append(profile_av)
                    fract_marker_up[condition] = fract_marker_up[
                        condition
                    ].append(class_up)

                if mode == "noisedaverage":
                    for i in range(class_difference):
                        # print(i)
                        ID_rnd_1 = random.choice(list(data_class.index))
                        ID_rnd_2 = random.choice(list(data_class.index))
                        ID_rnd_3 = random.choice(list(data_class.index))
                        name_up = (
                            "up_"
                            + str(k)
                            + "_"
                            + ID_rnd_1
                            + "_"
                            + ID_rnd_2
                            + "_"
                            + ID_rnd_3
                        )
                        k += 1

                        # class_std = data_class.std(axis = 1).to_frame().transpose()

                        profile_rnd_1 = data_class.loc[[ID_rnd_1]]
                        profile_rnd_1 = profile_rnd_1[
                            ~profile_rnd_1.index.duplicated(keep="first")
                        ]
                        profile_rnd_2 = data_class.loc[[ID_rnd_2]]
                        profile_rnd_2 = profile_rnd_2[
                            ~profile_rnd_2.index.duplicated(keep="first")
                        ]
                        profile_rnd_3 = data_class.loc[[ID_rnd_3]]
                        profile_rnd_3 = profile_rnd_3[
                            ~profile_rnd_3.index.duplicated(keep="first")
                        ]

                        profile_av = (
                            pd.concat(
                                [profile_rnd_1, profile_rnd_2, profile_rnd_3]
                            )
                            .median(axis=0)
                            .to_frame()
                            .transpose()
                        )
                        # print(len(profile_av))
                        profile_av_flat = profile_av.values.flatten()

                        list_noised = []

                        # print(class_std)
                        #
                        # for j in range(len(class_std.columns)):
                        #     # sigma = 2*class_std[class_std.columns[j]]
                        #     nv = np.random.normal(profile_av_flat, 2* class_std_flat, size = profile_av.shape)
                        #     nv = 0. if nv < 0 else 1. if nv > 1 else nv[0]
                        #     # print(sigma)
                        nv = np.random.normal(
                            profile_av_flat,
                            2 * class_std_flat,
                            size=profile_av.shape,
                        )
                        nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))

                        # print(nv)

                        # values = np.where(values > 1, 1, np.where(values < 0, 0, values))

                        profile_noised = pd.DataFrame(
                            nv, columns=profile_av.columns
                        )
                        profile_noised.index = [name_up]
                        profile_noised["class"] = [classname]
                        class_up = class_up.append(profile_noised)

                        # profile_av.index = [name_up]
                        # profile_av['class'] = [classname]
                        # class_up = class_up.append(profile_av)
                    fract_marker_up[condition] = fract_marker_up[
                        condition
                    ].append(class_up)

        # noised_df = pd.DataFrame([noised_values], columns=df_mean.columns)

    return fract_marker_up, class_maxsize


def marker_mix(fract_marker_up):
    fract_mixed_up = {}

    for condition in fract_marker_up:
        class_list = list(set(list(fract_marker_up[condition]["class"])))
        combinations = [
            (a, b)
            for idx, a in enumerate(class_list)
            for b in class_list[idx + 1 :]
        ]

        fract_mixed_up[condition] = pd.DataFrame(
            columns=fract_marker_up[condition].drop(columns=["class"]).columns
        )
        for classname in class_list:
            fract_mixed_up[condition][classname] = 0.0

        cur = 1
        for comb in combinations:
            profiles_own = (
                fract_marker_up[condition]
                .copy()
                .loc[fract_marker_up[condition]["class"] == comb[0]]
                .drop(columns=["class"])
            )
            profiles_other = (
                fract_marker_up[condition]
                .copy()
                .loc[fract_marker_up[condition]["class"] == comb[1]]
                .drop(columns=["class"])
            )

            new_index = [
                f"{i}_{j}"
                for i, j in zip(profiles_own.index, profiles_other.index)
            ]
            mix_steps = [i / 10 for i in range(0, 11)]

            for part in mix_steps:
                new_index_part = [
                    f"{i + cur}_{value}" for i, value in enumerate(new_index)
                ]
                own_part = profiles_own.multiply(part)
                other_part = profiles_other.multiply(1 - part)

                own_part.index = new_index_part
                other_part.index = new_index_part

                # profiles_mixed = profiles_own.add(profiles_other, fill_value = 0)
                profiles_mixed = own_part + other_part

                for classname in class_list:
                    if classname == comb[0]:
                        profiles_mixed[classname] = part
                    elif classname == comb[1]:
                        profiles_mixed[classname] = 1 - part
                    else:
                        profiles_mixed[classname] = 0.0

                fract_mixed_up[condition] = pd.concat(
                    [fract_mixed_up[condition], profiles_mixed]
                )

                cur += len(profiles_mixed)

        # fract_marker_dummies = pd.get_dummies(fract_marker_up['Class'])
        # fract_mixed_up = pd.concat([])

        # df = pd.concat([df, class_dummies], axis=1).drop(columns=['Class'])

        # fract_mixed_up[condition] = pd.concat([fract_mixed_up[condition], fract_marker_up[condition]])
    return fract_mixed_up


# res = [(a, b) for idx, a in enumerate(test_list) for b in test_list[idx + 1:]]


# df_add = df1.add(df2, fill_value=0)

# df.loc[df['column_name'] == some_value]

# def marker_mix (fract_marker, class_maxsize):
#     mix_steps = [i/5 for i in range(1, 5)]
#     mix_steps_reverse = [1-i for i in mix_steps]

#     fract_marker_mixed = {}

#     for condition in fract_marker:
#         # print('mixing classes for', condition)

#         fract_marker_mixed[condition] = pd.DataFrame(columns = fract_marker[condition].columns)
#         classlist = list(set(fract_marker[condition]['class']))
#         for newclass in classlist:
#             fract_marker_mixed[condition][newclass] = 0

#         k = 1
#         for class_own in list(set(fract_marker[condition]['class'])):
#             classlist.remove(class_own)
#             data_class_temp_own = fract_marker[condition].loc[fract_marker[condition]['class'] == class_own]
#             data_class_own = data_class_temp_own.drop(columns = ['class'])

#             for class_other in classlist:
#                 data_class_temp_other = fract_marker[condition].loc[fract_marker[condition]['class'] == class_other]
#                 data_class_other = data_class_temp_other.drop(columns = ['class'])

#                 print('mixing classes for', condition, class_own, class_other)

#                 for i in range(int(0.1*class_maxsize)):
#                     ID_rnd_own = random.choice(list(data_class_own.index))
#                     ID_rnd_other = random.choice(list(data_class_other.index))
#                     # print('own ID:', ID_rnd_own)
#                     # print('other ID:', ID_rnd_other)

#                     profile_rnd_own = data_class_own.loc[[ID_rnd_own]]
#                     profile_rnd_own = profile_rnd_own[~profile_rnd_own.index.duplicated(keep = 'first')]
#                     profile_rnd_other = data_class_other.loc[[ID_rnd_other]]
#                     profile_rnd_other = profile_rnd_other[~profile_rnd_other.index.duplicated(keep = 'first')]

#                     for r in range(len(mix_steps)):
#                         ID_mixed = 'mixed_' + str(k) + '_' + str(mix_steps[r]) + '*' + ID_rnd_own + '_' + str(mix_steps_reverse[r]) + '*' + ID_rnd_other

#                         profile_rnd_own *= mix_steps[r]
#                         profile_rnd_other *= mix_steps_reverse[r]

#                         profile_rnd_own.index = [ID_mixed]
#                         profile_rnd_other.index = [ID_mixed]

#                         profile_rnd_combined = profile_rnd_own.add(profile_rnd_other, fill_value = 0)
#                         # profile_rnd_combined = profile_rnd_combined.div(profile_rnd_combined.sum(axis = 1), axis = 0)

#                         for newclass in list(set(fract_marker[condition]['class'])):
#                             profile_rnd_combined[newclass] = 0.

#                         # print(profile_rnd_combined)
#                         profile_rnd_combined[class_own] = [mix_steps[r]]
#                         profile_rnd_combined[class_other] = [mix_steps_reverse[r]]
#                         profile_rnd_combined.index = [ID_mixed]

#                         fract_marker_mixed[condition] = pd.merge(fract_marker_mixed[condition], profile_rnd_combined, left_index = True, right_index = True, how = 'outer')

#                         k +=1


#     return fract_marker_mixed


# concatenated_df = pd.concat([df1, df2, df3])

# # Calculate the median for each column
# median_values = concatenated_df.median()

# # Create a new dataframe with the median values
# median_df = pd.DataFrame(median_values).transpose()


def create_fullprofiles(fract_marker, fract_test, fract_marker_up):
    fract_full = {}
    fract_full_up = {}
    for condition in fract_test:
        fract_full[condition] = pd.concat(
            [fract_test[condition], fract_marker[condition]]
        )
        fract_full_up[condition] = pd.concat(
            [fract_test[condition], fract_marker_up[condition]]
        )

    # fract_full = pd.concat([fract_test, fract_marker])
    # fract_full_up = pd.concat([fract_test, fract_marker_up])

    return fract_full, fract_full_up


# suffix = sample[sample.rfind('_')+1:]

# df.columns[df.columns.str.endswith("_Emp")]


def enable_markersettings(window, is_enabled):
    for element in [
        "-marker_add-",
        "-marker_remove-",
        "-marker_key-",
        "-marker_class-",
        "-marker_fractkey-",
        "-marker_parameters-",
        "-marker_manage-",
        "-marker_accept-",
        "-marker_profiles-",
        "-marker_test-",
    ]:
        window[element].Update(disabled=not is_enabled)

    for element in ["-marker_reset-"]:
        window[element].Update(disabled=is_enabled)

    # window['-marker_reset-'].Update(disabled = is_enabled)
    if is_enabled:
        window["-status_marker-"].Update(value="missing", text_color="white")
    else:
        window["-status_marker-"].Update(
            value="ready!", text_color="dark green"
        )
    return


# data_up = data_learning
# n = 1
# for organelle in list(set(data_learning['class'])):
#     data_class_temp = data_learning.loc[data_learning['class'] == organelle]
#     data_class = data_class_temp.drop(columns=['class'])


#     #data_class = data_learning.loc[data_learning['class'] == organelle].drop(columns=['class'])


#     class_difference = abs(batch_target - class_sizes[organelle])
#     print(class_difference)

#     if class_sizes[organelle] > batch_target:
#         ID_rnd = random.sample(list(data_class.index), class_difference-1)
#         data_up.drop(ID_rnd, inplace = True)

#     if class_sizes[organelle] < batch_target:
#         class_up = pd.DataFrame(columns = data_class.columns)
#         for i in range(class_difference):
#             ID_rnd = random.choice(list(data_class.index))
#             name_up = 'up_' + str(n) + '_' + ID_rnd
#             n = n +1

#             std_rnd = data_std.loc[[ID_rnd]]

#             profile_rnd = data_class.loc[[ID_rnd]]
#             profile_rnd = profile_rnd[~profile_rnd.index.duplicated(keep='first')]

#             list_noised = []

#             for j in range(len(profile_rnd.columns)):
#                 col_val = profile_rnd.columns[j]
#                 col_std = std_rnd.columns[j]
#                 sigma = 0.5* std_rnd[col_std][0]
#                 nv = np.random.normal(profile_rnd[col_val][0], sigma)
#                 nv = 0. if nv < 0 else 1. if nv > 1 else nv
#                 list_noised.append(nv)

#             profile_noised = pd.DataFrame([list_noised], columns= list(profile_rnd.columns))
#             profile_noised.index = [name_up]
#             profile_noised['class'] = [organelle]
#             class_up = class_up.append(profile_noised)

#         data_up = data_up.append(class_up)

# class_sizes = {}
# for organelle in list(set(data_learning['class'])):
#     class_sizes[organelle] = list(data_learning['class']).count(organelle)
# class_maxsize = max(class_sizes.values())


def refresh_window(window: sg.Window, status: dict[str, bool]):
    for element in ["-fractionation_reset-", "-fractionation_summary-"]:
        window[element].Update(disabled=not status["fractionation_data"])
    for element in [
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
    ]:
        window[element].Update(disabled=status["fractionation_data"])

    for element in ["-tp_reset-", "-tp_summary-", "-tp_export-"]:
        window[element].Update(disabled=not status["tp_data"])
    for element in [
        "-tp_add-",
        "-tp_remove-",
        "-tp_edit_remove-",
        "-tp_edit_keep-",
        "-tp_edit_condition-",
        "-tp_edit_identifier-",
        "-tp_parameters-",
        "-tp_start-",
    ]:
        window[element].Update(disabled=status["tp_data"])

    for element in ["-marker_remove-", "-marker_manage-", "-marker_accept-"]:
        if status["marker_matched"]:
            window[element].Update(disabled=True)
        else:
            window[element].Update(disabled=not status["marker_file"])
    for element in ["-marker_test-", "-marker_profiles-"]:
        window[element].Update(disabled=not status["marker_file"])
    for element in ["-marker_reset-"]:
        window[element].Update(disabled=not status["marker_matched"])
    for element in ["-marker_add-", "-marker_parameters-", "-marker_preset-"]:
        window[element].Update(disabled=status["marker_matched"])

    for element in ["-statistic_import-"]:
        window[element].Update(disabled=status["comparison_global"])

    if status["fractionation_data"]:
        window["-status_fract-"].Update("ready", text_color="dark green")
    else:
        window["-status_fract-"].Update("none", text_color="black")
    if status["tp_data"]:
        window["-status_tp-"].Update("ready", text_color="dark green")
    else:
        window["-status_tp-"].Update("none", text_color="black")
    if status["lipidome_data"]:
        window["-status_fract_lipid-"].Update("ready", text_color="dark green")
    else:
        window["-status_fract_lipid-"].Update("none", text_color="black")
    if status["marker_matched"]:
        window["-status_marker-"].Update("ready", text_color="dark green")
    else:
        window["-status_marker-"].Update("none", text_color="black")
    if status["lipidome_total"]:
        window["-status_total_lipid-"].Update("ready", text_color="dark green")
    else:
        window["-status_total_lipid-"].Update("none", text_color="black")

    for element in ["-classification_MOP-"]:
        if status["fractionation_data"] and status["marker_matched"]:
            window[element].Update(disabled=status["training"])
        else:
            window[element].Update(disabled=True)

    for element in ["-classification_validation-", "-classification_reset-"]:
        window[element].Update(disabled=not status["training"])

    if status["training"]:
        for element in ["-statistic_predict-"]:
            window[element].Update(disabled=status["proteome_prediction"])
        for element in [
            "-statistic_export-",
            "-statistic_report-",
            "-statistic_reset-",
            "-statistic_heatmap-",
            "-statistic_distribution-",
        ]:
            window[element].Update(disabled=not status["proteome_prediction"])

        if status["proteome_prediction"]:
            for element in ["-global_run-"]:
                window[element].Update(disabled=status["comparison_global"])
            for element in [
                "-global_heatmap-",
                "-global_distance-",
                "-global_report-",
                "-global_reset-",
            ]:
                window[element].Update(
                    disabled=not status["comparison_global"]
                )

            if status["comparison_global"] and status["tp_data"]:
                for element in ["-class_run-"]:
                    window[element].Update(disabled=status["comparison_class"])
                for element in [
                    "-class_heatmap-",
                    "-class_reorganization-",
                    "-class_report-",
                    "-class_reset-",
                ]:
                    window[element].Update(
                        disabled=not status["comparison_class"]
                    )
            else:
                for element in [
                    "-class_run-",
                    "-class_heatmap-",
                    "-class_reorganization-",
                    "-class_report-",
                    "-class_reset-",
                ]:
                    window[element].Update(disabled=True)

                if status["lipidome_data"]:
                    for element in ["-lipidome_predict-"]:
                        window[element].Update(
                            disabled=status["lipidome_prediction"]
                        )
                    for element in [
                        "-lipidome_report-",
                        "-lipidome_reset-",
                        "-lipidome_heatmap-",
                        "-lipidome_reorganization-",
                        "-lipidome_density-",
                        "-lipidome_composition-",
                    ]:
                        window[element].Update(
                            disabled=not status["lipidome_prediction"]
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
                        window[element].Update(disabled=True)

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
                window[element].Update(disabled=True)

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
            window[element].Update(disabled=True)

    return
