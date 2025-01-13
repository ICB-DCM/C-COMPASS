"""C-COMPASS main window."""

from . import CCMPS_actions as action
from . import PPMS
from . import FDP
from . import TPP
from . import CM
from . import TM
from . import MOP
from . import MOA
from . import SM
from . import RP
import copy
import pickle
from tkinter import messagebox
import FreeSimpleGUI as sg
import pandas as pd
import os
from pathlib import Path

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"


def get_data_import_frame(fract_paths, tp_paths) -> sg.Frame:
    """Create the "Data Import" frame."""

    # The fractionation tab
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
                            sg.Tab(
                                " - Fractionation - ",
                                layout_fractionation,
                            ),
                            sg.Tab(
                                " - TotalProteomes - ",
                                layout_TP,
                                key="-total_tab-",
                            ),
                        ]
                    ],
                    tab_location="topleft",
                    tab_background_color="grey",
                    size=(600, 450),
                )
            ]
        ],
        title="Data Import",
        size=(620, 480),
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
    )


def main():
    """The entry point for the C-COMPASS application."""
    fract_paths, fract_tables, fract_data, fract_pos = action.resetinput()
    fract_data, fract_std, fract_intermediate, fract_info, fract_conditions = (
        action.reset_fract()
    )

    tp_paths, tp_tables, tp_data, tp_pos = action.resetinput()

    fract_indata, fract_identifiers = action.reset_infract()
    tp_indata, tp_identifiers = action.reset_intp()
    tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr = (
        action.reset_tp()
    )

    fract_preparams = PPMS.fract_default()
    tp_preparams = PPMS.tp_default()

    marker_sets = marker_conv = fract_full = fract_full_up = fract_marker = (
        fract_marker_vis
    ) = fract_marker_up = fract_mixed_up = fract_test = svm_marker = (
        svm_test
    ) = svm_metrics = learning_xyz = results = comparison = {}
    marker_params = {"how": "exclude", "what": "unite"}
    marker_list = pd.DataFrame
    NN_params = PPMS.NN_default()

    status = action.default_status()

    # status = {'fractionation_data' : False,
    #           'tp_data' : False,
    #           'lipidome_data' : False,
    #           'lipidome_total' : False,
    #           'marker_file' : False,
    #           'marker_matched' : False,
    #           'training' : False,
    #           'proteome_prediction' : False,
    #           'lipidome_prediction' : False,
    #           'comparison_global' : False,
    #           'comparison_class' : False}

    # tp_preparams = {}
    # marker_sets = {}
    # marker_conv = {}
    # fract_full = {}
    # fract_full_up = {}
    # fract_marker = {}
    # fract_marker_vis = {}
    # fract_marker_up = {}
    # fract_mixed_up = {}
    # fract_test = {}
    # svm_marker = {}
    # svm_test = {}
    # svm_metrics = {}
    # learning_xyz = {}
    # results = {}
    # comparison = {}

    sg.theme("Dark Blue 3")

    # The main menu
    menu_def = [
        ["Session", ["Save...", "Open...", "New", "Exit"]],
        ["Help", ["About...", "Open Website", "Manual"]],
    ]

    layout_CCMPS = [
        [sg.Menu(menu_def, tearoff=False)],
        [
            get_data_import_frame(fract_paths=fract_paths, tp_paths=tp_paths),
            get_spatial_prediction_frame(),
        ],
        [
            get_marker_selection_frame(),
            get_conditional_comparison_frame(),
        ],
    ]

    main_window = sg.Window("C-COMPASS", layout_CCMPS, size=(1260, 720))

    # The event loop
    while True:
        event, values_CCMPS = main_window.read()
        # action.refresh_window(window_CCMPS, status)

        # if status['fractionation_data']:
        #     window_CCMPS['-status_fract-'].Update('ready')
        # else:
        #     window_CCMPS['-status_fract-'].Update('missing')

        if event == "-fractionation_add-":
            action.fract_add(
                values_CCMPS,
                main_window,
                fract_paths,
                fract_tables,
                fract_indata,
                fract_pos,
                fract_identifiers,
            )
        elif event == "-fractionation_remove-":
            if values_CCMPS["-fractionation_path-"]:
                action.fract_rem(
                    values_CCMPS,
                    main_window,
                    fract_paths,
                    fract_tables,
                    fract_data,
                )
        elif event == "-fractionation_path-":
            action.fract_refreshtable(
                main_window,
                fract_tables[values_CCMPS["-fractionation_path-"]],
            )
        elif event == "-fractionation_edit_remove-":
            if values_CCMPS["-fractionation_table-"]:
                action.fract_defrem(values_CCMPS, main_window, fract_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-fractionation_edit_keep-":
            if values_CCMPS["-fractionation_table-"]:
                action.fract_defkeep(values_CCMPS, main_window, fract_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-fractionation_edit_condition-":
            if values_CCMPS["-fractionation_table-"]:
                action.fract_defcon(values_CCMPS, main_window, fract_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-fractionation_edit_replicate-":
            if values_CCMPS["-fractionation_table-"]:
                action.fract_defrep(values_CCMPS, main_window, fract_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-fractionation_edit_fractions-":
            if values_CCMPS["-fractionation_table-"]:
                action.fract_deffract(values_CCMPS, main_window, fract_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-fractionation_edit_identifier-":
            if values_CCMPS["-fractionation_table-"]:
                fract_identifiers = action.fract_defident(
                    values_CCMPS,
                    main_window,
                    fract_tables,
                    fract_pos,
                    fract_identifiers,
                )
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-fractionation_parameters-":
            fract_preparams = PPMS.PPMS_exec_fract(fract_preparams)
        elif event == "-fractionation_reset-":
            sure = sg.popup_yes_no(
                "Reset Fractionation Pre-Processing? You have to run it again to use your data."
            )
            if sure == "Yes":
                (
                    fract_data,
                    fract_std,
                    fract_intermediate,
                    fract_info,
                    fract_conditions,
                ) = action.reset_fract()
                action.fract_buttons(main_window, False)
                # window_CCMPS['-marker_fractkey-'].Update(values = ['[IDENTIFIER]'] + list(fract_info))
                marker_list = pd.DataFrame()
                fract_marker = {}
                fract_marker_vis = {}
                fract_test = {}
                fract_full = {}
                # action.enable_markersettings(window_CCMPS, True)
                status["fractionation_data"] = False

                main_window["-marker_fractkey-"].Update(
                    values=["[IDENTIFIER]"], value=""
                )
                status["marker_matched"] = False
                # window_CCMPS['-classification_SVM-'].Update(disabled = True)

                svm_marker = svm_test = svm_metrics = learning_xyz = (
                    results
                ) = comparison = {}
                status["marker_matched"] = status["marker_matched"] = status[
                    "training"
                ] = status["proteome_prediction"] = status[
                    "lipidome_prediction"
                ] = status["comparison_global"] = status[
                    "comparison_class"
                ] = False

            else:
                pass
        elif event == "-fractionation_start-":
            if fract_paths:
                (
                    fract_data,
                    fract_std,
                    fract_intermediate,
                    fract_info,
                    fract_conditions,
                ) = FDP.FDP_exec(
                    main_window,
                    fract_tables,
                    fract_preparams,
                    fract_identifiers,
                    fract_data,
                    fract_std,
                    fract_intermediate,
                    fract_info,
                    fract_conditions,
                    fract_indata,
                )
                main_window["-marker_fractkey-"].Update(
                    values=["[IDENTIFIER]"] + list(fract_info)
                )
                if fract_data["class"]:
                    status["fractionation_data"] = True
                #     action.fract_buttons(window_CCMPS, True)
            else:
                messagebox.showerror(
                    "No dataset!", "Please import a fractionation dataset."
                )
        elif event == "-fractionation_summary-":
            RP.RP_gradient_heatmap(fract_data)
            # FSD.FSD_exec(fract_preparams, fract_data)
        # if event_CCMPS == '-fractionation_export-':
        #     action.fract_export(values_CCMPS, fract_data, fract_info)

        elif event == "-tp_add-":
            action.tp_add(
                values_CCMPS,
                main_window,
                tp_paths,
                tp_tables,
                tp_indata,
                tp_pos,
                tp_identifiers,
            )
        elif event == "-tp_remove-":
            if values_CCMPS["-tp_path-"]:
                action.tp_rem(
                    values_CCMPS, main_window, tp_paths, tp_tables, tp_data
                )
        elif event == "-tp_path-":
            action.tp_refreshtable(
                main_window, tp_tables[values_CCMPS["-tp_path-"]]
            )
        elif event == "-tp_edit_remove-":
            if values_CCMPS["-tp_table-"]:
                action.tp_defrem(values_CCMPS, main_window, tp_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-tp_edit_keep-":
            if values_CCMPS["-tp_table-"]:
                action.tp_defkeep(values_CCMPS, main_window, tp_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-tp_edit_condition-":
            if values_CCMPS["-tp_table-"]:
                action.tp_defcon(values_CCMPS, main_window, tp_tables)
            else:
                messagebox.showerror("Error", "Select (a) row(s).")
        elif event == "-tp_edit_identifier-":
            if values_CCMPS["-tp_table-"]:
                tp_identifiers = action.tp_defident(
                    values_CCMPS,
                    main_window,
                    tp_tables,
                    tp_pos,
                    tp_identifiers,
                )
        elif event == "-tp_parameters-":
            tp_preparams = PPMS.PPMS_exec_TP(tp_preparams)
        elif event == "-tp_reset-":
            sure = sg.popup_yes_no(
                "Reset TotalProteome Pre-Processing? You have to run it again to use your data."
            )
            if sure == "Yes":
                tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr = (
                    action.reset_tp()
                )
                status["tp_data"] = False
                if status["comparison_class"]:
                    results, comparison = MOA.class_reset(results, comparison)
                    status["comparison_class"] = False
            else:
                pass
        elif event == "-tp_start-":
            if tp_paths:
                tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr = (
                    TPP.TPP_exec(
                        main_window,
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

                if tp_data:
                    status["tp_data"] = True
                    # action.tp_buttons(window_CCMPS, True)
            else:
                messagebox.showerror(
                    "No dataset!", "Please import a TP dataset."
                )
        elif event == "-tp_export-":
            action.tp_export(tp_data, tp_info)

        elif event == "-marker_add-":
            action.marker_add(main_window, values_CCMPS, marker_sets)
            event, values_CCMPS = main_window.read(timeout=50)
            if marker_sets:
                status["marker_file"] = True
            else:
                status["marker_file"] = False
        elif event == "-marker_remove-":
            try:
                action.marker_remove(main_window, values_CCMPS, marker_sets)
            except Exception:
                pass
            if marker_sets:
                status["marker_file"] = True
            else:
                status["marker_file"] = False
        elif event == "-marker_list-":
            action.refresh_markercols(main_window, values_CCMPS, marker_sets)
        elif event == "-marker_key-":
            action.marker_setkey(values_CCMPS, marker_sets)
        elif event == "-marker_class-":
            marker_conv = action.marker_setclass(values_CCMPS, marker_sets)

        elif event == "-marker_parameters-":
            marker_params = PPMS.PPMS_exec_marker(marker_params)

        elif event == "-marker_manage-":
            if action.check_markers(marker_sets):
                marker_conv = CM.CM_exec(marker_sets, marker_conv)
            else:
                messagebox.showerror(
                    "Error", "Please define key and class column."
                )

        elif event == "-marker_test-":
            if action.check_markers(marker_sets):
                try:
                    marker_list = action.create_markerlist(
                        marker_sets, marker_conv, marker_params
                    )
                    TM.TM_exec(
                        fract_data,
                        fract_info,
                        marker_list,
                        values_CCMPS["-marker_fractkey-"],
                    )
                except Exception:
                    messagebox.showerror(
                        "Error", "Something is wrong with your marker list."
                    )
            else:
                messagebox.showerror(
                    "Error", "Please define key and class column."
                )
        elif event == "-marker_profiles-":
            if action.check_markers(marker_sets):
                try:
                    marker_list = action.create_markerlist(
                        marker_sets, marker_conv, marker_params
                    )
                    SM.SM_exec(
                        fract_data,
                        fract_info,
                        marker_list,
                        values_CCMPS["-marker_fractkey-"],
                    )
                except Exception:
                    messagebox.showerror(
                        "Error", "Something is wrong with your marker list."
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
            elif fract_data["class"] == []:
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
                try:
                    # print('Starting try block')
                    marker_list = action.create_markerlist(
                        marker_sets, marker_conv, marker_params
                    )
                    # print('check1: marker_list created')
                    fract_marker, fract_marker_vis, fract_test, classnames = (
                        action.create_markerprofiles(
                            fract_data,
                            values_CCMPS["-marker_fractkey-"],
                            fract_info,
                            marker_list,
                        )
                    )
                    # print('check2: marker profiles created')
                    fract_full = MOP.create_fullprofiles(
                        fract_marker, fract_test
                    )
                    status["marker_matched"] = True
                    # print('check3: full profiles created')
                    # action.enable_markersettings(window_CCMPS, False)
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
            marker_list = pd.DataFrame()
            fract_marker = {}
            fract_marker_vis = {}
            fract_test = {}
            fract_full = {}
            svm_marker = svm_test = svm_metrics = learning_xyz = results = (
                comparison
            ) = {}
            status["marker_matched"] = status["training"] = status[
                "proteome_prediction"
            ] = status["lipidome_prediction"] = status[
                "comparison_global"
            ] = status["comparison_class"] = False
            # status['marker_matched'] = False

            # action.enable_markersettings(window_CCMPS, True)
            # window_CCMPS['-classification_MOP-'].Update(disabled = True)
            # window_CCMPS['-classification_SVM-'].Update(disabled = True)

        elif event == "-classification_parameters-":
            NN_params = PPMS.PPMS_exec_NN(NN_params)

        elif event == "-classification_MOP-":
            (
                learning_xyz,
                NN_params,
                fract_full_up,
                fract_marker_up,
                fract_mixed_up,
                fract_unmixed_up,
                svm_marker,
                svm_test,
                svm_metrics,
            ) = MOP.MOP_exec(
                fract_conditions,
                fract_full,
                fract_marker,
                fract_test,
                fract_std,
                fract_info,
                values_CCMPS["-marker_fractkey-"],
                NN_params,
            )
            # window_CCMPS['-classification_statistics-'].Update(disabled = False)
            # window_CCMPS['-status_comparison-'].Update('done!')
            status["training"] = True

        elif event == "-classification_reset-":
            svm_marker = svm_test = svm_metrics = learning_xyz = results = (
                comparison
            ) = {}
            status["training"] = status["proteome_prediction"] = status[
                "lipidome_prediction"
            ] = status["comparison_global"] = status["comparison_class"] = (
                False
            )

        elif event == "-statistic_predict-":
            results = MOA.stats_proteome(
                learning_xyz, NN_params, fract_data, fract_conditions
            )
            status["proteome_prediction"] = True

        elif event == "-statistic_export-":
            filename = sg.popup_get_file(
                "Export Statistics",
                no_window=True,
                file_types=(("Pickle", "*.pkl"),),
                save_as=True,
            )
            if filename:
                with open(filename, "wb") as file:
                    pickle.dump(results, file)

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
                        if condition in results:
                            messagebox.showerror(
                                "Error",
                                "There are already statistics for "
                                + condition
                                + " in your current session.",
                            )
                        else:
                            results[condition] = copy.deepcopy(
                                results_new[condition]
                            )
                    status["proteome_prediction"] = status["training"] = True
                except Exception:
                    messagebox.showerror("Error", "Incompatible file type!")

        elif event == "-statistic_report-":
            export_folder = sg.popup_get_folder("Statistics Report")
            if export_folder:
                for condition in results:
                    fname = Path(
                        export_folder, f"CCMPS_statistics_{condition}.xlsx"
                    )
                    selected_columns = [
                        col
                        for col in results[condition]["metrics"].columns
                        if col.startswith("fCC_")
                    ] + ["SVM_winner", "fNN_winner", "marker"]
                    df_out = results[condition]["metrics"][selected_columns]
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
                for comb in comparison:
                    fname = Path(
                        export_folder,
                        f"CCMPS_comparison_{comb[0]}_{comb[1]}.xlsx",
                    )
                    selected_columns = [
                        col
                        for col in comparison[comb]["metrics"].columns
                        if col.startswith("fRL_")
                    ] + ["fRLS", "DS", "P(t)_RLS"]
                    df_out = comparison[comb]["metrics"][selected_columns]
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
            export_folder = sg.popup_get_folder("Class-centric Changes Report")
            if export_folder:
                for condition in results:
                    fname = Path(
                        export_folder,
                        f"CCMPS_ClassComposition_{condition}.xlsx",
                    )
                    selected_columns = [
                        col
                        for col in results[condition]["metrics"].columns
                        if col.startswith("nCPA")
                    ] + ["TPA"]
                    df_out = results[condition]["metrics"][selected_columns]
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
                for comb in comparison:
                    fname = Path(
                        export_folder,
                        f"CCMPS_ClassComparison_{comb[0]}_{comb[1]}.xlsx",
                    )
                    selected_columns = [
                        col
                        for col in comparison[comb]["metrics"].columns
                        if col.startswith("nCFC_")
                    ]
                    df_out = comparison[comb]["metrics"][selected_columns]
                    df_out.columns = [
                        col.replace(
                            "nCFC_", "nCFC_normalizedClasscentricFoldChange_"
                        )
                        if col.startswith("nCFC_")
                        else col
                        for col in df_out.columns
                    ]
                    df_out.to_excel(fname, index=True)

        elif event == "-statistic_reset-":
            results = comparison = {}
            status["proteome_prediction"] = status[
                "lipidome_prediction"
            ] = status["comparison_global"] = status["comparison_class"] = (
                False
            )

        elif event == "-statistic_heatmap-":
            RP.RP_stats_heatmap(results)

        elif event == "-statistic_distribution-":
            RP.RP_stats_distribution(results)

        elif event == "-global_heatmap-":
            RP.RP_global_heatmap(comparison)

        elif event == "-global_distance-":
            RP.RP_global_distance(comparison)

        elif event == "-class_heatmap-":
            RP.RP_class_heatmap(results)

        elif event == "-class_reorganization-":
            RP.RP_class_reorganization(comparison)

        elif event == "-global_run-":
            comparison = MOA.global_comparison(results)
            status["comparison_global"] = True

        elif event == "-global_reset-":
            comparison = {}
            status["comparison_global"] = status["comparison_class"] = False

        elif event == "-class_run-":
            comparison = MOA.class_comparison(
                tp_data, fract_conditions, results, comparison
            )
            status["comparison_class"] = True

        elif event == "-class_reset-":
            results, comparison = MOA.class_reset(results, comparison)
            status["comparison_class"] = False

        # if event_CCMPS == '-classification_comparison-':
        #     # results = MOP_stats.comp_exec(learning_xyz, results)
        #     comparison = MOP_stats.comp_exec3('deep', results, learning_xyz)

        # if event_CCMPS == '-classification_comparison_rough-':
        #     comparison = MOP_stats.comp_exec3('rough', results, learning_xyz)

        elif event == "-export_statistics-":
            export_folder = sg.popup_get_folder("Export Statistics")
            if export_folder:
                for condition in results:
                    fname = Path(
                        export_folder, f"CCMPS_statistics_{condition}.tsv"
                    )
                    df_out = pd.merge(
                        fract_data["vis"][condition + "_median"],
                        results[condition]["metrics"],
                        left_index=True,
                        right_index=True,
                        how="outer",
                    )
                    for colname in fract_info:
                        df_out = pd.merge(
                            df_out,
                            fract_info[colname],
                            left_index=True,
                            right_index=True,
                            how="left",
                        )
                    df_out.to_csv(
                        fname, sep="\t", index=True, index_label="Identifier"
                    )
                    # results[condition]['metrics'].to_csv(fname, sep='\t', index=True, index_label='Identifier')

        elif event == "-export_comparison-":
            export_folder = sg.popup_get_folder("Export Statistics")
            if export_folder:
                for comb in comparison:
                    fname = Path(
                        export_folder,
                        f"CCMPS_comparison_{comb[0]}_{comb[1]}.tsv",
                    )
                    df_out = pd.DataFrame(
                        index=comparison[comb]["intersection_data"].index
                    )
                    df_out = pd.merge(
                        df_out,
                        comparison[comb]["metrics"],
                        left_index=True,
                        right_index=True,
                        how="left",
                    )
                    for colname in fract_info:
                        df_out = pd.merge(
                            df_out,
                            fract_info[colname],
                            left_index=True,
                            right_index=True,
                            how="left",
                        )
                    df_out.to_csv(
                        fname, sep="\t", index=True, index_label="Identifier"
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
            action.session_save(
                fract_paths,
                fract_tables,
                fract_pos,
                fract_indata,
                fract_data,
                fract_std,
                fract_intermediate,
                fract_identifiers,
                fract_info,
                fract_preparams,
                tp_paths,
                tp_tables,
                tp_pos,
                tp_indata,
                tp_data,
                tp_intermediate,
                tp_identifiers,
                tp_info,
                tp_preparams,
                marker_sets,
                marker_params,
                marker_conv,
                fract_conditions,
                fract_full,
                fract_full_up,
                fract_marker,
                fract_marker_vis,
                fract_marker_up,
                fract_mixed_up,
                fract_test,
                svm_marker,
                svm_test,
                svm_metrics,
                marker_list,
                learning_xyz,
                NN_params,
                results,
                comparison,
                values_CCMPS["-marker_fractkey-"],
                status,
            )
        elif event == "Open...":
            filename = sg.popup_get_file(
                "Open Session",
                no_window=True,
                file_types=(("Numpy", "*.npy"),),
            )
            if filename:
                (
                    fract_paths,
                    fract_tables,
                    fract_pos,
                    fract_indata,
                    fract_data,
                    fract_std,
                    fract_intermediate,
                    fract_identifiers,
                    fract_info,
                    fract_preparams,
                    tp_paths,
                    tp_tables,
                    tp_pos,
                    tp_indata,
                    tp_data,
                    tp_intermediate,
                    tp_identifiers,
                    tp_info,
                    tp_preparams,
                    marker_sets,
                    marker_params,
                    marker_conv,
                    fract_conditions,
                    svm_marker,
                    svm_test,
                    svm_metrics,
                    fract_full,
                    fract_full_up,
                    fract_marker,
                    fract_marker_up,
                    fract_marker_vis,
                    fract_mixed_up,
                    fract_test,
                    marker_list,
                    learning_xyz,
                    results,
                    NN_params,
                    comparison,
                    marker_fractkey,
                    status,
                ) = action.session_open(main_window, values_CCMPS, filename)
                # window_CCMPS['-marker_tpkey-'].Update(values = ['[IDENTIFIER]'] + tp_info.columns.tolist())
                main_window["-marker_fractkey-"].Update(
                    values=["[IDENTIFIER]"] + list(fract_info),
                    value=marker_fractkey,
                )
        elif event == "New":
            sure = sg.popup_yes_no(
                "Are you sure to close the session and start a new one?"
            )
            if sure == "Yes":
                fract_paths, fract_tables, fract_data, fract_pos = (
                    action.resetinput()
                )
                action.fract_clearinput(main_window)
                action.tp_clearinput(main_window)
                fract_indata, fract_identifiers = action.reset_infract()
                (
                    fract_data,
                    fract_std,
                    fract_intermediate,
                    fract_info,
                    fract_conditions,
                ) = action.reset_fract()
                tp_paths, tp_tables, tp_data, tp_pos = action.resetinput()
                tp_indata, tp_identifiers = action.reset_intp()
                tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr = (
                    action.reset_tp()
                )
                fract_preparams = PPMS.fract_default()
                tp_preparams = PPMS.tp_default()
                marker_sets = marker_conv = fract_full = fract_full_up = (
                    fract_marker
                ) = fract_marker_vis = fract_marker_up = fract_mixed_up = (
                    fract_test
                ) = svm_marker = svm_test = svm_metrics = learning_xyz = (
                    results
                ) = comparison = {}
                main_window["-marker_list-"].Update(values=[])
                main_window["-marker_key-"].Update(values=[])
                main_window["-marker_class-"].Update(values=[])
                marker_params = {"how": "exclude", "what": "unite"}
                marker_list = pd.DataFrame
                NN_params = PPMS.NN_default()

                status = action.default_status()

                main_window["-marker_fractkey-"].Update(
                    values=["[IDENTIFIER]"] + list(fract_info)
                )
            else:
                pass

        # -----------------------------------------------------------------------------------------------------------------------------

        elif event == sg.WIN_CLOSED or event == "Exit":
            break

        action.refresh_window(main_window, status)

    main_window.close()

    # import dill
    # filepath = 'session.pkl'
    # dill.dump_session(filepath) # Save the session
    # dill.load_session(filepath) # Load the session
