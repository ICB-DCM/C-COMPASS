"""Show markers."""

import os
from pathlib import Path

import FreeSimpleGUI as sg
import matplotlib.pyplot as plt
import pandas as pd

from .marker_correlation_dialog import (
    draw_figure,
    update_class_info,
    update_figure,
)
from .visualize import plot_marker_profiles


def show_marker_profiles_dialog(fract_data, fract_info, marker_list, key):
    profiles_dict = {}
    class_info_dict = {}
    distinct_profiles_dict = {}

    for condition in fract_data["vis"]:
        data = pd.merge(
            fract_data["vis"][condition],
            fract_info[key],
            left_index=True,
            right_index=True,
            how="left",
        )
        data.set_index(key, inplace=True)

        median_classprofiles = {}
        distinct_profiles = {}

        classnames = list(set(marker_list["class"]))
        for classname in classnames:
            marker_class = marker_list[marker_list["class"] == classname]
            data_class = data[data.index.isin(marker_class.index)]
            median_classprofiles[classname] = data_class.median()
            distinct_profiles[classname] = data_class

        profiles_df = pd.DataFrame(median_classprofiles)
        profiles_dict[condition] = profiles_df
        class_info_dict[condition] = update_class_info(
            marker_list, classnames, data
        )
        distinct_profiles_dict[condition] = distinct_profiles

    condition = list(profiles_dict.keys())[0]
    layout = [
        [
            sg.Text("Select Condition:"),
            sg.Combo(
                list(profiles_dict.keys()),
                key="-condition-",
                enable_events=True,
                default_value=condition,
                readonly=True,
                size=(25, 1),
            ),
        ],
        [
            sg.Canvas(key="-CANVAS-", expand_x=True, expand_y=True),
            sg.Table(
                values=class_info_dict[condition],
                headings=["Class", "n"],
                key="-CLASSINFO-",
                col_widths=[20, 6],
                auto_size_columns=False,
                cols_justification=["l", "r"],
                num_rows=26,
            ),
        ],
        [sg.Button("Export all Conditions...", key="-EXPORT-", size=(20, 1))],
    ]

    window = sg.Window(
        "Marker profiles",
        layout,
        finalize=True,
        size=(1100, 520),
        modal=True,
        resizable=True,
    )

    # Initial drawing
    fig = plot_marker_profiles(profiles_dict[condition], title=condition)
    figure_agg = draw_figure(window["-CANVAS-"].TKCanvas, fig)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED:
            break

        if event == "-condition-":
            condition = values["-condition-"]
            fig = plot_marker_profiles(
                profiles_dict[condition], title=condition
            )
            figure_agg = update_figure(
                window["-CANVAS-"].TKCanvas, figure_agg, fig
            )
            window["-CLASSINFO-"].update(values=class_info_dict[condition])
        elif event == "-EXPORT-":
            folder_path = sg.popup_get_folder("Select Folder")
            if folder_path:
                Path(folder_path).mkdir(parents=True, exist_ok=True)
                # Save the main Excel file with all conditions and median profiles
                with pd.ExcelWriter(
                    os.path.join(folder_path, "markerprofiles_combined.xlsx")
                ) as writer:
                    for cond, df in profiles_dict.items():
                        df.to_excel(writer, sheet_name=cond)

                # Save individual Excel files for each condition with distinct profiles for each class
                for cond, distinct_profiles in distinct_profiles_dict.items():
                    with pd.ExcelWriter(
                        os.path.join(
                            folder_path, f"markerprofiles_{cond}.xlsx"
                        )
                    ) as writer:
                        for classname, df in distinct_profiles.items():
                            df.to_excel(writer, sheet_name=classname)

                # Save the plot
                for cond, df in profiles_dict.items():
                    fig = plot_marker_profiles(df, title=cond)
                    fig.savefig(
                        os.path.join(
                            folder_path, f"markerprofiles_{cond}.pdf"
                        ),
                        format="pdf",
                    )
                    plt.close(fig)

    window.close()
