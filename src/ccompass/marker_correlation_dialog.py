"""Test marker"""

import os

import FreeSimpleGUI as sg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def draw_figure(canvas: sg.Canvas, figure: plt.Figure) -> FigureCanvasTkAgg:
    """Draw a figure on a canvas."""
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)
    return figure_canvas_agg


def update_figure(
    canvas: sg.Canvas, figure_agg: FigureCanvasTkAgg, figure: plt.Figure
) -> FigureCanvasTkAgg:
    """Update a figure on a canvas."""
    figure_agg.get_tk_widget().forget()
    plt.close("all")
    return draw_figure(canvas, figure)


def create_heatmap(dataframe: pd.DataFrame, title: str = None) -> plt.Figure:
    """Create correlation heatmap."""
    fig, ax = plt.subplots(figsize=(8, 8))  # Adjust the figure size as needed
    cax = ax.matshow(dataframe, cmap="coolwarm", vmin=-1, vmax=1)
    fig.colorbar(cax)
    ax.set_xticks(range(len(dataframe.columns)))
    ax.set_xticklabels(
        dataframe.columns, rotation=90, fontsize=8
    )  # Rotate x-axis labels 90 degrees
    ax.set_yticks(range(len(dataframe.index)))
    ax.set_yticklabels(dataframe.index, fontsize=8)
    plt.subplots_adjust(
        top=0.8, bottom=0.1, left=0.2
    )  # Adjust the top and bottom margins
    if title:
        plt.title(title)
    return fig


def update_class_info(
    marker_list: pd.DataFrame, classnames: list[str], data: pd.DataFrame
) -> list[tuple[str, int]]:
    """Compute the number of markers in each class."""
    class_info = []
    for classname in classnames:
        count = data[
            data.index.isin(
                marker_list[marker_list["class"] == classname].index
            )
        ].shape[0]
        class_info.append(
            (classname, count),
        )
    return class_info


def _create_window(
    condition: str, correlation_matrices, class_info_dict
) -> sg.Window:
    """Create the window for correlation plots."""
    layout = [
        [
            sg.Text("Select Condition:"),
            sg.Combo(
                list(correlation_matrices.keys()),
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
                col_widths=[20, 10],
                auto_size_columns=False,
                justification="left",
                num_rows=35,
            ),
        ],
        [sg.Button("Export all Conditions...", key="-EXPORT-", size=(20, 1))],
    ]

    window = sg.Window(
        "Marker correlations",
        layout,
        finalize=True,
        size=(900, 650),
        modal=True,
        resizable=True,
    )

    return window


def compute_correlation_and_class_info(
    fract_data, fract_info, marker_list, key
) -> tuple[dict[str, pd.DataFrame], dict[str, list[tuple[str, int]]]]:
    correlation_matrices = {}
    class_info_dict = {}

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

        classnames = list(set(marker_list["class"]))
        for classname in classnames:
            marker_class = marker_list[marker_list["class"] == classname]
            data_class = data[data.index.isin(marker_class.index)]
            median_classprofiles[classname] = data_class.median().to_numpy()

        correlation_matrix = np.zeros((len(classnames), len(classnames)))

        for i, class1 in enumerate(classnames):
            for j, class2 in enumerate(classnames):
                if i == j:
                    correlation_matrix[i, j] = 1.0
                else:
                    correlation_matrix[i, j] = np.corrcoef(
                        median_classprofiles[class1],
                        median_classprofiles[class2],
                    )[0, 1]

        correlation_df = pd.DataFrame(
            correlation_matrix, index=classnames, columns=classnames
        )
        correlation_matrices[condition] = correlation_df
        class_info_dict[condition] = update_class_info(
            marker_list, classnames, data
        )

    return correlation_matrices, class_info_dict


def show_marker_correlation_dialog(
    fract_data, fract_info, marker_list, key
) -> None:
    correlation_matrices, class_info_dict = compute_correlation_and_class_info(
        fract_data, fract_info, marker_list, key
    )
    # Initially selected condition
    condition = list(correlation_matrices.keys())[0]

    window = _create_window(condition, correlation_matrices, class_info_dict)

    # Initial drawing
    fig = create_heatmap(correlation_matrices[condition], title=condition)
    figure_agg = draw_figure(window["-CANVAS-"].TKCanvas, fig)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED:
            break

        if event == "-condition-":
            condition = values["-condition-"]
            fig = create_heatmap(
                correlation_matrices[condition], title=condition
            )
            figure_agg = update_figure(
                window["-CANVAS-"].TKCanvas, figure_agg, fig
            )
            window["-CLASSINFO-"].update(values=class_info_dict[condition])
        elif event == "-EXPORT-":
            folder_path = sg.popup_get_folder("Select Folder")
            if folder_path:
                for cond, df in correlation_matrices.items():
                    # Save the plot
                    fig = create_heatmap(df, title=cond)
                    fig.savefig(
                        os.path.join(folder_path, f"{cond}.pdf"), format="pdf"
                    )
                    plt.close(fig)

                # Save all data to an Excel file
                with pd.ExcelWriter(
                    os.path.join(folder_path, "correlation_matrices.xlsx")
                ) as writer:
                    for cond, df in correlation_matrices.items():
                        df.to_excel(writer, sheet_name=cond)

    window.close()
