"""Class manager dialog.

Dialog window for selecting and mapping annotations to classes.
"""

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd


def _refresh_conversion(
    conversion: dict[str, str | float], values: dict
) -> dict[str, str | float]:
    """Refresh the conversion dictionary with the new values from the GUI."""
    for o in conversion:
        conversion[o] = (
            values[f"--{o}_class--"] if values[f"--{o}--"] else np.nan
        )
    return conversion


def _create_class_manager_window(
    marker_conv: dict[str, str | float], num_names: int
) -> sg.Window:
    """Create the class manager window."""
    layout_column = [
        [sg.Text("Annotation", size=(25, 1)), sg.Text("Class", size=(20, 1))],
        [sg.Text("-" * 80)],
        *[
            [
                sg.Checkbox(
                    o,
                    default=not pd.isnull(marker_conv[o]),
                    enable_events=True,
                    size=(20, 5),
                    key=f"--{o}--",
                ),
                sg.InputText(
                    str(marker_conv[o]),
                    visible=not pd.isnull(marker_conv[o]),
                    size=(20, 5),
                    key=f"--{o}_class--",
                ),
            ]
            for o in marker_conv
        ],
    ]

    layout_CM = [
        [
            sg.Frame(
                layout=[
                    [
                        sg.Column(
                            layout=layout_column,
                            size=(380, 340),
                            scrollable=True,
                            vertical_scroll_only=True,
                        )
                    ]
                ],
                title="Classes",
                size=(400, 380),
            ),
            sg.Column(
                layout=[
                    [
                        sg.Text("Initial annotations: "),
                        sg.Text(str(len(marker_conv))),
                    ],
                    [
                        sg.Text("Used annotations: "),
                        sg.Text(str(num_names), key="-num_anno-"),
                    ],
                    [
                        sg.Button(
                            "Accept",
                            size=(15, 1),
                            enable_events=True,
                            button_color="dark green",
                            key="--accept--",
                        )
                    ],
                    [
                        sg.Button(
                            "Cancel",
                            size=(15, 1),
                            enable_events=True,
                            button_color="black",
                            key="--cancel--",
                        )
                    ],
                ],
                size=(200, 340),
            ),
        ]
    ]

    return sg.Window("Classes", layout_CM, size=(600, 400), modal=True)


def show_class_manager_dialog(
    marker_conv: dict[str, str | float],
) -> dict[str, str | float]:
    """Show the class managar dialog."""
    conv_old = marker_conv
    num_names = sum(not pd.isnull(v) for v in marker_conv.values())

    window = _create_class_manager_window(marker_conv, num_names)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "--cancel--":
            marker_conv = conv_old
            break

        for k in marker_conv:
            # checkbox clicked?
            if event == f"--{k}--":
                window[f"--{k}_class--"].Update(visible=values[f"--{k}--"])
                if values[f"--{k}--"]:
                    window[f"--{k}_class--"].Update(value=k)
                    num_names += 1
                else:
                    window[f"--{k}_class--"].Update(value=False)
                    num_names -= 1
                window["-num_anno-"].Update(value=str(num_names))

        if event == "--accept--":
            marker_conv = _refresh_conversion(marker_conv, values)
            break

    window.close()
    return marker_conv
