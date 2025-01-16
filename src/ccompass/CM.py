"""Class manager.

Window for mapping annotations to classes.
"""

import FreeSimpleGUI as sg
import numpy as np
import pandas as pd


def refresh_conversion(conversion, values):
    """Refresh the conversion dictionary with the new values from the GUI."""
    for o in conversion:
        conversion[o] = (
            values[f"--{o}_class--"] if values[f"--{o}--"] else np.nan
        )
    return conversion


def CM_exec(marker_conv: dict[str, str | float]) -> dict[str, str | float]:
    conv_old = marker_conv
    num_names = sum(not pd.isnull(v) for v in marker_conv.values())

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

    window_CM = sg.Window("Classes", layout_CM, size=(600, 400))

    while True:
        event_CM, values_CM = window_CM.read()
        for k in marker_conv:
            # checkbox clicked?
            if event_CM == f"--{k}--":
                window_CM[f"--{k}_class--"].Update(
                    visible=values_CM[f"--{k}--"]
                )
                if values_CM[f"--{k}--"]:
                    window_CM[f"--{k}_class--"].Update(value=k)
                    num_names += 1
                else:
                    window_CM[f"--{k}_class--"].Update(value=False)
                    num_names -= 1
                window_CM["-num_anno-"].Update(value=str(num_names))

        if event_CM == sg.WIN_CLOSED or event_CM == "--cancel--":
            marker_conv = conv_old
            break
        if event_CM == "--accept--":
            marker_conv = refresh_conversion(marker_conv, values_CM)
            break

    window_CM.close()
    return marker_conv
