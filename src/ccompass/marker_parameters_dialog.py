"""Dialog for setting marker processing parameters."""

import copy
from typing import Any

import FreeSimpleGUI as sg


def _create_window(marker_params: dict[str, Any]) -> sg.Window:
    if marker_params["how"] == "exclude":
        how_exclude = True
        how_majority = False
    elif marker_params["how"] == "majority":
        how_exclude = False
        how_majority = True
    else:
        raise ValueError("Invalid value for 'how' in marker parameters.")

    if marker_params["what"] == "unite":
        what_unite = True
        what_intersect = False
    elif marker_params["what"] == "intersect":
        what_unite = False
        what_intersect = True
    else:
        raise ValueError("Invalid value for 'what' in marker parameters.")

    layout = [
        [
            sg.Text("discrepancies:\t"),
            sg.Checkbox(
                "exclude\t",
                key="--PMM_exclude--",
                default=how_exclude,
                enable_events=True,
            ),
            sg.Checkbox(
                "take majority",
                key="--PMM_majority--",
                default=how_majority,
                enable_events=True,
            ),
        ],
        [
            sg.Text("selection:\t"),
            sg.Checkbox(
                "unite\t",
                key="--PMM_unite--",
                default=what_unite,
                enable_events=True,
            ),
            sg.Checkbox(
                "intersect",
                key="--PMM_intersect--",
                default=what_intersect,
                enable_events=True,
            ),
        ],
        [
            sg.Button(
                "OK",
                key="--PMM_accept--",
                disabled=False,
                enable_events=True,
                button_color="dark green",
            ),
            sg.Button(
                "Cancel",
                key="--PMM_cancel--",
                disabled=False,
                enable_events=True,
                button_color="black",
            ),
        ],
    ]

    return sg.Window("Marker Parameters", layout, size=(400, 100), modal=True)


def show_dialog(params_old: dict[str, Any]) -> dict[str, Any]:
    """Show dialog for setting marker processing parameters."""
    marker_params = copy.deepcopy(params_old)
    window = _create_window(marker_params)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "--PMM_cancel--":
            marker_params = copy.deepcopy(params_old)
            break

        if event == "--PMM_accept--":
            break

        if event == "--PMM_exclude--":
            marker_params["how"] = "exclude"
            window["--PMM_exclude--"].Update(value=True)
            window["--PMM_majority--"].Update(value=False)
        elif event == "--PMM_majority--":
            marker_params["how"] = "majority"
            window["--PMM_exclude--"].Update(value=False)
            window["--PMM_majority--"].Update(value=True)
        elif event == "--PMM_unite--":
            marker_params["what"] = "unite"
            window["--PMM_unite--"].Update(value=True)
            window["--PMM_intersect--"].Update(value=False)
        elif event == "--PMM_intersect--":
            marker_params["what"] = "intersect"
            window["--PMM_unite--"].Update(value=False)
            window["--PMM_intersect--"].Update(value=True)

    window.close()
    return marker_params
