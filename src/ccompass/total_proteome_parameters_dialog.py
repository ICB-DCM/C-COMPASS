"""
The total proteome parameters dialog.

This dialog allows the user to specify the parameters for the total proteome
data processing.
"""

import copy

import FreeSimpleGUI as sg


def _create_window(params_old) -> sg.Window:
    """Create the total proteome parameters dialog window."""
    is_normal = params_old["imputation"] == "normal"
    is_constant = params_old["imputation"] == "constant"

    layout = [
        [
            sg.Text("found in at least"),
            sg.Spin(
                values=(list(range(1, 10))),
                size=(10, 2),
                key="--tp_mincount--",
                disabled=False,
                enable_events=False,
                readonly=True,
                initial_value=params_old["minrep"],
                text_color="black",
            ),
        ],
        [
            sg.Text("Imputation of Missing Values:"),
            sg.Checkbox(
                "by normal distribution",
                key="--normal--",
                enable_events=True,
                disabled=False,
                default=is_normal,
            ),
            sg.Checkbox(
                "by constant (0)",
                key="--constant--",
                enable_events=True,
                disabled=False,
                default=is_constant,
            ),
        ],
        [
            sg.Button(
                "Accept",
                button_color="dark green",
                key="--TPPM_accept--",
                enable_events=True,
                disabled=False,
            ),
            sg.Button(
                "Cancel",
                key="--TPPM_cancel--",
                disabled=False,
                enable_events=True,
                button_color="black",
            ),
        ],
    ]

    return sg.Window("TP Parameters", layout, size=(500, 100), modal=True)


def show_dialog(params_old) -> dict:
    """
    Show the total proteome parameters dialog.
    """
    tp_params = copy.deepcopy(params_old)
    window = _create_window(tp_params)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "--TPPM_cancel--":
            tp_params = params_old
            break

        if event == "--TPPM_accept--":
            break

        if event == "--normal--":
            tp_params["imputation"] = "normal"
            window["--normal--"].update(value=True)
            window["--constant--"].update(value=False)
        elif event == "--constant--":
            tp_params["imputation"] = "constant"
            window["--constant--"].update(value=True)
            window["--normal--"].update(value=False)

    window.close()
    return tp_params
