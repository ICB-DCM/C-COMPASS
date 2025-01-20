"""
The fractionation parameters dialog.

This dialog allows the user to specify the parameters for the fractionation
data processing.
"""

import FreeSimpleGUI as sg

from .core import fract_default


def _set_combinations(params):
    combination_both = {}
    if params["class"]["combination"] == "median":
        combination_both["class_median"] = True
        combination_both["class_concat"] = False
        combination_both["class_separate"] = False
    elif params["class"]["combination"] == "concat":
        combination_both["class_median"] = False
        combination_both["class_concat"] = True
        combination_both["class_separate"] = False
    elif params["class"]["combination"] == "separate":
        combination_both["class_median"] = False
        combination_both["class_concat"] = False
        combination_both["class_separate"] = True
    if params["vis"]["combination"] == "median":
        combination_both["vis_median"] = True
        combination_both["vis_concat"] = False
        combination_both["vis_separate"] = False
    elif params["vis"]["combination"] == "concat":
        combination_both["vis_median"] = False
        combination_both["vis_concat"] = True
        combination_both["vis_separate"] = False
    elif params["vis"]["combination"] == "separate":
        combination_both["vis_median"] = False
        combination_both["vis_concat"] = False
        combination_both["vis_separate"] = True
    return combination_both


def _accept(values):
    params = {
        "class": {
            "scale1": [
                values["--class_scale1--"],
                values["--class_scale1_how--"],
            ],
            "corrfilter": values["--class_corrfilter--"],
            "scale2": [
                values["--class_scale2--"],
                values["--class_scale2_how--"],
            ],
            "zeros": values["--class_zeros--"],
        },
        "vis": {
            "scale1": [
                values["--vis_scale1--"],
                values["--vis_scale1_how--"],
            ],
            "corrfilter": values["--vis_corrfilter--"],
            "scale2": [values["--vis_scale2--"], values["--vis_scale2_how--"]],
            "zeros": values["--vis_zeros--"],
        },
        "global": {
            "missing": [values["--glob_missing--"], values["--glob_minval--"]],
            "minrep": [values["--glob_minrep--"], values["--glob_mincount--"]],
            "outcorr": values["--glob_outcorr--"],
        },
    }
    if values["--class_median--"]:
        params["class"]["combination"] = "median"
    elif values["--class_concat--"]:
        params["class"]["combination"] = "concat"
    elif values["--class_separate--"]:
        params["class"]["combination"] = "separate"
    if values["--vis_median--"]:
        params["vis"]["combination"] = "median"
    elif values["--vis_concat--"]:
        params["vis"]["combination"] = "concat"
    elif values["--vis_separate--"]:
        params["vis"]["combination"] = "separate"
    return params


def _reset(window, params):
    window["--class_scale1--"].Update(value=params["class"]["scale1"][0])
    window["--class_scale1_how--"].Update(value=params["class"]["scale1"][1])
    window["--class_corrfilter--"].Update(value=params["class"]["corrfilter"])
    window["--class_scale2--"].Update(value=params["class"]["scale2"][0])
    window["--class_scale2_how--"].Update(value=params["class"]["scale2"][1])
    window["--class_zeros--"].Update(value=params["class"]["zeros"])

    window["--vis_scale1--"].Update(value=params["vis"]["scale1"][0])
    window["--vis_scale1_how--"].Update(value=params["vis"]["scale1"][1])
    window["--vis_corrfilter--"].Update(value=params["vis"]["corrfilter"])
    window["--vis_scale2--"].Update(value=params["vis"]["scale2"][0])
    window["--vis_scale2_how--"].Update(value=params["vis"]["scale2"][1])
    window["--vis_zeros--"].Update(value=params["vis"]["zeros"])

    window["--glob_missing--"].Update(value=params["global"]["missing"][0])
    window["--glob_minval--"].Update(value=params["global"]["missing"][1])
    window["--glob_minrep--"].Update(value=params["global"]["minrep"][0])
    window["--glob_mincount--"].Update(value=params["global"]["minrep"][1])
    window["--glob_outcorr--"].Update(value=params["global"]["outcorr"])

    if params["class"]["combination"] == "median":
        window["--class_median--"].Update(value=True)
        window["--class_concat--"].Update(value=False)
        window["--class_separate--"].Update(value=False)
    elif params["class"]["combination"] == "concat":
        window["--class_median--"].Update(value=False)
        window["--class_concat--"].Update(value=True)
        window["--class_separate--"].Update(value=False)
    elif params["class"]["combination"] == "separate":
        window["--class_median--"].Update(value=False)
        window["--class_concat--"].Update(value=False)
        window["--class_separate--"].Update(value=True)
    if params["vis"]["combination"] == "median":
        window["--vis_median--"].Update(value=True)
        window["--vis_concat--"].Update(value=False)
        window["--vis_separate--"].Update(value=False)
    elif params["vis"]["combination"] == "concat":
        window["--vis_median--"].Update(value=False)
        window["--vis_concat--"].Update(value=True)
        window["--vis_separate--"].Update(value=False)
    elif params["vis"]["combination"] == "separate":
        window["--vis_median--"].Update(value=False)
        window["--vis_concat--"].Update(value=False)
        window["--vis_separate--"].Update(value=True)


def _create_window(params_old) -> sg.Window:
    combination_both = _set_combinations(params_old)

    layout_class = [
        [
            sg.Checkbox(
                "pre-scaling",
                key="--class_scale1--",
                disabled=False,
                enable_events=True,
                default=params_old["class"]["scale1"][0],
            ),
            sg.Combo(
                ["minmax", "area"],
                key="--class_scale1_how--",
                size=(10, 1),
                disabled=not params_old["class"]["scale1"][0],
                enable_events=False,
                readonly=True,
                default_value=params_old["class"]["scale1"][1],
            ),
        ],
        [
            sg.Checkbox(
                "exclude proteins from worst correlated replicate",
                key="--class_corrfilter--",
                disabled=False,
                enable_events=False,
                default=params_old["class"]["corrfilter"],
            )
        ],
        [
            sg.Checkbox(
                "median profile",
                key="--class_median--",
                disabled=False,
                enable_events=True,
                default=combination_both["class_median"],
                tooltip=" for consistent/reproducible replicates ",
            ),
            sg.Checkbox(
                "concatenated profiles",
                key="--class_concat--",
                disabled=False,
                enable_events=True,
                default=combination_both["class_concat"],
                tooltip=" for variations over replicate ",
            ),
            sg.Checkbox(
                "process separately",
                key="--class_separate--",
                disabled=False,
                enable_events=True,
                default=combination_both["class_separate"],
                tooltip=" for other purposes ",
            ),
        ],
        [
            sg.Checkbox(
                "post-scaling",
                key="--class_scale2--",
                disabled=False,
                enable_events=True,
                default=params_old["class"]["scale2"][0],
            ),
            sg.Combo(
                ["minmax", "area"],
                key="--class_scale2_how--",
                size=(10, 1),
                disabled=not params_old["class"]["scale2"][1],
                enable_events=False,
                readonly=True,
                default_value=params_old["class"]["scale2"][1],
            ),
        ],
        [
            sg.Checkbox(
                "remove baseline profiles (zeroes)",
                key="--class_zeros--",
                disabled=False,
                enable_events=False,
                default=params_old["class"]["zeros"],
            )
        ],
    ]
    layout_vis = [
        [
            sg.Checkbox(
                "pre-scaling",
                key="--vis_scale1--",
                disabled=False,
                enable_events=True,
                default=params_old["vis"]["scale1"][0],
            ),
            sg.Combo(
                ["minmax", "area"],
                key="--vis_scale1_how--",
                size=(10, 1),
                disabled=not params_old["vis"]["scale1"][0],
                enable_events=False,
                readonly=True,
                default_value=params_old["vis"]["scale1"][1],
            ),
        ],
        [
            sg.Checkbox(
                "exclude proteins from worst correlated replicate",
                key="--vis_corrfilter--",
                disabled=False,
                enable_events=False,
                default=params_old["vis"]["corrfilter"],
            )
        ],
        [
            sg.Checkbox(
                "median profile",
                key="--vis_median--",
                disabled=False,
                enable_events=True,
                default=combination_both["vis_median"],
                tooltip=" for consistent/reproducible replicates ",
            ),
            sg.Checkbox(
                "concatenated profiles",
                key="--vis_concat--",
                disabled=False,
                enable_events=True,
                default=combination_both["vis_concat"],
                tooltip=" for variations over replicate ",
            ),
            sg.Checkbox(
                "process separately",
                key="--vis_separate--",
                disabled=False,
                enable_events=True,
                default=combination_both["vis_separate"],
                tooltip=" for other purposes ",
            ),
        ],
        [
            sg.Checkbox(
                "post-scaling",
                key="--vis_scale2--",
                disabled=False,
                enable_events=True,
                default=params_old["vis"]["scale2"][0],
            ),
            sg.Combo(
                ["minmax", "area"],
                key="--vis_scale2_how--",
                size=(10, 1),
                disabled=not params_old["vis"]["scale2"][0],
                enable_events=False,
                readonly=True,
                default_value=params_old["vis"]["scale2"][1],
            ),
        ],
        [
            sg.Checkbox(
                "remove baseline profiles (zeroes)",
                key="--vis_zeros--",
                disabled=False,
                enable_events=False,
                default=params_old["class"]["zeros"],
            )
        ],
    ]

    layout_PPMS = [
        [
            sg.TabGroup(
                [
                    [
                        sg.Tab(" - Classification - ", layout_class),
                        sg.Tab(" - Visualization - ", layout_vis),
                    ]
                ],
                tab_location="top",
                tab_background_color="grey",
                size=(450, 225),
            ),
            sg.Frame(
                layout=[
                    [
                        sg.Checkbox(
                            "min. valid fractions:  ",
                            key="--glob_missing--",
                            disabled=False,
                            enable_events=True,
                            default=params_old["global"]["missing"][0],
                        )
                    ],
                    [
                        sg.Spin(
                            values=(list(range(1, 100))),
                            size=(10, 2),
                            key="--glob_minval--",
                            disabled=not params_old["global"]["missing"][0],
                            enable_events=False,
                            readonly=True,
                            initial_value=params_old["global"]["missing"][1],
                            text_color="black",
                        )
                    ],
                    [sg.HorizontalSeparator()],
                    [
                        sg.Checkbox(
                            "found in at least: ",
                            key="--glob_minrep--",
                            disabled=False,
                            enable_events=True,
                            default=params_old["global"]["minrep"][0],
                        )
                    ],
                    [
                        sg.Spin(
                            values=(list(range(1, 10))),
                            size=(10, 2),
                            key="--glob_mincount--",
                            disabled=not params_old["global"]["minrep"][0],
                            enable_events=False,
                            readonly=True,
                            initial_value=params_old["global"]["minrep"][1],
                            text_color="black",
                        ),
                        sg.Text(" replicates", text_color="light grey"),
                    ],
                    [sg.HorizontalSeparator()],
                    [
                        sg.Checkbox(
                            "calculate outer correlations",
                            key="--glob_outcorr--",
                            disabled=False,
                            enable_events=False,
                            default=params_old["global"]["outcorr"],
                        )
                    ],
                ],
                title="Global Parameters",
                size=(200, 220),
            ),
            sg.Column(
                layout=[
                    [
                        sg.Button(
                            "OK",
                            size=(10, 1),
                            key="--accept--",
                            button_color="darkgreen",
                        )
                    ],
                    [
                        sg.Button(
                            "Default",
                            size=(10, 1),
                            key="--default--",
                            button_color="grey",
                        )
                    ],
                    [
                        sg.Button(
                            "Cancel",
                            size=(10, 1),
                            key="--cancel--",
                            button_color="darkred",
                        )
                    ],
                ],
                size=(100, 250),
            ),
        ]
    ]

    return sg.Window(
        "Parameters for Pre-Processing",
        layout_PPMS,
        size=(800, 250),
        modal=True,
    )


def show_dialog(params_old):
    """Show dialog for fractionation parameters."""
    window = _create_window(params_old)

    while True:
        event, values = window.read()

        if event == "--cancel--" or event == sg.WIN_CLOSED:
            params = params_old
            window.close()
            break

        if event == "--accept--":
            params = _accept(values)
            break

        if event == "--class_scale1--":
            window["--class_scale1_how--"].Update(
                disabled=not values["--class_scale1--"]
            )
        elif event == "--vis_scale1--":
            window["--vis_scale1_how--"].Update(
                disabled=not values["--vis_scale1--"]
            )
        elif event == "--class_median--":
            window["--class_median--"].Update(value=True)
            window["--class_concat--"].Update(value=False)
            window["--class_separate--"].Update(value=False)
        elif event == "--class_concat--":
            window["--class_median--"].Update(value=False)
            window["--class_concat--"].Update(value=True)
            window["--class_separate--"].Update(value=False)
        elif event == "--class_separate--":
            window["--class_median--"].Update(value=False)
            window["--class_concat--"].Update(value=False)
            window["--class_separate--"].Update(value=True)
        elif event == "--vis_median--":
            window["--vis_median--"].Update(value=True)
            window["--vis_concat--"].Update(value=False)
            window["--vis_separate--"].Update(value=False)
            window["--glob_outcorr--"].Update(disabled=False)
        elif event == "--vis_concat--":
            window["--vis_median--"].Update(value=False)
            window["--vis_concat--"].Update(value=True)
            window["--vis_separate--"].Update(value=False)
            window["--glob_outcorr--"].Update(disabled=True, value=False)
        elif event == "--vis_separate--":
            window["--vis_median--"].Update(value=False)
            window["--vis_concat--"].Update(value=False)
            window["--vis_separate--"].Update(value=True)
            window["--glob_outcorr--"].Update(disabled=True, value=False)
        elif event == "--class_scale2--":
            window["--class_scale2_how--"].Update(
                disabled=not values["--class_scale2--"]
            )
        elif event == "--vis_scale2--":
            window["--vis_scale2_how--"].Update(
                disabled=not values["--vis_scale2--"]
            )
        elif event == "--glob_missing--":
            window["--glob_minval--"].Update(
                disabled=not values["--glob_missing--"]
            )
        elif event == "--glob_minrep--":
            window["--glob_mincount--"].Update(
                disabled=not values["--glob_minrep--"]
            )
        elif event == "--default--":
            _reset(window, fract_default())

    window.close()
    return params
