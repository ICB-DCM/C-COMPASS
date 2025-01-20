"""The neural network training parameters dialog."""

import copy

import FreeSimpleGUI as sg

from .core import NeuralNetworkParametersModel


def _create_window(nn_params: NeuralNetworkParametersModel) -> sg.Window:
    """Create the window for the neural network training parameters dialog."""
    isadam = "adam" in nn_params.optimizers
    isrms = "rmsprop" in nn_params.optimizers
    issgd = "sgd" in nn_params.optimizers

    if nn_params.class_loss == "mean_squared_error":
        loss_default = "mean squared error"
    elif nn_params.class_loss == "binary_crossentropy":
        loss_default = "binary cross-entropy"
    else:
        raise ValueError(f"Unknown loss function: {nn_params.class_loss}")

    layout = [
        [
            sg.Checkbox(
                "upsampling\t\t",
                default=nn_params.upsampling,
                key="--upsampling--",
                disabled=False,
                enable_events=True,
            ),
            sg.Combo(
                ["noised", "average", "noisedaverage"],
                default_value=nn_params.upsampling_method,
                key="--upsampling_method--",
                size=(23, 1),
                disabled=not nn_params.upsampling,
                enable_events=True,
                readonly=True,
            ),
        ],
        [
            sg.Text("\t\t\t\tnoise:\t"),
            sg.Spin(
                list(range(0, 6, 1)),
                initial_value=nn_params.upsampling_noise,
                size=(10, 1),
                key="--noise--",
                disabled=not nn_params.upsampling,
                enable_events=True,
                readonly=True,
                text_color="black",
            ),
            sg.Text("* std"),
        ],
        [
            sg.Checkbox(
                "SVM-filter",
                default=nn_params.svm_filter,
                key="--svm_filter--",
                disabled=False,
                enable_events=True,
            )
        ],
        [
            sg.Text("dividing-step for profile mixing:\t"),
            sg.Spin(
                ["-", 2, 4, 5, 10],
                initial_value=nn_params.mixed_part,
                size=(15, 1),
                key="--mixstep--",
                disabled=False,
                enable_events=True,
                readonly=True,
                text_color="black",
            ),
        ],
        [sg.HSep()],
        [
            sg.Text("NN optimization method:\t\t"),
            sg.Combo(
                ["long", "short"],
                default_value=nn_params.NN_optimization,
                size=(25, 1),
                key="--optimization_method--",
                disabled=False,
                enable_events=True,
                readonly=True,
            ),
        ],
        [
            sg.Text("dense layer activation:\t\t"),
            sg.Combo(
                ["relu", "leakyrelu"],
                default_value=nn_params.NN_activation,
                size=(25, 1),
                key="--dense_activation--",
                disabled=False,
                enable_events=True,
                readonly=True,
            ),
        ],
        [
            sg.Text("output activation:\t\t\t"),
            sg.Combo(
                ["linear", "sigmoid", "softmax"],
                default_value=nn_params.class_activation,
                size=(25, 1),
                key="--out_activation--",
                disabled=False,
                enable_events=True,
                readonly=True,
            ),
        ],
        [
            sg.Text("loss function:\t\t\t"),
            sg.Combo(
                ["mean squared error", "binary cross-entropy"],
                default_value=loss_default,
                size=(25, 1),
                key="--loss--",
                disabled=False,
                enable_events=True,
                readonly=True,
            ),
        ],
        [
            sg.Text("optimizers:\t\t\t"),
            sg.Checkbox(
                "adam",
                default=isadam,
                key="--adam--",
                disabled=False,
                enable_events=True,
            ),
            sg.Checkbox(
                "rmsprop",
                default=isrms,
                key="--rmsprop--",
                disabled=False,
                enable_events=True,
            ),
            sg.Checkbox(
                "sgd",
                default=issgd,
                key="--sgd--",
                disabled=False,
                enable_events=True,
            ),
        ],
        [
            sg.Text("max. training epochs:\t\t\t"),
            sg.Spin(
                list(range(10, 110, 10)),
                initial_value=nn_params.NN_epochs,
                size=(15, 1),
                key="--epochs--",
                disabled=False,
                enable_events=True,
                readonly=True,
                text_color="black",
            ),
        ],
        [
            sg.Text("tuning runs:\t\t\t"),
            sg.Spin(
                list(range(1, 11, 1)),
                initial_value=nn_params.rounds,
                size=(15, 1),
                key="--rounds--",
                disabled=False,
                enable_events=True,
                readonly=True,
                text_color="black",
            ),
        ],
        [
            sg.Text("ensemble runs:\t\t\t"),
            sg.Spin(
                list(range(5, 55, 5)),
                initial_value=nn_params.subrounds,
                size=(15, 1),
                key="--subrounds--",
                disabled=False,
                enable_events=True,
                readonly=True,
                text_color="black",
            ),
        ],
        [sg.HSep()],
        [
            sg.Text("TP filter cut-off:\t\t\t"),
            sg.Spin(
                list(range(0, 100, 1)),
                initial_value=nn_params.reliability,
                size=(15, 1),
                key="--filter_threshold--",
                disabled=False,
                enable_events=True,
                readonly=True,
                text_color="black",
            ),
        ],
        [sg.HSep()],
        [
            sg.Button(
                "Accept",
                key="--NNP_accept--",
                disabled=False,
                enable_events=True,
                button_color="dark green",
            ),
            sg.Button(
                "Cancel",
                key="--NNP_cancel--",
                disabled=False,
                enable_events=True,
                button_color="black",
            ),
        ],
    ]

    return sg.Window("NN Parameters", layout, size=(470, 420), modal=True)


def show_dialog(
    params_old: NeuralNetworkParametersModel,
) -> NeuralNetworkParametersModel:
    nn_params = copy.deepcopy(params_old)
    window = _create_window(nn_params)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "--NNP_cancel--":
            nn_params = params_old
            break

        if event == "--NNP_accept--":
            break

        if event == "--upsampling--":
            nn_params.upsampling = values["--upsampling--"]
            window["--upsampling_method--"].Update(
                disabled=not values["--upsampling--"]
            )
            window["--noise--"].Update(disabled=not values["--upsampling--"])
        elif event == "--upsampling_method--":
            nn_params.upsampling_method = values["--upsampling_method--"]
        elif event == "--noise--":
            nn_params.upsampling_noise = values["--noise--"]
        elif event == "--svm_filter--":
            nn_params.svm_filter = values["--svm_filter--"]
        elif event == "--mixstep--":
            nn_params.mixed_part = values["--mixstep--"]
        elif event == "--optimization_method--":
            nn_params.NN_optimization = values["--optimization_method--"]
        elif event == "--dense_activation--":
            nn_params.NN_activation = values["--dense_activation--"]
        elif event == "--out_activation--":
            nn_params.class_activation = values["--out_activation--"]
        elif event == "--loss--":
            if values["--loss--"] == "mean squared error":
                nn_params.class_loss = "mean_squared_error"
            elif values["--loss--"] == "binary cross-entropy":
                nn_params.class_loss = "binary_crossentropy"
        elif event == "--adam--":
            if values["--adam--"] == True:
                nn_params.optimizers.append("adam")
            elif values["--adam--"] == False:
                nn_params.optimizers.remove("adam")
        elif event == "--rmsprop--":
            if values["--rmsprop--"] == True:
                nn_params.optimizers.append("rmsprop")
            elif values["--rmsprop--"] == False:
                nn_params.optimizers.remove("rmsprop")
        elif event == "--sgd--":
            if values["--sgd--"] == True:
                nn_params.optimizers.append("sgd")
            elif values["--sgd--"] == False:
                nn_params.optimizers.remove("sgd")
        elif event == "--epochs--":
            nn_params.NN_epochs = values["--epochs--"]
        elif event == "--rounds--":
            nn_params.rounds = values["--rounds--"]
        elif event == "--subrounds--":
            nn_params.subrounds = values["--subrounds--"]
        elif event == "--filter_threshold--":
            nn_params.reliability = values["--filter_threshold--"]

    window.close()
    return nn_params
