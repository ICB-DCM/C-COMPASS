"""Fractionation summary display."""

# FIXME: this module is currently unused
import FreeSimpleGUI as sg


def FSD_exec(preparams, data_ways):
    class_string = ""
    for condition in data_ways["class"]:
        class_string = class_string + condition + ", "
    class_string = class_string[:-2]

    layout_SD = [
        [
            sg.Frame(
                layout=[
                    [sg.Text("pre-scaling:", text_color="black")],
                    [
                        sg.Text(
                            "\t"
                            + str(preparams["class"]["scale1"][0])
                            + " ("
                            + preparams["class"]["scale1"][1]
                            + ")",
                            text_color="dark grey",
                        )
                    ],
                    [
                        sg.Text(
                            "filtered by InnerCorrelation:", text_color="black"
                        )
                    ],
                    [
                        sg.Text(
                            "\t" + str(preparams["class"]["corrfilter"]),
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("combination:", text_color="black")],
                    [
                        sg.Text(
                            "\t" + str(preparams["class"]["combination"]),
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("baselines removed:", text_color="black")],
                    [
                        sg.Text(
                            "\t" + str(preparams["class"]["zeros"]),
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("post-scaling:", text_color="black")],
                    [
                        sg.Text(
                            "\t"
                            + str(preparams["class"]["scale2"][0])
                            + " ("
                            + preparams["class"]["scale2"][1]
                            + ")",
                            text_color="dark grey",
                        )
                    ],
                ],
                size=(200, 300),
                title="Classification",
            ),
            sg.Frame(
                layout=[
                    [sg.Text("pre-scaling:", text_color="black")],
                    [
                        sg.Text(
                            "\t"
                            + str(preparams["vis"]["scale1"][0])
                            + " ("
                            + preparams["vis"]["scale1"][1]
                            + ")",
                            text_color="dark grey",
                        )
                    ],
                    [
                        sg.Text(
                            "filtered by InnerCorrelation:", text_color="black"
                        )
                    ],
                    [
                        sg.Text(
                            "\t" + str(preparams["vis"]["corrfilter"]),
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("combination:", text_color="black")],
                    [
                        sg.Text(
                            "\t" + str(preparams["vis"]["combination"]),
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("baselines removed:", text_color="black")],
                    [
                        sg.Text(
                            "\t" + str(preparams["vis"]["zeros"]),
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("post-scaling:", text_color="black")],
                    [
                        sg.Text(
                            "\t"
                            + str(preparams["vis"]["scale2"][0])
                            + " ("
                            + preparams["vis"]["scale2"][1]
                            + ")",
                            text_color="dark grey",
                        )
                    ],
                ],
                size=(200, 300),
                title="Visualization",
            ),
            sg.Frame(
                layout=[
                    [sg.Text("min. replicates:", text_color="black")],
                    [
                        sg.Text(
                            "\t"
                            + str(preparams["global"]["minrep"][0])
                            + " ("
                            + preparams["global"]["minrep"][1]
                            + ")",
                            text_color="dark grey",
                        )
                    ],
                    [sg.Text("max. MissingValues:", text_color="black")],
                    [
                        sg.Text(
                            "\t"
                            + str(preparams["global"]["missing"][0])
                            + " ("
                            + preparams["global"]["missing"][1]
                            + ")",
                            text_color="dark grey",
                        )
                    ],
                    [
                        sg.Text(
                            "OuterCorrelations calculated:", text_color="black"
                        )
                    ],
                    [
                        sg.Text(
                            "\t" + str(preparams["global"]["outcorr"]),
                            text_color="dark grey",
                        )
                    ],
                ],
                size=(200, 300),
                title="Global",
            ),
        ],
        [
            sg.Text("Conditions: ", text_color="black"),
            sg.Text(class_string, text_color="dark grey"),
        ],
    ]

    window_SD = sg.Window("Summary", layout_SD, size=(640, 400))

    while True:
        event_SD, values_PPMS = window_SD.read()

        if event_SD == sg.WIN_CLOSED:
            window_SD.close()
            break
    window_SD.close()
    return
