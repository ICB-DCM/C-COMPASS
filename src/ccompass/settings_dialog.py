"""The C-COMPASS application settings dialog."""

import os

import FreeSimpleGUI as sg

from .core import AppSettings, app_name


def show_settings_dialog(settings: AppSettings) -> None:
    """Show the application settings dialog.

    :param settings: The application settings that may be modified.
    """
    layout = [
        [
            sg.Text("Number of processes:"),
            sg.Push(),
            sg.Spin(
                list(range(1, os.cpu_count() + 1)),
                initial_value=settings.max_processes,
                key="-NUM_PROCESSES-",
                size=(2, 1),
                tooltip="The maximum number of processes to use for parallel "
                "processing",
            ),
        ],
        [sg.VPush()],
        [sg.Button("Save"), sg.Button("Cancel")],
    ]

    # create the window
    window = sg.Window(
        f"{app_name} Settings", layout, modal=True, resizable=True
    )

    # event loop
    while True:
        event, values = window.read()

        if event in (sg.WIN_CLOSED, "Cancel"):
            break

        if event == "Save":
            settings.max_processes = int(values["-NUM_PROCESSES-"])
            settings.save()
            break

    window.close()
