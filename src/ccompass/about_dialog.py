"""The Help->About dialog."""

from __future__ import annotations

import webbrowser
from importlib.metadata import Distribution, distributions, version

import FreeSimpleGUI as sg

from . import is_frozen
from .core import app_name, repository_url


def get_license_id(dist: Distribution) -> str:
    """Get the license identifier for a package."""
    return (
        dist.metadata.get("License", "").lstrip()
        or license_from_classifier(dist)
        or "Unknown"
    )


def get_license_text(dist: Distribution) -> str:
    """Get the license text for a package."""
    if license_text := dist.metadata.get("License"):
        return license_text

    if license_file := dist.metadata.get("License-File"):
        return dist.read_text(license_file)

    return license_from_classifier(dist) or "No license text available."


def license_from_classifier(dist: Distribution) -> str:
    """Get the license from the package classifier."""
    licenses = []
    for classifier in dist.metadata.get_all("Classifier", []):
        if classifier.startswith("License"):
            licenses.append(classifier.split("::", 1)[1].strip())
    return "; ".join(licenses)


def show_about_dialog():
    """Show the "about" dialog."""
    installed_packages = sorted(
        [
            (
                dist.metadata["Name"],
                dist.version,
                get_license_id(dist),
                get_license_text(dist),
            )
            for dist in distributions()
        ]
    )
    package_data = [
        [name, ver, lic] for name, ver, lic, _ in installed_packages
    ]

    default_font = sg.DEFAULT_FONT
    default_font_size = default_font[1]
    heading_font = (default_font[0], default_font_size + 4, "bold")

    # Currently, for pyinstaller executables, we have to manually include the
    #  metadata for each package we want to show here ...
    package_table_caption = (
        "Installed Packages (this list may be incomplete):"
        if is_frozen
        else "Installed Packages:"
    )

    layout = [
        [sg.Text(app_name, font=heading_font)],
        [sg.Text(f"Version: {version('ccompass')}")],
        [sg.Text(f"Website: {repository_url}", enable_events=True, key="URL")],
        [sg.HSeparator()],
        [sg.Text(package_table_caption)],
        [
            sg.Table(
                values=package_data,
                headings=["Package", "Version", "License"],
                display_row_numbers=False,
                auto_size_columns=True,
                num_rows=min(25, len(package_data)),
                cols_justification=["left", "right", "left"],
                key="PACKAGE_TABLE",
                enable_events=True,
                expand_x=True,
                expand_y=True,
            )
        ],
        [
            sg.Text("Package:"),
            sg.Text(
                "<select a package>", key="PACKAGE_NAME", enable_events=True
            ),
        ],
        [sg.Text("License:")],
        [
            sg.Multiline(
                size=(80, 20),
                key="LICENSE_TEXT",
                default_text="<select a package>",
                font=("Courier", 10),
                disabled=True,
                expand_x=True,
                expand_y=True,
            ),
        ],
        [sg.Push(), sg.Button("OK"), sg.Push()],
    ]
    window = sg.Window(f"About {app_name}", layout, resizable=True, modal=True)

    while True:
        event, values = window.read()
        if event == sg.WINDOW_CLOSED or event == "OK":
            break

        if event == "URL":
            webbrowser.open(repository_url)
        elif event == "PACKAGE_TABLE" and values.get("PACKAGE_TABLE"):
            selected_row = values["PACKAGE_TABLE"][0]
            package_name = installed_packages[selected_row][0]
            license_text = installed_packages[selected_row][3]
            homepage = None
            for dist in distributions():
                if dist.metadata["Name"] == package_name:
                    homepage = dist.metadata.get("Home-page")
                    break
            if homepage:
                package_name = f"{package_name} ({homepage})"
            window["PACKAGE_NAME"].update(package_name)
            window["LICENSE_TEXT"].update(license_text)
        elif event == "PACKAGE_NAME":
            selected_row = values["PACKAGE_TABLE"][0]
            package_name = installed_packages[selected_row][0]
            for dist in distributions():
                if dist.metadata["Name"] == package_name:
                    homepage = dist.metadata.get("Home-page")
                    if homepage:
                        webbrowser.open(homepage)
                    break

    window.close()
