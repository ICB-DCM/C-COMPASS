"""Package entry-point."""

import argparse
import logging

from rich.logging import RichHandler

from . import app_name


def main():
    """The entry point for the C-COMPASS application."""
    from importlib.metadata import version

    parser = argparse.ArgumentParser(description="C-COMPASS application")
    parser.add_argument(
        "--version",
        action="version",
        version=f"{app_name} {version('ccompass')}",
    )
    parser.parse_args()

    launch_gui()


def launch_gui():
    """Launch the C-COMPASS GUI."""
    import FreeSimpleGUI as sg

    from .core import SessionModel
    from .main_gui import MainController

    # initialize logging
    log_format = "%(message)s"
    logging.basicConfig(
        level="WARNING",
        format=log_format,
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True)],
    )
    logger = logging.getLogger(__package__)
    logger.setLevel(logging.DEBUG)
    logger.info(f"Launching {app_name} GUI")

    # GUI theme
    sg.theme("Dark Blue 3")

    model = SessionModel()
    controller = MainController(model=model)
    controller.run()


if __name__ == "__main__":
    main()
