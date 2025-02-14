"""Package entry-point."""

import argparse
import logging
import sys

from rich.console import Console
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


def init_logging() -> logging.Logger:
    """Initialize logging."""
    # Write to sys.__stdout__ instead of sys.stdout to avoid infinite
    #  recursion, because we will redirect sys.stdout to the logger later on.
    console = Console(file=sys.__stdout__)
    log_format = "%(message)s"
    logging.basicConfig(
        level="WARNING",
        format=log_format,
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, console=console)],
    )
    logger = logging.getLogger(__package__)
    logger.setLevel(logging.DEBUG)

    return logger


def launch_gui():
    """Launch the C-COMPASS GUI."""
    import FreeSimpleGUI as sg

    from .core import SessionModel
    from .main_gui import MainController

    logger = init_logging()
    logger.info(f"Launching {app_name} GUI")

    # GUI theme
    sg.theme("Dark Blue 3")

    model = SessionModel()
    controller = MainController(model=model)
    controller.run()


if __name__ == "__main__":
    main()
