"""Package entry-point."""

import argparse
import logging
import os
import sys

from .core import app_name


def main():
    """The entry point for the C-COMPASS application."""
    from importlib.metadata import version

    parser = argparse.ArgumentParser(description="C-COMPASS application")
    parser.add_argument(
        "--version",
        action="version",
        version=f"{app_name} {version('ccompass')}",
    )
    parser.add_argument(
        "--test-run",
        action="store_true",
        help="For testing purposes only. "
        "Perform a test run of the application.",
    )
    args = parser.parse_args()

    if args.test_run:
        from ._testing import do_test_run

        do_test_run()
        return

    launch_gui()


def init_logging() -> logging.Logger:
    """Initialize logging."""
    from rich.console import Console
    from rich.logging import RichHandler

    console = Console(
        # Write to sys.__stdout__ instead of sys.stdout to avoid infinite
        #  recursion, because we will redirect sys.stdout to the logger later
        #  on.
        file=sys.__stdout__,
        # When running in PyCharm, `rich` can't determine the terminal width
        #  correctly, so we set it manually, since the default width (80) is
        #  too narrow.
        width=160 if os.environ.get("PYCHARM_HOSTED", False) else None,
    )
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

    logger.info(f"Launching {app_name} GUI")

    # GUI theme
    sg.theme("Dark Blue 3")

    model = SessionModel()
    controller = MainController(model=model)
    controller.run()


logger = init_logging()


if __name__ == "__main__":
    main()
