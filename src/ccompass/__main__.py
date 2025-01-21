"""Package entry-point."""

import argparse

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

    from .CCMPS import MainController
    from .core import SessionModel

    sg.theme("Dark Blue 3")

    model = SessionModel()
    controller = MainController(model=model)
    controller.run()

    # import dill
    # filepath = 'session.pkl'
    # dill.dump_session(filepath) # Save the session
    # dill.load_session(filepath) # Load the session


if __name__ == "__main__":
    main()
