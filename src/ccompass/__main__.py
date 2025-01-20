"""Package entry-point."""

import FreeSimpleGUI as sg

from .CCMPS import MainController
from .core import SessionModel


def main():
    """The entry point for the C-COMPASS application."""
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
