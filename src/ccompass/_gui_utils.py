"""GUI utilities for internal use."""

from contextlib import contextmanager

import FreeSimpleGUI as sg


@contextmanager
def wait_cursor(window: sg.Window):
    """Context manager to set the cursor to a wait cursor."""
    try:
        window.set_cursor("watch")
        window.read(timeout=0)
        yield
    finally:
        window.set_cursor("arrow")
