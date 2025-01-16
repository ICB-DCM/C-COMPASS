"""Test related to SessionModel."""

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from ccompass.CCMPS import SessionModel


def test_serialization():
    """Test serialization of SessionModel."""
    # test round trip
    session = SessionModel()

    with TemporaryDirectory() as tempdir:
        fpath = Path(tempdir, "session.npy")
        session.to_numpy(fpath)
        session2 = SessionModel.from_numpy(fpath)
    assert_session_equal(session, session2)


def assert_equal(obj1, obj2):
    """Check if two objects are equal."""
    if isinstance(obj1, dict):
        for key in obj1:
            assert key in obj2
            assert_equal(obj1[key], obj2[key])
    elif isinstance(obj1, list):
        for i in range(len(obj1)):
            assert_equal(obj1[i], obj2[i])
    elif isinstance(obj1, pd.DataFrame):
        pd.testing.assert_frame_equal(obj1, obj2)
    else:
        assert obj1 == obj2


def assert_session_equal(session, session2):
    """Check if two SessionModel objects are equal."""
    for attr in session.__dict__:
        assert attr in session2.__dict__
        assert_equal(getattr(session, attr), getattr(session2, attr))
    for attr in session2.__dict__:
        assert attr in session.__dict__
