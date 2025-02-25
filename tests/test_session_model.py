"""Test related to SessionModel."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pydantic

from ccompass.main_gui import SessionModel


def test_serialization():
    """Test serialization of SessionModel."""
    # test round trip
    session = SessionModel()

    with TemporaryDirectory() as tempdir:
        fpath = Path(tempdir, "session.npy")
        session.to_numpy(fpath)
        session2 = SessionModel.from_numpy(fpath)
    assert_session_equal(session, session2)


def assert_equal(obj1, obj2) -> None:
    """Check if two objects are equal."""
    if isinstance(obj1, pydantic.BaseModel):
        assert isinstance(obj2, pydantic.BaseModel)
        assert_equal(obj1.model_dump(), obj2.model_dump())
        return

    if isinstance(obj1, dict):
        for key in obj1:
            assert key in obj2
            assert_equal(obj1[key], obj2[key])
    elif isinstance(obj1, list):
        for i in range(len(obj1)):
            assert_equal(obj1[i], obj2[i])
    elif isinstance(obj1, pd.DataFrame):
        assert isinstance(obj2, pd.DataFrame)
        if obj1.empty and obj2.empty:
            # if both are empty, we don't compare columns
            return
        pd.testing.assert_frame_equal(obj1, obj2, check_dtype=False)
    elif isinstance(obj1, pd.Series):
        assert isinstance(obj2, pd.Series), f"{obj1} != {obj2}"
        pd.testing.assert_series_equal(
            obj1, obj2, atol=1e-14, rtol=1e-14, check_dtype=False
        )
    elif isinstance(obj1, np.ndarray):
        np.testing.assert_almost_equal(obj1, obj2)
    elif isinstance(obj1, float) and pd.isna(obj1):
        assert pd.isna(obj2)
    else:
        assert obj1 == obj2


def assert_session_equal(session, session2):
    """Check if two SessionModel objects are equal."""
    for attr in session.__dict__:
        assert attr in session2.__dict__
        assert_equal(getattr(session, attr), getattr(session2, attr))
    for attr in session2.__dict__:
        assert attr in session.__dict__


def test_serialize_zip():
    """Test serialization of SessionModel to zip."""
    session = SessionModel()

    # round trip
    with TemporaryDirectory() as tempdir:
        fpath = Path(tempdir, "session.zip")
        session.to_zip(fpath)
        session2 = SessionModel.from_zip(fpath)

    assert_session_equal(session, session2)
