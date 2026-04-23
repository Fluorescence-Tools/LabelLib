import pytest
import numpy as np


def test_import():
    import labellib

    assert labellib.__version__ == "1.0.0"


def test_version():
    import labellib

    assert hasattr(labellib, "__version__")
