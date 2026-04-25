import pytest
import numpy as np
from packaging.version import Version


def test_import():
    import labellib

    assert isinstance(labellib.__version__, str)
    assert Version(labellib.__version__).release[:2] == (1, 0)


def test_version():
    import labellib

    assert hasattr(labellib, "__version__")
