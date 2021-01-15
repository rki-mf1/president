import pytest
from president import __main__ as pm


def test_isavailable():
    with pytest.raises(ValueError):
        pm.is_available("asdhasdasdasd_tool")
