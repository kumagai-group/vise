# -*- coding: utf-8 -*-

import pytest
from pathlib import Path


@pytest.fixture(scope="session", autouse=True)
def test_data_files():
    return Path(__file__).parent / "test_data_files"

