# -*- coding: utf-8 -*-

import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_files():
    return Path(__file__).parent / "test_data_files"

