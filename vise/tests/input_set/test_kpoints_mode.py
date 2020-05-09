# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.kpoints_mode import KpointsMode


def test_kpoints_mode_from_string():
    assert KpointsMode.from_string("band") == KpointsMode.band


def test_kpoints_mode_from_string_raise_error():
    with pytest.raises(AttributeError):
        KpointsMode.from_string("fail_string")


def test_kpoints_mode_name_list():
    expected = "band, primitive, uniform"
    assert KpointsMode.names_string() == expected
