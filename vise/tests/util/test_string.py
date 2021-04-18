# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from vise.util.string import numbers_to_lowercases


def test_numbers_to_lowercases():
    assert numbers_to_lowercases("Mg2") == "Mg₂"