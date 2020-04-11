# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.datasets.element_parameters import unoccupied_bands


def test_unoccupied_bands():
    assert unoccupied_bands["H"] == 4

