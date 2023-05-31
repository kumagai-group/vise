# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.

import pytest
from vise.util.structure_handler import sanitize_matrix


def test_sanitize_matrix_9_input_values():
    actual = sanitize_matrix(list(range(9)))
    expected = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    assert actual == expected


def test_sanitize_matrix_3_input_values():
    actual = sanitize_matrix(list(range(1, 4)))
    expected = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
    assert actual == expected


def test_sanitize_matrix_1_input_value():
    actual = sanitize_matrix([4])
    expected = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
    assert actual == expected


