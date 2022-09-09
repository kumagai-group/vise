# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.

from vise.tests.helpers.assertion import assert_msonable


def test_irrep_msonable(irreps):
    assert_msonable(irreps)


def test_irrep_info_set(irreps):
    assert list(irreps.irreps["GM"].irrep_info_set) == [("GM1+", 0.1, 1)]
