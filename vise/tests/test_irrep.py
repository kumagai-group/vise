# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.

import pytest
from vise.analyzer.irrep import Irreps, Irrep
from vise.tests.helpers.assertion import assert_msonable


@pytest.fixture()
def irreps():
    irrep = Irrep(frac_coords=[0.0, 0.0, 0.0], symbols=["GM1+"], energies=[0.1],
                  degeneracies=[1])
    return Irreps(sg_num=225, irreps={"GM": irrep})


def test_irrep_msonable(irreps):
    assert_msonable(irreps)


def test_irrep_info_set(irreps):
    assert list(irreps.irreps["GM"].irrep_info_set) == [("GM1+", 0.1, 1)]
