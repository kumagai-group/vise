# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.util.bravais_lattice import BravaisLattice


@pytest.mark.parametrize("sg_num, bravais_lattice",
                         [(1, BravaisLattice.aP),
                          (3, BravaisLattice.mP),
                          (5, BravaisLattice.mC),
                          (16, BravaisLattice.oP),
                          (20, BravaisLattice.oC),
                          (22, BravaisLattice.oF),
                          (23, BravaisLattice.oI),
                          (38, BravaisLattice.oA),
                          (75, BravaisLattice.tP),
                          (79, BravaisLattice.tI),
                          (143, BravaisLattice.hP),
                          (146, BravaisLattice.hR),
                          (195, BravaisLattice.cP),
                          (196, BravaisLattice.cF),
                          (197, BravaisLattice.cI),
                          ])
def test_bravais_lattice_from_sg_num(sg_num, bravais_lattice):
    assert BravaisLattice.from_sg_num(sg_num) == bravais_lattice


def test_kpt_centering():
    assert BravaisLattice.oF.kpt_centering == [0.0, 0.0, 0.0]
    assert BravaisLattice.hP.kpt_centering == [0.0, 0.0, 0.5]
    assert BravaisLattice.cF.kpt_centering == [0.0, 0.0, 0.0]


def test_need_same_num_kpt():
    assert BravaisLattice.oI.need_same_num_kpt is True
    assert BravaisLattice.cP.need_same_num_kpt is False
