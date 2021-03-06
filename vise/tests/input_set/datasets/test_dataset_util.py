# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from math import ceil

from pymatgen.core import Composition
from pymatgen.io.vasp import Potcar

from vise.input_set.datasets.dataset_util import (
    LDAU, unoccupied_bands, num_bands, npar_kpar)


def test_ldau_3d_transition_metal():
    actual = LDAU(symbol_list=["Mn", "Cu"], )
    assert actual.ldauu == [3, 5]
    assert actual.ldaul == [2, 2]
    assert actual.lmaxmix == 4


def test_ldau_rare_earth():
    actual = LDAU(symbol_list=["Mn", "Ho"], )
    assert actual.ldauu == [3, 5]
    assert actual.ldaul == [2, 3]
    assert actual.lmaxmix == 6


def test_ldau_3d_override():
    actual = LDAU(symbol_list=["Mn", "Cu"])
    assert actual.ldauu == [3, 5]
    assert actual.ldaul == [2, 2]
    assert actual.lmaxmix == 4


def test_unoccupied_bands():
    assert unoccupied_bands["H"] == 4


def test_nbands():
    composition = Composition("KH2")
    potcar = Potcar(symbols=["Mn", "O"])
    nelect_mn_pv = 7
    nelect_o = 6
    unoccupied_band_mn = 9
    unoccupied_band_o = 4
    expected = (ceil((nelect_mn_pv + nelect_o * 2) / 2)
                + unoccupied_band_mn + unoccupied_band_o * 2)
    assert num_bands(composition, potcar) == expected


def test_npar_kpar():
    assert npar_kpar(num_kpoints=10, num_nodes=5) == (2, 1)
    assert npar_kpar(num_kpoints=100, num_nodes=4) == (16, 1)

