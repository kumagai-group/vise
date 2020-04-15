# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pymatgen import Composition
from pymatgen.io.vasp import Potcar

from vise.input_set.datasets.dataset_util import (
    unoccupied_bands, nbands, npar_kpar)


def test_unoccupied_bands():
    assert unoccupied_bands["H"] == 4


def test_nbands():
    composition = Composition("KH2")
    potcar = Potcar(symbols=["K_pv", "H"], functional="PBE")
    # Nelect of K_pv is 7, so (7 + 1 * 2) / 2 + (9 + 4 * 2) = 22
    assert nbands(composition, potcar) == 22


def test_npar_kpar():
    assert npar_kpar(num_kpoints=10, num_nodes=5) == (2, 1)
    assert npar_kpar(num_kpoints=100, num_nodes=4) == (16, 1)

