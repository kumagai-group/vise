# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from vise.input_set.fft_grids import vasp_grid


def test_vasp_grids():
    actual = vasp_grid(encut=250, lattice_length=10.0, symprec="Normal")
    assert actual == 39


