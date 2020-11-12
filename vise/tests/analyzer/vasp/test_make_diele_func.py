# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.vasp.make_diele_func import make_diele_func


def test_make_diele_func(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    actual = make_diele_func(v, o)
    assert actual.energies[1] == 0.0104

#
# def test_make_diele_func_2(test_data_files):
#     v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
#     o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
#     actual = make_diele_func(v, o, corrected_band_gap=2.530)
#     print(actual.diele_func_imag[1000])
#
#     actual = make_diele_func(v, o)
#     print(actual.diele_func_imag[1000])
#
# #    assert actual.energies[1] == 0.0104
