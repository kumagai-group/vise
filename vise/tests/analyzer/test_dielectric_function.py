# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from math import pi, sqrt
from pathlib import Path

import numpy as np
import pytest
from monty.serialization import loadfn
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import DieleFuncData, \
    eV_to_inv_cm
from vise.analyzer.vasp.make_diele_func import make_diele_func
from vise.tests.helpers.assertion import assert_msonable

try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


@pytest.fixture
def diele_func_data():
    array_real = [1, 2, 3, 0, 0, 0]
    array_imag = [4, 5, 6, 0, 0, 0]
    return DieleFuncData(energies=list(np.linspace(0.0, 10.0, num=11)),
                         diele_func_real=[array_real]*11,
                         diele_func_imag=[array_imag]*11,
                         band_gap=1.0)


def test_msonable(diele_func_data):
    assert_msonable(diele_func_data)


def test_json_file_mixin(diele_func_data, tmpdir):
    tmpdir.chdir()
    diele_func_data.to_json_file()
    actual = loadfn("diele_func_data.json")
    assert actual.diele_func_real == diele_func_data.diele_func_real


def test_ave_absorption_coeff(diele_func_data):
    actual = diele_func_data.ave_absorption_coeff
    expected = [2 * sqrt(2) * pi * sqrt(sqrt(2 ** 2 + 5 ** 2) - 2)
                * i * eV_to_inv_cm for i in range(0, 11)]
    assert actual == expected


def test_to_csv_file(tmpdir):
    diele = DieleFuncData(energies=[0.0],
                          diele_func_real=[[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]],
                          diele_func_imag=[[1.1, 1.2, 1.3, 1.4, 1.5, 1.6]],
                          band_gap=1.0)
    tmpdir.chdir()
    diele.to_csv_file()
    actual = Path("diele_func_data.csv").read_text()
    expected = """energies(eV),real_xx,real_yy,real_zz,real_xy,real_yz,real_xz,imag_xx,imag_yy,imag_zz,imag_xy,imag_yz,imag_xz
0.0,0.1,0.2,0.3,0.4,0.5,0.6,1.1,1.2,1.3,1.4,1.5,1.6
"""
    assert actual == expected


@pytest.fixture
def actual_diele_func_data(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    diele_func = make_diele_func(v, o)
    print(diele_func)
    return diele_func


def test_target_coeff_e_from_band_gap(actual_diele_func_data):
    actual = actual_diele_func_data.target_coeff_min_e()
    np.testing.assert_almost_equal(actual, 2.9304)


# def test_actual_diele_func_data_with_kk_trans():
#     diele_func_data = DieleFuncData(
#         energies=[0.0, 1.0],
#         diele_func_real=[[]])
#     diele_func = make_diele_func(v, o, use_vasp_real=False)
#     print(diele_func.diele_func_real)


