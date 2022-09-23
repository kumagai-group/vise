# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from math import pi, sqrt
from pathlib import Path

import numpy as np
import pytest
from monty.serialization import loadfn
from vise.analyzer.dielectric_function import DieleFuncData, \
    eV_to_inv_cm, min_e_w_target_coeff
from vise.tests.helpers.assertion import assert_msonable, \
    assert_dataclass_almost_equal

try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


@pytest.fixture
def diele_func_data():
    return DieleFuncData(energies=[0.0, 10.0],
                         directions=["xx", "yy", "ave"],
                         diele_func_real=[[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]],
                         diele_func_imag=[[1.1, 1.2], [1.3, 1.4], [1.5, 1.6]],
                         band_gap=1.0)


def test_msonable(diele_func_data):
    assert_msonable(diele_func_data)


def test_json_file_mixin(diele_func_data, tmpdir):
    tmpdir.chdir()
    diele_func_data.to_json_file()
    actual = loadfn("diele_func_data.json")
    assert actual.diele_func_real == diele_func_data.diele_func_real


def test_absorption_coeff(diele_func_data):
    actual = diele_func_data.absorption_coeff[2][1]
    real, imag = 0.6, 1.6
    expected = (2 * sqrt(2) * pi * sqrt(sqrt(real ** 2 + imag ** 2) - real)
                * 10.0 * eV_to_inv_cm)
    assert actual == expected


def test_refractive_idx(diele_func_data):
    e_real, e_imag = 0.6, 1.6  # values at xx

    actual = diele_func_data.refractive_idx_real[2][1]
    expected = sqrt(e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)
    assert actual == expected

    actual = diele_func_data.refractive_idx_imag[2][1]
    expected = sqrt(-e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)
    assert actual == expected


def test_reflectivity(diele_func_data):
    e_real, e_imag = 0.6, 1.6  # values at xx
    n = sqrt(e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)
    k = sqrt(-e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)
    expected = ((n - 1)**2 + k**2) / ((n + 1)**2 + k**2)
    actual = diele_func_data.reflectivity[2][1]
    assert actual == expected


def test_to_csv_file(tmpdir):
    diele = DieleFuncData(energies=[0.0, 3.0],
                          diele_func_real=[[0.1, 0.1], [0.2, 0.2], [0.3, 0.3]],
                          diele_func_imag=[[1.1, 1.1], [1.2, 1.2], [1.3, 1.3]],
                          directions=["xx", "yy", "ave"],
                          band_gap=1.0)
    tmpdir.chdir()
    diele.to_csv_file()
    actual = Path("diele_func_data.csv").read_text()
    expected = """energies(eV),real_xx,real_yy,real_ave,imag_xx,imag_yy,imag_ave,band_gap
0.0,0.1,0.2,0.3,1.1,1.2,1.3,1.0
3.0,0.1,0.2,0.3,1.1,1.2,1.3,
"""
    assert actual == expected

    actual = DieleFuncData.from_csv_file("diele_func_data.csv")
    assert_dataclass_almost_equal(actual, diele, digit=3)


def test_target_coeff_e_from_band_gap():
    actual = min_e_w_target_coeff([0.0, 0.1, 0.2], [1, 2, 3], 1.5)
    assert actual == 0.1


# def test_actual_diele_func_data_with_kk_trans():
#     diele_func_data = DieleFuncData(
#         energies=[0.0, 1.0],
#         diele_func_real=[[]])
#     diele_func = make_diele_func(v, o, use_vasp_real=False)
#     print(diele_func.diele_func_real)


