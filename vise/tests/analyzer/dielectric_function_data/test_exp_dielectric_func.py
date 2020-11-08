# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from vise.analyzer.dielectric_function_data.exp_dielectric_func import \
    ExpDieleFunc
import numpy as np

import pytest


def test_gaas():
    diele_gaas = ExpDieleFunc("GaAs")
    assert diele_gaas.band_gap == 1.41
    assert diele_gaas.energies[10] == 1.41
    np.testing.assert_almost_equal(diele_gaas.dielectric_real[10], 13.16)
    np.testing.assert_almost_equal(diele_gaas.dielectric_imag[10], 0.013)
    np.testing.assert_almost_equal(diele_gaas.absorption_coeff[10], 2.58e+02)


def test_others():
    assert ExpDieleFunc("Si").band_gap == 1.17
    assert ExpDieleFunc("Si").energies[21] == 0.85


def test_file_not_found_error():
    with pytest.raises(FileNotFoundError):
        ExpDieleFunc("Au")
