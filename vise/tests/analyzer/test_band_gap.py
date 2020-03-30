# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from pathlib import Path
from unittest import mock
import pytest

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin

from vise.analyzer.band_gap import band_gap_properties


parent_dir = Path(__file__).parent


@pytest.fixture
def actual_kpt():
    return [[10.1, 10.2, 10.3], [10.4, 10.5, 10.6]]


expected_insulator = ({'energy': 0.9, 'direct': False},
                      {'energy': 1.1, 'spin': None,
                       'band_index': 1, 'kpoints': [10.4, 10.5, 10.6]},
                      {'energy': 2.0, 'spin': None,
                       'band_index': 2, 'kpoints': [10.1, 10.2, 10.3]})

expected_metal = {'energy': 0.0, 'direct': None, 'transition': None}, None, None


@pytest.mark.parametrize("frac_val,expected", [(1.0, expected_insulator),
                                               (0.901, expected_insulator),
                                               (0.899, expected_metal)])
def test_nonmag_insulator(actual_kpt, frac_val, expected):
    vasprun = mock.Mock(Vasprun)
    # k-point indices run fast.
    vasprun.eigenvalues = \
        {Spin.up: np.array([[[0.0, 1.0], [1.0, 1.0], [2.0, 0.0]],
                            [[0.1, 1.0], [1.1, frac_val], [2.1, 0.0]]])}
    vasprun.actual_kpoints = actual_kpt

    assert expected == band_gap_properties(vasprun)


expected_mag_insulator = ({'energy': 0.3, 'direct': False},
                          {'energy': 1.1, 'spin': 1,
                           'band_index': 1, 'kpoints': [10.4, 10.5, 10.6]},
                          {'energy': 1.4, 'spin': -1,
                           'band_index': 1, 'kpoints': [10.1, 10.2, 10.3]})


@pytest.mark.parametrize("frac_val,expected", [(1.0, expected_mag_insulator),
                                               (0.901, expected_mag_insulator),
                                               (0.899, expected_metal)])
def test_mag_insulator(actual_kpt, frac_val, expected):
    vasprun = mock.Mock(Vasprun)
    vasprun.eigenvalues = \
        {Spin.up: np.array([[[0.0, 1.0], [1.0, 1.0], [10.0, 0.0]],
                            [[0.0, 1.0], [1.1, frac_val], [10.0, 0.0]]]),
         Spin.down: np.array([[[0.0, 1.0], [1.4, 1 - frac_val], [10.0, 0.0]],
                             [[0.0, 1.0], [1.5, 0.0], [10.0, 0.0]]])}
    vasprun.actual_kpoints = actual_kpt

    assert expected == band_gap_properties(vasprun)

