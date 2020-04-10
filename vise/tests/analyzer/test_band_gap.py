# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
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
                      {'energy': 1.1,
                       'spin': None,
                       'band_index': 1,
                       'kpoints': [10.4, 10.5, 10.6]},
                      {'energy': 2.0,
                       'spin': None,
                       'band_index': 2,
                       'kpoints': [10.1, 10.2, 10.3]})

expected_metal = {'energy': 0.0, 'direct': None, 'transition': None}, None, None


def test_nonmag_insulator(actual_kpt):
    # k-point indices run fast.
    eigenvalues = {Spin.up: np.array([[0.0, 1.0, 2.0], [0.1, 1.1, 2.1]])}

    assert expected_insulator == band_gap_properties(eigenvalues=eigenvalues,
                                                     nelect=2.0,
                                                     magnetization=0.0,
                                                     kpoints=actual_kpt)


@pytest.mark.parametrize("cbm_energy,exp_gap,exp_cbm",
                         [(2.0001, 0.9001, 2.0001),
                          (2.00001, 0.9, 2.0)])
def test_round_ndigits(actual_kpt, cbm_energy, exp_gap, exp_cbm):
    eigenvalues = {Spin.up: np.array([[0.0, 1.0, cbm_energy], [0.1, 1.1, 2.1]])}

    actual = band_gap_properties(eigenvalues=eigenvalues,
                                 nelect=2.0,
                                 magnetization=0.0,
                                 kpoints=actual_kpt)

    assert exp_gap == actual[0]["energy"]
    assert exp_cbm == actual[2]["energy"]


expected_mag_insulator = ({'energy': 0.3, 'direct': False},
                          {'energy': 1.1, 'spin': 1,
                           'band_index': 1, 'kpoints': [10.4, 10.5, 10.6]},
                          {'energy': 1.4, 'spin': -1,
                           'band_index': 1, 'kpoints': [10.1, 10.2, 10.3]})


def test_mag_insulator(actual_kpt):
    eigenvalues = {Spin.up:   np.array([[0.0, 1.0, 10.0], [0.0, 1.1, 10.0]]),
                   Spin.down: np.array([[0.0, 1.4, 10.0], [0.0, 1.5, 10.0]])}
    assert expected_mag_insulator == band_gap_properties(eigenvalues=eigenvalues,
                                                         nelect=3.0,
                                                         magnetization=1.0,
                                                         kpoints=actual_kpt)

