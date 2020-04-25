# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
import numpy as np
from pathlib import Path
from unittest import mock
import pytest

from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.electronic_structure.core import Spin

from vise.analyzer.band_gap import (
    BandEdge, BandEdgeProperties, VaspBandEdgeProperties,
    eigenvalues_from_vasprun)

parent_dir = Path(__file__).parent


def test_band_edge_equal():
    band_edge_1 = BandEdge(0.0, spin=Spin.up, band_index=0,
                           kpoint_coords=[0.0, 0.0, 0.0])
    band_edge_2 = BandEdge(0.0, spin=Spin.up, band_index=0,
                           kpoint_coords=[0.0, 0.0, 0.0])
    band_edge_3 = BandEdge(0.0, spin=Spin.up, band_index=0,
                           kpoint_coords=[0.1, 0.0, 0.0])
    assert band_edge_1.is_direct(band_edge_2) is True
    assert band_edge_1.is_direct(band_edge_3) is False


@pytest.fixture
def actual_kpt():
    return [[10.1, 10.2, 10.3], [10.4, 10.5, 10.6]]


expected_metal = {'energy': 0.0, 'direct': None, 'transition': None}, None, None


def test_metal_by_occupation(actual_kpt):
    eigenvalues = {Spin.up: np.array([[0.0, 0.1, 0.2], [0.2, 0.3, 0.4]])}
    band_edge = BandEdgeProperties(eigenvalues=eigenvalues,
                                   nelect=4.0,
                                   magnetization=0.0,
                                   kpoints=actual_kpt)

    assert band_edge.band_gap is None
    assert band_edge.is_direct is None
    assert band_edge.vbm is None
    assert band_edge.cbm is None


def test_metal_by_frac_nelect(actual_kpt):
    eigenvalues = {Spin.up: np.array([[0.0, 1.0, 2.0], [0.1, 1.1, 2.1]])}
    integer_criterion = 0.1
    band_edge = BandEdgeProperties(eigenvalues=eigenvalues,
                                   nelect=4.0 + integer_criterion + 1e-5,
                                   magnetization=0.0,
                                   kpoints=actual_kpt,
                                   integer_criterion=0.1)

    assert band_edge.band_gap is None
    assert band_edge.is_direct is None
    assert band_edge.vbm is None
    assert band_edge.cbm is None


def test_metal_by_frac_magnetization(actual_kpt):
    eigenvalues = {Spin.up:   np.array([[0.0, 1.0, 10.0], [0.0, 1.1, 10.0]]),
                   Spin.down: np.array([[0.0, 1.4, 10.0], [0.0, 1.5, 10.0]])}
    integer_criterion = 0.1
    band_edge = BandEdgeProperties(eigenvalues=eigenvalues,
                                   nelect=3.0,
                                   magnetization=1.0 + integer_criterion + 1e-5,
                                   kpoints=actual_kpt,
                                   integer_criterion=0.1)

    assert band_edge.band_gap is None
    assert band_edge.is_direct is None
    assert band_edge.vbm is None
    assert band_edge.cbm is None


def test_nonmag_insulator(actual_kpt):
    # k-point indices run fast.
    eigenvalues = {Spin.up: np.array([[0.0, 1.0, 2.0], [0.1, 1.1, 2.1]])}
    integer_criterion = 0.1
    band_edge = BandEdgeProperties(eigenvalues=eigenvalues,
                                   nelect=4.0 + integer_criterion - 1e-5,
                                   magnetization=0.0,
                                   kpoints=actual_kpt,
                                   integer_criterion=0.1)

    assert pytest.approx(band_edge.band_gap) == 0.90

    assert band_edge.is_direct is False

    assert pytest.approx(band_edge.vbm.energy) == 1.1
    assert pytest.approx(band_edge.cbm.energy) == 2.0

    assert band_edge.vbm.spin == Spin.up
    assert band_edge.cbm.spin == Spin.up

    assert band_edge.vbm.band_index == 1
    assert band_edge.cbm.band_index == 2

    assert band_edge.vbm.kpoint_coords == [10.4, 10.5, 10.6]
    assert band_edge.cbm.kpoint_coords == [10.1, 10.2, 10.3]


def test_mag_insulator(actual_kpt):
    eigenvalues = {Spin.up:   np.array([[0.0, 1.0, 10.0], [0.0, 1.1, 10.0]]),
                   Spin.down: np.array([[0.0, 1.4, 10.0], [0.0, 1.5, 10.0]])}
    band_edge = BandEdgeProperties(eigenvalues=eigenvalues,
                                   nelect=3.0,
                                   magnetization=1.0,
                                   kpoints=actual_kpt)

    assert pytest.approx(band_edge.band_gap) == 0.3

    assert band_edge.is_direct is False

    assert pytest.approx(band_edge.vbm.energy) == 1.1
    assert pytest.approx(band_edge.cbm.energy) == 1.4

    assert band_edge.vbm.spin == Spin.up
    assert band_edge.cbm.spin == Spin.down

    assert band_edge.vbm.band_index == 1
    assert band_edge.cbm.band_index == 1

    assert band_edge.vbm.kpoint_coords == [10.4, 10.5, 10.6]
    assert band_edge.cbm.kpoint_coords == [10.1, 10.2, 10.3]


"""
TODO:
+ Read vasprun.xml file and Outcar file
+ Construct eigenvalues, nelect, magnetization, kpoints
"""


def test_eigenvalues_from_vasprun():
    vasprun = Vasprun(parent_dir / "MnO_uniform_vasprun.xml")
    assert eigenvalues_from_vasprun(vasprun)[Spin.up][1, 2] == -0.1397


def test_band_edge_properties_from_vasp():
    vasprun = Vasprun(parent_dir / "MnO_uniform_vasprun.xml")
    outcar = Outcar(parent_dir / "MnO_uniform_OUTCAR")
    band_edge = VaspBandEdgeProperties(vasprun, outcar)
    assert pytest.approx(band_edge.band_gap) == 0.4702
