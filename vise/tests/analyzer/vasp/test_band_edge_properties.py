# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen import Spin
from pymatgen.io.vasp import Vasprun, Outcar

from vise.analyzer.vasp.band_edge_properties import eigenvalues_from_vasprun, \
    VaspBandEdgeProperties
from vise.tests.analyzer.test_band_edge_properties import parent_dir


def test_band_edge_properties_from_vasp(test_data_files):
    vasprun_file = str(test_data_files / "MnO_uniform_vasprun.xml")
    vasprun = Vasprun(vasprun_file)
    outcar_file = str(test_data_files / "MnO_uniform_OUTCAR")
    outcar = Outcar(outcar_file)
    band_edge = VaspBandEdgeProperties(vasprun, outcar)

    assert pytest.approx(band_edge.band_gap) == 0.4702


def test_band_edge_properties_from_vasp_non_mag(test_data_files):
    vasprun_file = str(test_data_files / "MgO_band_vasprun.xml")
    vasprun = Vasprun(vasprun_file)
    outcar_file = str(test_data_files / "MgO_band_OUTCAR")
    outcar = Outcar(outcar_file)
    band_edge = VaspBandEdgeProperties(vasprun, outcar)

    assert pytest.approx(band_edge.band_gap) == 4.6597


def test_eigenvalues_from_vasprun(test_data_files):
    vasprun_file = str(test_data_files / "MnO_uniform_vasprun.xml")
    vasprun = Vasprun(vasprun_file)
    assert eigenvalues_from_vasprun(vasprun)[Spin.up][1, 2] == -0.1397


