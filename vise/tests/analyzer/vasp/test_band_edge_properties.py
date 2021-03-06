# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen.core import Structure, Lattice
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Vasprun, Outcar, Procar
from vise.analyzer.band_edge_properties import BandEdge

from vise.analyzer.vasp.band_edge_properties import eigenvalues_from_vasprun, \
    VaspBandEdgeProperties, edge_orbital_contributions
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


def test_band_edge_properties_orbital_contributions(test_data_files):
    procar = Procar(test_data_files / "MgO_band_PROCAR")
    structure = Structure(
        Lattice.cubic(1), species=["Mg", "O"], coords=[[0]*3, [0.5]*3])
    vbm_info = BandEdge.from_dict({'energy': 3.0663, 'spin': 1, 'band_index': 3,
                                   'kpoint_index': 0,
                                   'kpoint_coords': [0.0, 0.0, 0.0]})
    cbm_info = BandEdge.from_dict({'energy': 7.726, 'spin': 1, 'band_index': 4,
                                   'kpoint_index': 0,
                                   'kpoint_coords': [0.0, 0.0, 0.0]})
    vbm, cbm = edge_orbital_contributions(procar, structure, vbm_info, cbm_info)
    assert vbm["O"]["p"] == 0.753
    assert cbm["Mg"]["s"] == 0.238


def test_eigenvalues_from_vasprun(test_data_files):
    vasprun_file = str(test_data_files / "MnO_uniform_vasprun.xml")
    vasprun = Vasprun(vasprun_file)
    assert eigenvalues_from_vasprun(vasprun)[Spin.up][1, 2] == -0.1397


