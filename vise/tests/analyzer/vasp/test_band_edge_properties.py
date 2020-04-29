# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen import Spin
from pymatgen.io.vasp import Vasprun, Outcar

from vise.analyzer.vasp.band_edge_properties import eigenvalues_from_vasprun, \
    VaspBandEdgeProperties
from vise.tests.analyzer.test_band_edge_properties import parent_dir


def test_eigenvalues_from_vasprun():
    vasprun = Vasprun(parent_dir / "MnO_uniform_vasprun.xml")
    assert eigenvalues_from_vasprun(vasprun)[Spin.up][1, 2] == -0.1397


def test_band_edge_properties_from_vasp():
    vasprun = Vasprun(parent_dir / "MnO_uniform_vasprun.xml")
    outcar = Outcar(parent_dir / "MnO_uniform_OUTCAR")
    band_edge = VaspBandEdgeProperties(vasprun, outcar)
    assert pytest.approx(band_edge.band_gap) == 0.4702