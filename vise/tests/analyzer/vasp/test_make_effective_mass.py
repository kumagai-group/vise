# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen.electronic_structure.boltztrap import BoltztrapError

from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties

try:
    from pymatgen.electronic_structure.boltztrap2 import VasprunBSLoader, \
        BztInterpolator, BztTransportProperties
    BOLTZTRAP2_NOT_PRESENT = False
except BoltztrapError:
    BOLTZTRAP2_NOT_PRESENT = True


@pytest.mark.skipif(BOLTZTRAP2_NOT_PRESENT, reason="boltztrap2 is not installed.")
def test_make_effective_mass(test_data_files):
#def test_make_effective_mass(mocker, test_data_files):
    from vise.analyzer.effective_mass import EffectiveMass
    from vise.analyzer.vasp.make_effective_mass import make_effective_mass
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    # mock_vl = mocker.patch("vise.analyzer.vasp.make_effective_mass.VasprunBSLoader")
    # mock_bi = mocker.patch("vise.analyzer.vasp.make_effective_mass.BztInterpolator")
    # mock_btp = mocker.patch("vise.analyzer.vasp.make_effective_mass.BztTransportProperties")

    p = [[[ 2.27995989e+00, -7.16015256e-17, -1.92884504e-16],
          [-7.16015256e-17,  2.27995989e+00, -1.17150055e-16],
          [-1.92884504e-16, -1.17150055e-16,  2.27995989e+00]]]
    n = [[[ 5.49195397e-01, -1.09903151e-17, -3.07622696e-17],
          [-1.09903151e-17,  5.49195397e-01, -5.50277544e-17],
          [-3.07622696e-17, -5.50277544e-17,  5.49195397e-01]]]
    # mock_btp.return_value.Effective_mass_doping = {"p": np.array([p]),
    #                                                "n": np.array([n])}
    # assume the situation where the efermi locates lower than CBM.
    v.efermi = 0.55
    actual = make_effective_mass(vasprun=v, temp=300, concentrations=[1e18],
                                 vbm=0.5614, cbm=3.0904)

    # mock_bi.assert_called_with(mock_vl.return_value, energy_range=3.0)
    # mock_btp.assert_called_with(mock_bi.return_value, temp_r=np.array([300]))
    # mock_btp.return_value.compute_properties_doping.assert_called_with([1e18])

    expected = EffectiveMass(p=p, n=n, temperature=300, concentrations=[1e+18])
    np.testing.assert_array_almost_equal(actual.p, expected.p)
    np.testing.assert_array_almost_equal(actual.n, expected.n)
    np.testing.assert_array_almost_equal(actual.concentrations,
                                         expected.concentrations)

