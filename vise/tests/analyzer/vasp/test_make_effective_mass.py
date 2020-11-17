# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pymatgen.io.vasp import Vasprun
from vise.analyzer.effective_mass import EffectiveMass
from vise.analyzer.vasp.make_effective_mass import make_effective_mass


def test_effective_mass_json_file_mixin(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    actual = make_effective_mass(vasprun=v, temp=300, concentrations=[1e18],
                                 band_gap=1.0)
    expected = EffectiveMass(
        p=[[[ 2.27995989e+00, -7.16015256e-17, -1.92884504e-16],
            [-7.16015256e-17,  2.27995989e+00, -1.17150055e-16],
            [-1.92884504e-16, -1.17150055e-16,  2.27995989e+00]]],
        n=[[[ 5.49195397e-01, -1.09903151e-17, -3.07622696e-17],
                       [-1.09903151e-17,  5.49195397e-01, -5.50277544e-17],
                       [-3.07622696e-17, -5.50277544e-17,  5.49195397e-01]]],
        concentrations=[1e+18])
    np.testing.assert_array_almost_equal(actual.p, expected.p)
    np.testing.assert_array_almost_equal(actual.n, expected.n)
    np.testing.assert_array_almost_equal(actual.concentrations,
                                         expected.concentrations)


