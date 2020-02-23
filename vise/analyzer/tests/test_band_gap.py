# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from pathlib import Path
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.electronic_structure.core import Spin

from vise.analyzer.band_gap import band_gap_properties, edge_info
from vise.util.testing import ViseTest


parent_dir = Path(__file__).parent


class BandGapPropertiesTest(ViseTest):
    def test_mgo(self):
        vasprun = Vasprun(parent_dir / "MgO_band_vasprun.xml")
        outcar = Outcar(parent_dir / "MgO_band_OUTCAR")
        expected = ({'energy': 4.6597,
                     'direct': True},
                    {'energy': 3.0663,
                     'spin': None,
                     'band_index': 3,
                     'kpoints': [0.0, 0.0, 0.0]},
                    {'energy': 7.726,
                     'spin': None,
                     'band_index': 4,
                     'kpoints': [0.0, 0.0, 0.0]})
        self.assertEqual(expected, band_gap_properties(vasprun, outcar))

    def test_bata6o16(self):
        vasprun = Vasprun(parent_dir / "BaTa6O16_band_vasprun.xml")
        outcar = Outcar(parent_dir / "BaTa6O16_band_OUTCAR")
        expected = ({'energy': 2.6454000000000004,
                     'direct': False},
                    {'energy': 1.8945000000000001,
                     'spin': None,
                     'band_index': 67,
                     'kpoints': [-0.5, 0.5, 0.46774194]},
                    {'energy': 4.5399000000000003,
                     'spin': None,
                     'band_index': 68,
                     'kpoints': [0.0, 0.0, 0.5]})
        self.assertEqual(expected, band_gap_properties(vasprun, outcar))

    def test_mno(self):
        vasprun = Vasprun(parent_dir / "MnO_uniform_vasprun.xml")
        outcar = Outcar(parent_dir / "MnO_uniform_OUTCAR")
        expected = ({'energy': 0.47020000000000017,
                     'direct': False},
                    {'energy': 4.6665999999999999,
                     'spin': 1,
                     'band_index': 8,
                     'kpoints': [0.42857143, 0.0, 0.0]},
                    {'energy': 5.1368,
                     'spin': 1,
                     'band_index': 9,
                     'kpoints': [0.0, 0.0, 0.0]})
        self.assertEqual(expected, band_gap_properties(vasprun, outcar))


class EdgeInfoTest(ViseTest):
    def setUp(self) -> None:
        values = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        # [[[ 1.  2.]   1st k: 1st band (hob)
        #   [ 3.  4.]]  1st k: 2nd band (lub)
        #
        #  [[ 5.  6.]   2nd k: 1st band (hob)
        #   [ 7.  8.]]] 2nd k: 2st band (lub)
        self.eigenvalues = {Spin.up: values}

    def test(self):
        actual = edge_info(self.eigenvalues,
                           hob_index=0,
                           spin=Spin.up)
        expected = 5.0, 1, 3.0, 0
        self.assertEqual(expected, actual)



