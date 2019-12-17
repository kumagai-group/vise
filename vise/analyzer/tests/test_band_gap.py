# -*- coding: utf-8 -*-
from pymatgen.io.vasp.outputs import Vasprun, Outcar

from vise.analyzer.band_gap import band_gap_properties
from vise.util.testing import ViseTest


class BandGapPropertiesTest(ViseTest):
    def setUp(self) -> None:
        pass

    def test_mgo(self):
        vasprun = Vasprun("MgO_band_vasprun.xml")
        outcar = Outcar("MgO_band_OUTCAR")
        expected = (4.6597,
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
        vasprun = Vasprun("BaTa6O16_band_vasprun.xml")
        outcar = Outcar("BaTa6O16_band_OUTCAR")
        expected = (2.6454000000000004,
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
        vasprun = Vasprun("MnO_uniform_vasprun.xml")
        outcar = Outcar("MnO_uniform_OUTCAR")
        expected = (0.47020000000000017,
                    {'energy': 4.6665999999999999,
                     'spin': 1,
                     'band_index': 8,
                     'kpoints': [0.42857143, 0.0, 0.0]},
                    {'energy': 5.1368,
                     'spin': 1,
                     'band_index': 9,
                     'kpoints': [0.0, 0.0, 0.0]})
        self.assertEqual(expected, band_gap_properties(vasprun, outcar))

