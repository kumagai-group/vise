# -*- coding: utf-8 -*-

import tempfile

from vise.input_set.prior_info import PriorInfo
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class PriorInfoTest(PydefectTest):

    def setUp(self):
        """ """
        energy_per_atom = -0.5
        band_gap = 1.2
        total_magnetization = 0.000001
        data_source = "Materials Project"
        is_molecule = False
        mag_threshold = 0.001
        band_gap_threshold = 0.1

        self.nonmagnetic_insulator = \
            PriorInfo(energy_per_atom, band_gap, total_magnetization,
                      data_source, is_molecule, mag_threshold,
                      band_gap_threshold)

    def test_dict(self):
        # round trip test of dict
        d = self.nonmagnetic_insulator.as_dict()
        prior_info_from_dict = PriorInfo.from_dict(d)
        self.assertTrue(d == prior_info_from_dict.as_dict())

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.nonmagnetic_insulator.dump_json(tmp_file.name)
        prior_info_from_json = PriorInfo.load_json(tmp_file.name)
#        self.nonmagnetic_insulator.dump_json("aaa")
#        prior_info_from_json = PriorInfo.load_json("aaa")
        self.assertEqual(prior_info_from_json.as_dict(),
                         self.nonmagnetic_insulator.as_dict())

    def test_msonable(self):
        self.assertMSONable(self.nonmagnetic_insulator)
