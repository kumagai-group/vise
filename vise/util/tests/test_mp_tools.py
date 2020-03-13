# -*- coding: utf-8 -*-

from glob import glob
import tempfile
import unittest
from unittest.mock import patch

from vise.util.mp_tools import get_mp_materials, make_poscars_from_mp
from vise.util.testing import ViseTest


class TestGetMpMaterials(ViseTest):

    @unittest.skipIf(not ViseTest.PMG_MAPI_KEY, ViseTest.no_mapi_key)
    def test_mprester(self) -> None:
        materials = get_mp_materials(elements=["Si"],
                                     properties=["band_gap"],
                                     e_above_hull=0.000001)
        print(materials)
        self.assertAlmostEqual(0.6244, materials[0]["band_gap"])

    @patch("pymatgen.MPRester.query")
    def test_with_mock(self, mock) -> None:
        materials = get_mp_materials(elements=["Si"],
                                     properties=["band_gap"],
                                     e_above_hull=0.000001)
        self.assertEqual(mock(), materials)
        mock.assert_called()


class TestMakePoscarsFromMp(ViseTest):
    @unittest.skipIf(not ViseTest.PMG_MAPI_KEY, ViseTest.no_mapi_key)
    def test_get_poscars(self):
        with tempfile.TemporaryDirectory() as tmp_dirname:
            make_poscars_from_mp(elements=["Mg", "O"],
                                 path=tmp_dirname)
            files = glob("/".join([tmp_dirname, "*"]))
            files = set([f.split("/")[-1] for f in files])
            expected = {"Mg_mp-1094122", "MgO_mp-1265", "mol_O2"}
            self.assertEqual(expected, files)

    @unittest.skipIf(not ViseTest.PMG_MAPI_KEY, ViseTest.no_mapi_key)
    def test_get_unary_poscars(self):
        with tempfile.TemporaryDirectory() as tmp_dirname:
            make_poscars_from_mp(elements=["Mg", "O"],
                                 path=tmp_dirname,
                                 only_unary=True)
            files = glob("/".join([tmp_dirname, "*"]))
            files = set([f.split("/")[-1] for f in files])
            expected = {"Mg_mp-1094122", "mol_O2"}
            self.assertEqual(expected, files)

    @patch("vise.util.mp_tools.get_mp_materials")
    def test_get_poscars_with_mock(self, mock):
        with tempfile.TemporaryDirectory() as tmp_dirname:
            make_poscars_from_mp(elements=["Mg", "O"],
                                 path=tmp_dirname)
        mock.assert_called()


