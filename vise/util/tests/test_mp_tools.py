# -*- coding: utf-8 -*-

import shutil
from unittest.mock import patch

from pymatgen import Element

from vise.util.mp_tools import get_mp_materials, make_poscars_from_mp
from vise.util.testing import ViseTest


class TestGetMpMaterials(ViseTest):

    def test(self) -> None:
        materials = get_mp_materials(elements=["Si"],
                                     properties=["band_gap"],
                                     e_above_hull=0.000001)
        print(materials)
        self.assertAlmostEqual(0.6244, materials[0]["band_gap"])


class TestMakePoscarsFromMp(ViseTest):
    def test(self):
        make_poscars_from_mp(elements=["Mg", "O"])

    # uncomment these if one wants to check the created directories.
    def tearDown(self) -> None:
        shutil.rmtree("mol_O2")
        shutil.rmtree("MgO_mp-1265")
        shutil.rmtree("Mg_mp-1094122")
