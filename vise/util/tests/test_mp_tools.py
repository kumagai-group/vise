# -*- coding: utf-8 -*-

import shutil

from vise.util.mp_tools import make_poscars_from_mp
from vise.util.testing import ViseTest


class TestGetMpMaterials(ViseTest):
    pass


class TestMakePoscarsFromMp(ViseTest):
    def test(self):
        make_poscars_from_mp(elements=["Mg", "O"])

    # uncomment these if one wants to check the created directories.
    def tearDown(self) -> None:
        shutil.rmtree("mol_O2")
        shutil.rmtree("mp-1265_MgO")
        shutil.rmtree("mp-1094122_Mg")