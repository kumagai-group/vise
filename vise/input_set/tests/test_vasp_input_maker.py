# -*- coding: utf-8 -*-
import filecmp
import unittest

from vise.input_set.input_set import ModIncar, make_potcar, \
    make_hpkot_primitive_poscar, make_band_kpoints, make_kpoints, MakeIncar

from vise.util.structure_handler import NotAppropriatePrimitiveError


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ModIncarTest(unittest.TestCase):

    def test_ModIncar(self):
        i = ModIncar.from_file("INCAR-ModIncar_before")
        i.pretty_write_file("INCAR-ModIncar_actual")
        self.assertTrue(
            filecmp.cmp("INCAR-ModIncar_expected", "INCAR-ModIncar_actual"))


class MakePotcarTest(unittest.TestCase):

    def test(self):
        elements = ("Mg", "O")
        make_potcar(elements=elements)
        self.assertTrue(
            filecmp.cmp("POTCAR-MgO", "POTCAR"))


class MakeHpkotPrimitiveTest(unittest.TestCase):

    def test(self):
        make_hpkot_primitive_poscar(poscar="PPOSCAR-YMnO3", pposcar="PPOSCAR")
        self.assertTrue(
            filecmp.cmp("PPOSCAR-YMnO3-hpkot", "PPOSCAR"))





class MakeIncarTest(unittest.TestCase):

    # TODO: write better tests
    def test_structure_opt(self):
        MakeIncar(task="structure_opt", functional="hse",
                  defect_in="defect.in", hfscreen=0.2, aexx=0.3,
                  is_magnetization=True, ldau=True,
                  my_incar_setting="my_INCAR_setting.yaml")

    def test_band(self):
        MakeIncar(task="band", functional="pbe", poscar="PPOSCAR-YMnO3",
                  ldau=True, my_incar_setting="my_INCAR_setting.yaml")

    def test_dos(self):
        MakeIncar(task="dos", functional="pbe", poscar="PPOSCAR-YMnO3",
                  ldau=True, is_magnetization=True)

    def test_dos(self):
        MakeIncar(task="dielectric_function", functional="pbe",
                  poscar="PPOSCAR-YMnO3", ldau=True, is_magnetization=True)


if __name__ == "__main__":
    unittest.main()
