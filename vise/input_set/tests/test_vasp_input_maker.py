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


class MakeBandKpointsTest(unittest.TestCase):

    def test_make_band_kpoints_mgo(self):
        make_band_kpoints(ibzkpt="IBZKPT", poscar="PPOSCAR-MgO",
                          num_split_kpoints=2)
        self.assertTrue(filecmp.cmp("KPOINTS-make_band-MgO-0", "KPOINTS-0"))
        self.assertTrue(filecmp.cmp("KPOINTS-make_band-MgO-1", "KPOINTS-1"))

#    def test_make_band_kpoints_ymo3_fail(self):
#        with self.assertRaises(NotAppropriatePrimitiveError):
#            make_band_kpoints(ibzkpt="IBZKPT", poscar="BPOSCAR-MgO")

    def test_make_band_kpoints_ymo3_succeed(self):
        make_band_kpoints(ibzkpt="IBZKPT", poscar="PPOSCAR-YMnO3-hpkot")
        self.assertTrue(filecmp.cmp("KPOINTS-make_band-YMnO3", "KPOINTS"))


class MakeKpointsTest(unittest.TestCase):

    def test_structure_opt1(self):
        make_kpoints(task="structure_opt", poscar="PPOSCAR-MgO")
        self.assertTrue(filecmp.cmp("KPOINTS-so_expected_MgO", "KPOINTS"))

    def test_structure_opt2(self):
        make_kpoints(task="structure_opt", poscar="PPOSCAR-YMnO3")
        self.assertTrue(filecmp.cmp("KPOINTS-so_expected_YMnO3", "KPOINTS"))

    def test_band(self):
        make_kpoints(task="band", poscar="PPOSCAR-YMnO3-hpkot")
        self.assertTrue(filecmp.cmp("KPOINTS-make_band-YMnO3", "KPOINTS"))

    def test_dos(self):
        make_kpoints(task="dos", poscar="PPOSCAR-YMnO3")
        self.assertTrue(filecmp.cmp("KPOINTS-dos_expected_YMnO3", "KPOINTS"))

    def test_dos2(self):
        make_kpoints(task="dos", poscar="PPOSCAR-ZnSnP2")
        self.assertTrue(filecmp.cmp("KPOINTS-dos_expected_ZnSnP2", "KPOINTS"))

    def test_dielectric(self):
        make_kpoints(task="dielectric", poscar="PPOSCAR-YMnO3-hpkot")
        self.assertTrue(filecmp.cmp("KPOINTS-dielectric-YMnO3", "KPOINTS"))

    def test_dielectric_function(self):
        make_kpoints(task="dielectric_function", poscar="PPOSCAR-YMnO3-hpkot")
        self.assertTrue(filecmp.cmp("KPOINTS-dielectric_function-YMnO3",
                                    "KPOINTS"))

    def test_competing_phase_metal(self):
        make_kpoints(task="competing_phase", poscar="PPOSCAR-Ru", is_metal=True)
        self.assertTrue(filecmp.cmp("KPOINTS-cp_expected_Ru", "KPOINTS"))

    def test_molecule(self):
        make_kpoints(task="competing_phase_molecule", poscar="POSCAR-F2")
        self.assertTrue(filecmp.cmp("KPOINTS-mol_expected_F2", "KPOINTS"))


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
