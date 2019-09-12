# -*- coding: utf-8 -*-

from pymatgen.io.vasp import Kpoints

from vise.input_set.make_kpoints import make_band_kpoints
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class MakeBandKpointsTest(ViseTest):

    def test_make_band_kpoints_mgo(self):
        structure = self.get_structure_by_name("MgO")
        kpoints = Kpoints.from_string("kpt\n 0\n M\n 4 4 4")
        k, primitive, sg = make_band_kpoints(kpoints=kpoints, structure=structure)
#        print(k)

#    def test_make_band_kpoints_ymo3_fail(self):
#        with self.assertRaises(NotAppropriatePrimitiveError):
#            make_band_kpoints(ibzkpt="IBZKPT", poscar="BPOSCAR-MgO")

#    def test_make_band_kpoints_ymo3_succeed(self):
#        make_band_kpoints(ibzkpt="IBZKPT", poscar="PPOSCAR-YMnO3-hpkot")
#       self.assertTrue(filecmp.cmp("KPOINTS-make_band-YMnO3", "KPOINTS"))


# class MakeKpointsTest(ViseTest):

    # def test_structure_opt1(self):
    #     make_kpoints(task="structure_opt", poscar="PPOSCAR-MgO")
    #     self.assertTrue(filecmp.cmp("KPOINTS-so_expected_MgO", "KPOINTS"))

    # def test_structure_opt2(self):
    #     make_kpoints(task="structure_opt", poscar="PPOSCAR-YMnO3")
    #     self.assertTrue(filecmp.cmp("KPOINTS-so_expected_YMnO3", "KPOINTS"))

    # def test_band(self):
    #     make_kpoints(task="band", poscar="PPOSCAR-YMnO3-hpkot")
    #     self.assertTrue(filecmp.cmp("KPOINTS-make_band-YMnO3", "KPOINTS"))

    # def test_dos(self):
    #     make_kpoints(task="dos", poscar="PPOSCAR-YMnO3")
    #     self.assertTrue(filecmp.cmp("KPOINTS-dos_expected_YMnO3", "KPOINTS"))

    # def test_dos2(self):
    #     make_kpoints(task="dos", poscar="PPOSCAR-ZnSnP2")
    #     self.assertTrue(filecmp.cmp("KPOINTS-dos_expected_ZnSnP2", "KPOINTS"))

    # def test_dielectric(self):
    #     make_kpoints(task="dielectric", poscar="PPOSCAR-YMnO3-hpkot")
    #     self.assertTrue(filecmp.cmp("KPOINTS-dielectric-YMnO3", "KPOINTS"))

    # def test_dielectric_function(self):
    #     make_kpoints(task="dielectric_function", poscar="PPOSCAR-YMnO3-hpkot")
    #     self.assertTrue(filecmp.cmp("KPOINTS-dielectric_function-YMnO3",
    #                                 "KPOINTS"))

    # def test_competing_phase_metal(self):
    #     make_kpoints(task="competing_phase", poscar="PPOSCAR-Ru", is_metal=True)
    #     self.assertTrue(filecmp.cmp("KPOINTS-cp_expected_Ru", "KPOINTS"))

    # def test_molecule(self):
    #     make_kpoints(task="competing_phase_molecule", poscar="POSCAR-F2")
    #     self.assertTrue(filecmp.cmp("KPOINTS-mol_expected_F2", "KPOINTS"))
