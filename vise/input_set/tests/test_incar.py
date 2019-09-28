# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.input_set.incar import incar_flags, ViseIncar

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class IncarFlagsTest(ViseTest):

    def test(self) -> None:
        expected = ["SMASS", "NBLOCK"]
        actual = incar_flags["md"]
        self.assertEqual(expected, actual)



class ViseIncarTest(ViseTest):
    pass
    # def test_ModIncar(self):
    #     i = ModIncar.from_file("INCAR-ModIncar_before")
    #     i.pretty_write_file("INCAR-ModIncar_actual")
    #     self.assertTrue(
#           filecmp.cmp("INCAR-ModIncar_expected", "INCAR-ModIncar_actual"))

# #
# class MakeIncarTest(unittest.TestCase):
#     #
#     # # TODO: write better tests
#     # def test_structure_opt(self):
#     #     MakeIncar(task="structure_opt", functional="hse",
#     #               defect_in="defect.in", hfscreen=0.2, aexx=0.3,
#     #               is_magnetization=True, ldau=True,
#                   my_incar_setting="my_INCAR_setting.yaml")
#     #
#     # def test_band(self):
#     #     MakeIncar(task="band", functional="pbe", poscar="PPOSCAR-YMnO3",
#                   ldau=True, my_incar_setting="my_INCAR_setting.yaml")
#     #
#     # def test_dos(self):
#     #     MakeIncar(task="dos", functional="pbe", poscar="PPOSCAR-YMnO3",
#                   ldau=True, is_magnetization=True)
#     #
#     # def test_dos(self):
#        MakeIncar(task="dielectric_function", functional="pbe",
#                   # poscar="PPOSCAR-YMnO3", ldau=True, is_magnetization=True)

