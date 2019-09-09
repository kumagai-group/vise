# -*- coding: utf-8 -*-

from pymatgen.io.vasp import Potcar
from vise.input_set.sets import (
    load_potcar_yaml, load_default_incar_settings, check_keys, nelect, nbands,
    XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS, TaskStructureKpoints, XcTaskPotcar,
    ALL_FLAGS, OTHER_FLAGS, Task, TaskIncarSettings, XcIncarSettings, Xc)
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class LoadPotcarYamlTest(ViseTest):
    def setUp(self):
        self.normal = load_potcar_yaml("normal")
        self.normal_override = load_potcar_yaml(
            set_name="normal", override_potcar_set={"Zr": "Zr"})

    def test(self):
        self.assertEqual("Zr_sv", self.normal["Zr"])

    def test_override(self):
        self.assertEqual("Zr", self.normal_override["Zr"])

    def test_raise_error(self):
        with self.assertRaises(KeyError):
            load_potcar_yaml("not_exist_set_name")


class LoadDefaultIncarSettingsTest(ViseTest):
    def setUp(self):
        self.hse = \
            load_default_incar_settings(yaml_filename="xc_incar_set.yaml",
                                        required_flags=XC_REQUIRED_FLAGS,
                                        optional_flags=XC_OPTIONAL_FLAGS,
                                        key_name="hse")

    def test(self):
        expected = {'LWAVE': True,
                    'ALGO': 'D',
                    'PRECFOCK': 'Fast',
                    'LHFSCREEN': 0.208,
                    'LHFCALC': True,
                    'TIME': 0.4}
        self.assertEqual(expected, self.hse)


class CheckKeysTest(ViseTest):
    def test_succeed(self):
        d = {"a": 1, "b": 2, "c": 3}
        self.assertTrue(check_keys(d=d, required={"a"}, optional={"b", "c"}))

    def test_raise_error1(self):
        d = {"a": 1, "b": 2, "c": 3}
        with self.assertRaises(KeyError):
            check_keys(d=d, required={"a", "d"}, optional={"b", "c"})

    def test_raise_error2(self):
        d = {"a": 1, "b": 2, "c": 3}
        with self.assertRaises(KeyError):
            check_keys(d=d, required={"a"}, optional={"b"})


class NelectTest(ViseTest):
    def setUp(self):
        self.mgo = self.get_structure_by_name("MgO")
        self.potcar = Potcar(["Mg", "O"], functional="PBE_54")

    def test_neutral(self):
        expected = 8
        self.assertEqual(expected, nelect(self.mgo, self.potcar))

    def test_charged(self):
        expected = 9
        self.assertEqual(expected, nelect(self.mgo, self.potcar, charge=-1))


class NbandsTest(ViseTest):
    def setUp(self):
        self.mgo = self.get_structure_by_name("MgO")
        self.potcar = Potcar(["Mg", "O"], functional="PBE_54")

    def test_neutral(self):
        expected = 8 / 2 + 4 + 4
        self.assertEqual(expected, nbands(self.mgo, self.potcar))


class TaskStructureKpointsTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO")
        self.structure_opt = \
            TaskStructureKpoints.from_options(task=Task.band,
                                              original_structure=mgo,
                                              standardize_structure=True,
                                              sort_structure=True,
                                              is_magnetization=False,
                                              kpt_mode="primitive",
                                              kpt_density=3,
                                              kpt_shift=[0, 0, 0],
                                              only_even=True,
                                              band_ref_dist=0.025,
                                              factor=1,
                                              symprec=0.01,
                                              angle_tolerance=3)

    def test(self):
        print(self.structure_opt.kpoints)


class XcTaskPotcarTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO")
        self.hse = \
            XcTaskPotcar.from_options(xc=Xc.hse, task=Task.band, structure=mgo)

    def test(self):
        print(self.hse.max_enmax)


class FlagsTest(ViseTest):
    def test(self):
        print(len(OTHER_FLAGS))
        print(len(ALL_FLAGS))


class TaskIncarSettingsTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO64atoms")
        potcar = Potcar(["Mg", "O"], functional="PBE_54")
        self.defect = TaskIncarSettings.from_options(task=Task.defect,
                                                     structure=mgo,
                                                     potcar=potcar,
                                                     num_kpoints=10,
                                                     max_enmax=400.0,
                                                     is_magnetization=False,
                                                     band_gap=5.0,
                                                     vbm_cbm=[2.0, 7.0],
                                                     npar_kpar=False,
                                                     encut=None,
                                                     structure_opt_encut_factor=1.3)

    def test(self):
        print(self.defect.settings)


class XcIncarSettingsTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO64atoms")
        self.pbe0 = XcIncarSettings.from_options(xc=Xc.pbe0,
                                                 structure=mgo,
                                                 factor=3,
                                                 aexx=0.25)

    def test(self):
        print(self.pbe0.settings)
