# -*- coding: utf-8 -*-

from pymatgen.io.vasp import Potcar

from vise.input_set.settings_util import (
    load_potcar_yaml, load_default_incar_settings, check_keys, nelect, nbands)
from vise.input_set.settings_incar import (
    XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
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
                    'HFSCREEN': 0.208,
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
        self.mgo = self.get_structure_by_name("MgO").composition
        self.potcar = Potcar(["Mg", "O"], functional="PBE_54")

    def test_neutral(self):
        expected = 8
        self.assertEqual(expected, nelect(self.mgo, self.potcar))

    def test_charged(self):
        expected = 9
        self.assertEqual(expected, nelect(self.mgo, self.potcar, charge=-1))


class NbandsTest(ViseTest):
    def setUp(self):
        self.mgo = self.get_structure_by_name("MgO").composition
        self.potcar = Potcar(["Mg", "O"], functional="PBE_54")

    def test_neutral(self):
        expected = 8 / 2 + 4 + 4
        self.assertEqual(expected, nbands(self.mgo, self.potcar))

