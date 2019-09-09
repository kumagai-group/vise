# -*- coding: utf-8 -*-

from vise.input_set.input_set import load_potcar_yaml
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


