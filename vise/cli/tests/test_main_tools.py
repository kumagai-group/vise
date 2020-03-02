# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import tempfile

from vise.cli.main_tools import (
    potcar_str2dict, list2dict, dict2list, get_user_settings, get_default_args)
from vise.util.testing import ViseTest


class PotcarStr2DictTest(ViseTest):
    def test_none(self):
        self.assertEqual({}, potcar_str2dict(None))

    def test_list(self):
        expected = {"Mg": "Mg_pv", "O": "O_h"}
        actual = potcar_str2dict(["Mg_pv", "O_h"])
        self.assertEqual(expected, actual)

    def test_str(self):
        expected = {"Mg": "Mg_pv"}
        actual = potcar_str2dict("Mg_pv")
        self.assertEqual(expected, actual)

    def test_mutliple_potcars_for_same_element_error(self):
        with self.assertRaises(ValueError):
            potcar_str2dict(["Mg_pv", "Mg"])

    def test_incorrect_element_potcar_error(self):
        with self.assertRaises(ValueError):
            potcar_str2dict(["MgHH", "Mg"])


class List2DictTest(ViseTest):
    def setUp(self) -> None:
        self.key_candidates = ["ENCUT", "MAGMOM", "LWAVE"]

    def test_dict2list(self):
        flattened_list = ["ENCUT", "500", "MAGMOM", "4", "4.0", "LWAVE", "F"]
        actual = list2dict(flattened_list, self.key_candidates)
        expected = {"ENCUT": 500, "MAGMOM": [4, 4.0], "LWAVE": False}
        self.assertEqual(expected, actual)

    def test_fail(self):
        flattened_list = ["ENCAT", "500"]
        with self.assertRaises(ValueError):
            list2dict(flattened_list, self.key_candidates)

    def test_fail2(self):
        flattened_list = ["ENCUT", "500", "MAGMOM"]
        with self.assertRaises(ValueError):
            list2dict(flattened_list, self.key_candidates)


class Dict2ListTest(ViseTest):
    def test_dict2list(self):
        actual = dict2list({"a": 1, "b": "2 3 4", "c": True})
        expected = ["a", "1", "b", "2", "3", "4", "c", "True"]
        self.assertEqual(expected, actual)


class GetUserSettingsTest(ViseTest):
    def setUp(self) -> None:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml") as tmp_file:
            tmp_file.write("""# This file is needed for the unittest.
symprec: 0.01
user_incar_settings:
  ALGO: False
ldauu:
  Mn: 5
  O: 3
xc: hse
potcar_set: Mg_pv O_h
dir: ./path""")
            # This is needed to actually write down the string.
            tmp_file.flush()
            self.setting_keys = ["symprec", "user_incar_settings", "ldauu",
                                 "xc", "potcar_set", "test_key", "dir"]
            self.user_setting, path = \
                get_user_settings(tmp_file.name, self.setting_keys)

    def test_dict(self):
        actual = self.user_setting['ldauu']
        expected = {'Mn': 5, 'O': 3}
        self.assertEqual(expected, actual)

    def test_bool(self):
        actual = self.user_setting['user_incar_settings']['ALGO']
        self.assertEqual(False, actual)

    def test_float(self):
        actual = self.user_setting['symprec']
        self.assertEqual(0.01, actual)

    def test_string(self):
        actual = self.user_setting['xc']
        self.assertEqual("hse", actual)

    def test_string_list(self):
        actual = self.user_setting['potcar_set']
        self.assertEqual("Mg_pv O_h", actual)

    def test_full_path(self):
        actual = self.user_setting['dir']
        cwd = Path.cwd()
        expected = cwd / "path"
        self.assertEqual(str(expected), actual)

    def test_fail(self):
        with self.assertRaises(ValueError):
            with tempfile.NamedTemporaryFile(mode="w",
                                             suffix=".yaml") as tmp_file:
                tmp_file.write("incorrect_key: 0.01")
                # This is needed to actually write down the string.
                tmp_file.flush()
                get_user_settings(tmp_file.name, self.setting_keys)


class GetDefaultArgsTest(ViseTest):
    def setUp(self) -> None:
        def test_func(a, b=1, c=True):
            pass

        class test_cls:
            def __init__(self, d, e=2, f=None):
                pass
            def test_method(self, g, h=3.0):
                pass

        self.default_func = get_default_args(test_func)
        self.default_cls = get_default_args(test_cls)
        self.default_method = get_default_args(test_cls.test_method)

    def test_func(self):
        expected = {"b": 1, "c": True}
        self.assertEqual(expected, self.default_func)

    def test_cls(self):
        expected = {"e": 2, "f": None}
        self.assertEqual(expected, self.default_cls)

    def test_method(self):
        expected = {"h": 3.0}
        self.assertEqual(expected, self.default_method)