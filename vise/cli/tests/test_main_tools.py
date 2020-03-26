# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
from pathlib import Path
import tempfile
from unittest.mock import patch

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

    @patch("pathlib.Path.home")
    def test(self, mock_home) -> None:
        with tempfile.TemporaryDirectory() as dirname:
            mock_home.return_value = Path(dirname)
            os.chdir(dirname)
            os.mkdir(".home_hidden_dir")
            with open('.home_hidden_dir/test_vise.yaml', mode="w") as f:
                f.write("""# This file is needed for the unittest.
a_float: 0.01
settings:
  b: False
c_dict:
  d: 5
  e: 3
f_str: g h
i_dir: ~/j""")

            with self.assertRaises(KeyError):
                get_user_settings(yaml_filename="test_vise.yaml",
                                  setting_keys=["nokey"],
                                  home_hidden_directory=".home_hidden_dir")

            setting_keys = ["a_float", "settings", "c_dict", "f_str", "i_dir"]

            user_setting, _ = \
                get_user_settings(yaml_filename="test_vise.yaml",
                                  setting_keys=setting_keys,
                                  home_hidden_directory=".home_hidden_dir")

            self.assertEqual(0.01, user_setting["a_float"])
            self.assertEqual({'b': False}, user_setting["settings"])
            self.assertEqual({'d': 5, 'e': 3}, user_setting["c_dict"])
            self.assertEqual("g h", user_setting["f_str"])
            expected = str(mock_home.return_value / ".home_hidden_dir" / "~/j")
            self.assertEqual(expected, user_setting["i_dir"])

            with open('test_vise.yaml', mode="w") as f:
                f.write("""a_float: 0.02,
f_str: x y
""")

            os.makedirs("dir_a/dir_b")

            with open('dir_a/dir_b/test_vise.yaml', mode="w") as f:
                f.write("""a_float: 0.03""")

            os.chdir("dir_a/dir_b")
            user_setting, yaml_files = \
                get_user_settings(yaml_filename="test_vise.yaml",
                                  setting_keys=setting_keys,
                                  home_hidden_directory=".home_hidden_dir")
            self.assertEqual(0.03, user_setting["a_float"])
            self.assertEqual("x y", user_setting["f_str"])


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