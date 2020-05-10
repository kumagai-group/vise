# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.cli.main_tools import (
    potcar_str2dict, list2dict)
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


#
# class GetDefaultArgsTest(ViseTest):
#     def setUp(self) -> None:
#         def test_func(a, b=1, c=True):
#             pass
#
#         class test_cls:
#             def __init__(self, d, e=2, f=None):
#                 pass
#             def test_method(self, g, h=3.0):
#                 pass
#
#         self.default_func = get_default_args(test_func)
#         self.default_cls = get_default_args(test_cls)
#         self.default_method = get_default_args(test_cls.test_method)
#
#     def test_func(self):
#         expected = {"b": 1, "c": True}
#         self.assertEqual(expected, self.default_func)
#
#     def test_cls(self):
#         expected = {"e": 2, "f": None}
#         self.assertEqual(expected, self.default_cls)
#
#     def test_method(self):
#         expected = {"h": 3.0}
#         self.assertEqual(expected, self.default_method)