# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose
from unittest.mock import patch

from vise.util.testing import ViseTest
from vise.cli.main_function import (
    vasp_setting_from_args, get_poscar_from_mp, vasp_set)
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger


# class TestGet_poscar_from_mp(TestCase):
#     def test_get_poscar_from_mp(self):
#         self.fail()

logger = get_logger(__name__)


prec = {"symprec": 0.1, "angle_tolerance": 5}
vasp = {"task": "band",
        "xc": "hse",
        "ldauu": ["Mn", "5", "O", "3"],
        "ldaul": ["Mn", "3", "O", "2"],
        "vise_opts": ["only_even", "True"],
        "user_incar_settings": ["POTIM", "0.4"],
        "additional_user_incar_settings": ["ALGO", "Fast"],
        "potcar_set": ["Mn_pv"],
        "potcar_set_name": "normal"}



class TestUpdateVaspSet(ViseTest):
    def setUp(self) -> None:
        self.args = Namespace(**vasp)

    def test(self):
        user_incar_settings, vis_base_kwargs, task, xc = \
            vasp_setting_from_args(self.args)

        expected = {"POTIM": 0.4, "ALGO": "Fast"}
        self.assertEqual(expected, user_incar_settings)

        expected = {"ldauu": {"Mn": 5, "O": 3},
                    "ldaul": {"Mn": 3, "O": 2},
                    "only_even": True,
                    "override_potcar_set": {"Mn": "Mn_pv"},
                    "potcar_set_name": "normal"}
        self.assertEqual(expected, vis_base_kwargs)
        self.assertEqual(task, Task.band)
        self.assertEqual(xc, Xc.hse)


class TestGetPoscarFromMp(ViseTest):
    def setUp(self) -> None:
        self.args_num = Namespace(number=1234)
        self.kwargs = {"elements": ["Mg", "O"],
                       "e_above_hull": 1.5,
                       "molecules": True}
        self.args_elements = Namespace(**self.kwargs)
        self.args_none = Namespace()

    # @patch('MPRester().get_structure_by_material_id')
    # def test_number(self):
    #     self.args_num

    @patch('vise.cli.main_function.make_poscars_from_mp')
    def test_mp_poscars(self, mock_make_poscars):
        get_poscar_from_mp(self.args_elements)
        mock_make_poscars.assert_called_with(**self.kwargs)

    def test_warning(self):
        expected = ['WARNING:vise.cli.main_function:Set mp number or elements']
        with self.assertLogs(level='WARNING') as cm:
            get_poscar_from_mp(self.args_none)
            self.assertEqual(expected, cm.output)


class TestVaspSet(ViseTest):
    def setUp(self) -> None:
        min_default_kwargs = {"poscar": "POSCAR",
                              "kpt_density": 2.5,
                              "standardize": True,
                              "charge": 0,

                              }
        self.args_minimum = Namespace(**min_default_kwargs, **prec, **vasp)

#        option_kwargs = {"prev_dir": "../"
#                         }

