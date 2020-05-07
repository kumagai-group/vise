# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os
import shutil

from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from vise.cli.main import parse_args
from vise.util.testing import ViseTest
from vise.defaults import defaults

parent_dir = Path(__file__).parent


def test_get_poscars_wo_options():
    actual = parse_args(["gp"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        poscar="POSCAR",
        mpid=None,
        func=actual.func)
    assert actual == expected


def test_get_poscars_w_options():
    parsed_args = parse_args(["gp",
                              "-p", "a",
                              "-m", "123"])

    expected = Namespace(
        poscar="a",
        mpid=123,
        func=parsed_args.func)
    assert parsed_args == expected


def test_vasp_set_wo_options():
    parsed_args = parse_args(["vs"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        additional_user_incar_settings=None,
        angle_tol=5.0,
        charge=0.0,
        dirs=["."],
        json="vise.json",
        kpt_density=defaults.kpoint_density,
        length_tol=defaults.symmetry_length_tolerance,
        ldaul=defaults.ldaul,
        ldauu=defaults.ldauu,
        potcar=defaults.potcar_dict,
        potcar_set_name=defaults.potcar_set_name,
        poscar="POSCAR",
        prev_dir=None,
        print=False,
        prior_info=True,
        standardize_structure=True,
        task=defaults.task,
        user_incar_settings=defaults.user_incar_settings,
        vise_opts={},
        xc=defaults.xc,
        func=parsed_args.func)
    assert parsed_args == expected


def test_vasp_set_w_options():
    parsed_args = parse_args(["vs",
                              "--print",
                              "--potcar", "Mg_pv", "O_h",
                              "--potcar_set_name", "gw",
                              "-x", "pbesol",
                              "-t", "band",
                              "--vise_opts", "encut", "800",
                              "--user_incar_settings", "LREAD", "F",
                              "-auis", "ALGO", "D",
                              "--ldauu", "Mg", "5",
                              "--ldaul", "Mg", "1",
                              "-c", "10",
                              "--length_tol", "0.0001",
                              "--angle_tol", "20.0",
                              "--json", "tmp.json",
                              "--poscar", "POSCAR-tmp",
                              "-k", "4.2",
                              "-s", "F",
                              "-pi", "T",
                              "--dirs", "a", "b",
                              "-d", "c"])

    expected = Namespace(
        additional_user_incar_settings=["ALGO", "D"],
        angle_tol=20.0,
        charge=10.0,
        dirs=["a", "b"],
        json="tmp.json",
        kpt_density=4.2,
        length_tol=0.0001,
        ldauu=["Mg", "5"],
        ldaul=["Mg", "1"],
        potcar=["Mg_pv", "O_h"],
        potcar_set_name="gw",
        poscar="POSCAR-tmp",
        prev_dir="c",
        print=True,
        prior_info=True,
        standardize_structure=False,
        task="band",
        user_incar_settings=["LREAD", "F"],
        vise_opts=["encut", "800"],
        xc="pbesol",
        func=parsed_args.func)

    assert parsed_args == expected


