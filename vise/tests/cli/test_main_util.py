# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from pathlib import Path

from pymatgen.core import Element
from vise.analyzer.vasp.handle_volumetric_data import default_border_fractions
from vise.cli.main_util import parse_args

parent_dir = Path(__file__).parent


def test_make_atom_poscars_wo_options():
    parsed_args = parse_args(["map"])
    expected = Namespace(
        dirname=Path.cwd(),
        elements=None,
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_atom_poscars_w_options():
    parsed_args = parse_args(["map", "-d", "a", "-e", "H", "He"])
    expected = Namespace(
        dirname=Path("a"),
        elements=[Element.H, Element.He],
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_spin_decomposed_volumetric_files():
    parsed_args = parse_args(["sdvf", "-c", "CHGCAR"])
    expected = Namespace(chgcar="CHGCAR", func=parsed_args.func)
    assert parsed_args == expected


def test_make_light_weight_vol_data_wo_options():
    parsed_args = parse_args(["lwvd", "-v", "CHGCAR"])
    expected = Namespace(volumetric_file="CHGCAR",
                         output_lw_volumetric_filename=None,
                         border_fractions=default_border_fractions,
                         output_vesta_filename=None,
                         original_vesta_file=None,
                         func=parsed_args.func)
    assert parsed_args == expected


def test_make_light_weight_vol_data_w_options():
    parsed_args = parse_args(["lwvd",
                             "-v", "AECCAR",
                              "-o", "AECCAR_lw",
                              "-b", "0.1",
                              "-out_vesta", "out.vesta",
                              "-orig_vesta", "orig.vesta"])
    expected = Namespace(volumetric_file="AECCAR",
                         output_lw_volumetric_filename=Path("AECCAR_lw"),
                         border_fractions=[0.1],
                         output_vesta_filename=Path("out.vesta"),
                         original_vesta_file=Path("orig.vesta"),
                         func=parsed_args.func)
    assert parsed_args == expected
