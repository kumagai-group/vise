# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from pathlib import Path

from pymatgen.core import Element
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


