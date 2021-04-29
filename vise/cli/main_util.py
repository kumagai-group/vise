#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys
import warnings
from pathlib import Path

from pymatgen.core import Element
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from vise.cli.main import description, epilog
from vise.cli.main_util_functions import make_atom_poscars
from vise.util.logger import get_logger

logger = get_logger(__name__)


warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers()

    # -- make_atom_poscars -----------------------------------------------------
    parser_make_atom_poscars = subparsers.add_parser(
        name="make_atom_poscars",
        description="Tools for generating POSCAR files for atom calculations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['map'])

    parser_make_atom_poscars.add_argument(
        "-d", "--dirname", type=Path, default=Path.cwd(),
        help="Directory name where atom calculation directories are created.")
    parser_make_atom_poscars.add_argument(
        "-e", "--elements", type=Element, default=None, nargs="+",
        help="Element names. When not set, all atom directories are created")

    parser_make_atom_poscars.set_defaults(func=make_atom_poscars)

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

