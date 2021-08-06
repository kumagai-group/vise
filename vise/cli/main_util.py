#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.core import Element, Structure
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from vise.analyzer.vasp.handle_volumetric_data import default_border_fractions
from vise.cli.main import description, epilog
from vise.cli.main_util_functions import make_atom_poscars, \
    make_spin_decomposed_volumetric_files, make_light_weight_vol_data, \
    make_phonon_poscars, make_phonon_figs
from vise.util.logger import get_logger

logger = get_logger(__name__)


warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers()

    # -- make_phonon_poscars ---------------------------------------------------
    parser_make_phonon_poscars = subparsers.add_parser(
        name="make_phonon_poscars",
        description="Tools for .",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mpp'])

    parser_make_phonon_poscars.add_argument(
        "-u", "--unitcell", type=Structure.from_file,
        help="")

    parser_make_phonon_poscars.set_defaults(func=make_phonon_poscars)

    # -- make_phonon_figs ---------------------------------------------------
    parser_make_phonon_figs = subparsers.add_parser(
        name="make_phonon_figs",
        description="Tools for .",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mpf'])

    parser_make_phonon_figs.add_argument(
        "-pi", "--phonopy_input", type=loadfn, help="")
    parser_make_phonon_figs.add_argument(
        "-vn", "--vasprun_name", type=str, help="")

    parser_make_phonon_figs.set_defaults(func=make_phonon_figs)

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

    # -- make_spin_decomposed_volumetric_files ---------------------------------
    parser_spin_decomposed_volumetric_files = subparsers.add_parser(
        name="spin_decomposed_volumetric_files",
        description="Tools for making spin decomposed volumetric files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sdvf'])

    parser_spin_decomposed_volumetric_files.add_argument(
        "-c", "--chgcar", type=str, required=True,
        help="CHGCAR-type file name with data for spin-down channel.")
    parser_spin_decomposed_volumetric_files.set_defaults(
        func=make_spin_decomposed_volumetric_files)

    # -- make_light_weight_vol_data --------------------------------------------
    parser_light_weight_vol_data = subparsers.add_parser(
        name="light_weight_vol_data",
        description="Tools for making light-weigh volumetric data file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['lwvd'])

    parser_light_weight_vol_data.add_argument(
        "-v", "--volumetric_file", type=str, required=True,
        help="CHGCAR-type volumetric file name.")
    parser_light_weight_vol_data.add_argument(
        "-o", "--output_lw_volumetric_filename", type=Path,
        help="Created file name.")
    parser_light_weight_vol_data.add_argument(
        "-b", "--border_fractions", type=float, nargs="*",
        default=default_border_fractions,
        help="Fractions that are used for depicting isosurfaces.")
    parser_light_weight_vol_data.add_argument(
        "-out_vesta", "--output_vesta_filename", type=Path,
        help="Output VESTA file name including volumetric file name and "
             "ISURFS.")
    parser_light_weight_vol_data.add_argument(
        "-orig_vesta", "--original_vesta_file", type=Path,
        help="Original VESTA file name to which volumetric file name and "
             "ISURFS are inserted.")

    parser_light_weight_vol_data.set_defaults(
        func=make_light_weight_vol_data)

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

