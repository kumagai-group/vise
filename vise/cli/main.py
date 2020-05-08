#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys
from pathlib import Path

from vise import __version__
from vise.cli.main_functions import get_poscar_from_mp, VaspSet
from vise.cli.main_tools import potcar_str2dict
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.tools import str2bool

logger = get_logger(__name__)


def parse_args(args):
    vise_yaml_files = '\n'.join(["* " + str(f) for f in defaults.yaml_files])
    parser = argparse.ArgumentParser(
        description=f"""
Vise is a package that helps researchers to do first-principles calculations
with the VASP code.

Author: Yu Kumagai
Version: {__version__}
    """,
        epilog=f"The parsed vise.yaml files are:\n{vise_yaml_files}",
        formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers()

    # -- get_poscars -----------------------------------------------------------
    parser_get_poscar = subparsers.add_parser(
        name="get_poscars",
        description="Tools for generating a POSCAR file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['gp'])

    parser_get_poscar.add_argument(
        "-p", "--poscar", type=str, default="POSCAR",
        help="POSCAR-type file name.", metavar="FILE")
    parser_get_poscar.add_argument(
        "-m", "--mpid", type=str,
        help="MP entry id without prefix, e.g., mp-1234.")

    parser_get_poscar.set_defaults(func=get_poscar_from_mp)

    # -- vasp_set ---------------------------------------------------------
    parser_vasp_set = subparsers.add_parser(
        name="vasp_set",
        description="Tools for constructing vasp input set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vs'])

    parser_vasp_set.add_argument(
        "-p", "--poscar",
        default="POSCAR",
        type=str,
        help="POSCAR-type input structure file name.")
    parser_vasp_set.add_argument(
        "-t", "--task",
        default=defaults.task,
        type=Task,
        choices=[t for t in Task],
        help=f"Task name from {Task.name_list()}.")
    parser_vasp_set.add_argument(
        "-x", "--xc",
        default=defaults.xc,
        type=Xc,
        choices=[xc for xc in Xc],
        help="Exchange-correlation (XC) interaction from {Xc.name_list()}.")
    parser_vasp_set.add_argument(
        "-k", "--kpt_density",
        default=defaults.kpoint_density,
        type=float,
        help="K-point density in Angstrom along each direction .")
    parser_vasp_set.add_argument(
        "--potcar",
        dest="overridden_potcar",
        default=defaults.overridden_potcar,
        type=str,
        nargs="+",
        help="Additional User specifying POTCAR set. E.g., Mg_pv O_h")
    parser_vasp_set.add_argument(
        "-c", "--charge",
        type=float,
        default=0.0,
        help="Charge state.")
    parser_vasp_set.add_argument(
        "-uis", "--user_incar_settings",
        type=str,
        nargs="+",
        help="Addtional user_incar_settings in make_input classmethod of ViseInputSet "
             "in vise. The default of this flag is set by the vise.yaml, "
             "so if one does not want to override the default, use "
             "additional_user_incar_setting instead. See also document in "
             "vise input_set for details.")
    parser_vasp_set.add_argument(
        "-d", "--prev_dir",
        type=Path,
        help="Parse the previous calculations for better input settings.")
    parser_vasp_set.add_argument(
        "-o", "--options",
        type=str,
        nargs="+",
        help="Additional Keyword args for options in make_input classmethod of "
             "ViseInputSet in vise. See document in vise for details.")
    parser_vasp_set.add_argument(
        "--uniform_kpt_mode",
        action="store_true",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive cell.")
    parser_vasp_set.add_argument(
        "--file_transfer_type",
        type=str,
        nargs="+",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive cell.")

    parser_vasp_set.set_defaults(func=VaspSet)


#     # # try:
#     # #     import argcomplete
#     # #     argcomplete.autocomplete(parser)
#     # #     # This supports bash autocompletion. To enable this, pip install
#     # #     # argcomplete, activate global completion, or add
#     # #     #      eval "$(register-python-argcomplete vise)"
#     # #     # into your .bash_profile or .bashrc
#     # # except ImportError:
#     # #     pass
#     #
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

