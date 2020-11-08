#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys
from pathlib import Path

from vise import __version__
from vise.analyzer.atom_grouping_type import AtomGroupingType
from vise.cli.main_functions import get_poscar_from_mp, VaspSet, plot_band, \
    plot_dos, band_edge_properties
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.str_related_tools import str2bool

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

    # -- parent parser: vasprun
    vasprun_parser = argparse.ArgumentParser(description="", add_help=False)
    vasprun_parser.add_argument(
        "-v", "--vasprun", type=Path, default=defaults.vasprun,
        help="vasprun.xml file name.")
    # -- parent parser: outcar
    outcar_parser = argparse.ArgumentParser(description="", add_help=False)
    outcar_parser.add_argument(
        "-o", "--outcar", type=Path, default=defaults.outcar,
        help="OUTCAR file name.")

    # -- get_poscar -----------------------------------------------------------
    parser_get_poscar = subparsers.add_parser(
        name="get_poscar",
        description="Tools for generating a POSCAR file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['gp'])

    parser_get_poscar.add_argument(
        "-p", "--poscar", type=str, default="POSCAR",
        help="POSCAR-type file name.", metavar="FILE")
    parser_get_poscar.add_argument(
        "-pi", "--prior_info", type=Path, default=Path("prior_info.yaml"),
        help="prior_info.yaml file name.", metavar="FILE")
    parser_get_poscar.add_argument(
        "-m", "--mpid", type=str, required=True,
        help="MP entry id with prefix, e.g., mp-1234.")

    parser_get_poscar.set_defaults(func=get_poscar_from_mp)

    # -- vasp_set ---------------------------------------------------------
    parser_vasp_set = subparsers.add_parser(
        name="vasp_set",
        description="Tools for constructing vasp input set",
        parents=[vasprun_parser, outcar_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vs'])

    parser_vasp_set.add_argument(
        "-p", "--poscar", default="POSCAR", type=Path,
        help="POSCAR-type input structure file name.")
    parser_vasp_set.add_argument(
        "-t", "--task", default=defaults.task, type=Task,
        choices=Task.name_list(),
        help=f"Task name from {Task.names_string()}.")
    parser_vasp_set.add_argument(
        "-x", "--xc", default=defaults.xc, type=Xc,
        choices=Xc.name_list(),
        help=f"Exchange-correlation (XC) interaction from {Xc.names_string()}.")
    parser_vasp_set.add_argument(
        # Default must be None to use the insulator_kpoint_density.
        # defaults.kpoints_density is used at structure_kpoints_generator.py
        "-k", "--kpt_density",
        type=float,
        help="K-point density in Angstrom along each direction. "
             f"Defaults of metal and insulators are {defaults.kpoint_density}, "
             f"and {defaults.insulator_kpoint_density}.")
    parser_vasp_set.add_argument(
        "--potcar",
        dest="overridden_potcar", default=defaults.overridden_potcar,
        type=str, nargs="+",
        help="Additional User specifying POTCAR set. E.g., Mg_pv O_h")
    parser_vasp_set.add_argument(
        "-uis", "--user_incar_settings", type=str, nargs="+",
        help="Additional user_incar_settings in make_input classmethod of "
             "ViseInputSet in vise. The default of this flag is set by the "
             "vise.yaml, so if one does not want to override the default, use "
             "additional_user_incar_setting instead. See also document in "
             "vise input_set for details.")
    parser_vasp_set.add_argument(
        "-d", "--prev_dir", type=Path,
        help="Parse the previous calculations for better input settings.")
    parser_vasp_set.add_argument(
        "--options", type=str, nargs="+",
        help="Additional Keyword args for options in make_input classmethod of "
             "ViseInputSet in vise. See document in vise for details.")
    parser_vasp_set.add_argument(
        "--uniform_kpt_mode", action="store_true",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive cell.")
    parser_vasp_set.add_argument(
        "--file_transfer_type", type=str, nargs="+",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive cell.")

    parser_vasp_set.set_defaults(func=VaspSet)

    # -- plot_band -------------------------------------------------------------
    parser_plot_band = subparsers.add_parser(
        name="plot_band",
        description="Tools for plotting band structures",
        parents=[vasprun_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pb'])

    parser_plot_band.add_argument(
        "-k", "--kpoints", dest="kpoints_filename", default="KPOINTS", type=str,
        help="Kpoints file name.")
    parser_plot_band.add_argument(
        "-y", "--y_range", nargs="+", default=[-10, 10], type=float,
        help="Energy range, requiring two values.")
    parser_plot_band.add_argument(
        "-f", "--filename", type=str, default="band.pdf", help="pdf file name.")

    parser_plot_band.set_defaults(func=plot_band)

    # -- plot_dos -----------------------------------------------------------
    parser_plot_dos = subparsers.add_parser(
        name="plot_dos",
        description="Tools for plotting density of states",
        parents=[vasprun_parser, outcar_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pd'])

    parser_plot_dos.add_argument(
        "-t", "--type",
        type=AtomGroupingType,
        default=AtomGroupingType.non_equiv_sites,
        choices=AtomGroupingType.name_list(),
        help="How to group atoms for pdos.")
    parser_plot_dos.add_argument(
        "-l", "--legend", type=str2bool, default=True,
        help="Whether to show the figure legend.")
    parser_plot_dos.add_argument(
        "-c", "--crop_first_value", type=str2bool, default=True,
        help="Whether to crop the first value in DOS.")
    parser_plot_dos.add_argument(
        "-x", "--x_range", nargs="+", type=float,
        help="Set minimum and maximum energies.")
    parser_plot_dos.add_argument(
        "-y", "--y_max_ranges", nargs="+", type=float,
        help="Set max energies for each plot. Assertion error is raised when"
             "the number of ranges does not fit to the number of figures.")
    parser_plot_dos.add_argument(
        "--target", type=str, nargs="+", help="""
        Show specific PDOS. The input depends on AtomGroupingType.\n
        AtomGroupingType.atoms: ["1", "2"] \n
        AtomGroupingType.elements: ["Mg", "O"] 
        """)
    parser_plot_dos.add_argument(
        "-f", "--filename", type=str, default="dos.pdf", help="Pdf file name.")
    parser_plot_dos.add_argument(
        "-b", "--base_energy", type=float,
        help="Set when showing the figure in the absolute energies scale.")
    parser_plot_dos.set_defaults(func=plot_dos)

    # -- band_edge -------------------------------------------------------------
    parser_band_edge = subparsers.add_parser(
        name="band_edge",
        description="Calculate the band edge properties",
        parents=[vasprun_parser, outcar_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['be'])

    parser_band_edge.set_defaults(func=band_edge_properties)

    try:
        import argcomplete
        argcomplete.autocomplete(parser)
        # https://github.com/kislyuk/argcomplete#activating-global-completion
        # This supports bash/zsh autocompletion. To enable this, pip install
        # argcomplete, activate global completion, or add
        #      eval "$(register-python-argcomplete vise)"
        # into your .bash_profile or .bashrc
    except ImportError:
        pass

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

