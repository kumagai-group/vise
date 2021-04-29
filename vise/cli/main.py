#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys
import warnings
from pathlib import Path

from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from vise import __version__
from vise.analyzer.atom_grouping_type import AtomGroupingType
from vise.cli.main_functions import get_poscar_from_mp, VaspSet, plot_band, \
    plot_dos, band_edge_properties, plot_absorption, \
    calc_effective_mass, structure_info
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.str_related_tools import str2bool

logger = get_logger(__name__)


warnings.simplefilter('ignore', UnknownPotcarWarning)


description = """pydefect is a package that helps researchers to 
do first-principles point defect calculations with the VASP code."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
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

    # -- structure_info ---------------------------------------------------
    parser_structure_info = subparsers.add_parser(
        name="structure_info",
        description="Tools showing structure info or conventional cell.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['si'])

    parser_structure_info.add_argument(
        "-p", "--poscar", type=str, default="POSCAR",
        help="POSCAR-type file name to be read.", metavar="FILE")
    parser_structure_info.add_argument(
        "-s", "--symprec", type=float,
        default=defaults.symmetry_length_tolerance,
        help="Length tolerance for symmetry analysis.")
    parser_structure_info.add_argument(
        "-a", "--angle_tolerance", type=float,
        default=defaults.symmetry_angle_tolerance,
        help="Angle tolerance in degree for symmetry analysis.")
    parser_structure_info.add_argument(
        "-c", "--show_conventional", action="store_true",
        help="Set when showing the conventional cell.")

    parser_structure_info.set_defaults(func=structure_info)

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
        "-t", dest="task", default=defaults.task, type=Task,
        choices=Task.name_list(),
        help=f"Task name.")
    parser_vasp_set.add_argument(
        "-x", dest="xc", default=defaults.xc, type=Xc, choices=Xc.name_list(),
        help=f"Exchange-correlation (XC) interaction.")
    parser_vasp_set.add_argument(
        # Default must be None to use the insulator_kpoint_density.
        # defaults.kpoints_density is used at structure_kpoints_generator.py
        "-k", "--kpt_density", type=float,
        help="K-point density in Angstrom along each direction. "
             f"Defaults of metal and insulators are {defaults.kpoint_density}, "
             f"and {defaults.insulator_kpoint_density}.")
    parser_vasp_set.add_argument(
        "--potcar", dest="overridden_potcar",
        default=defaults.overridden_potcar, type=str, nargs="+",
        help="User specifying POTCAR set. E.g., Mg_pv O_h")
    parser_vasp_set.add_argument(
        "-uis", "--user_incar_settings", type=str, nargs="+",
        help="""Overwritten incar settings described with pairs of INCAR tag 
        names and values. These can also be set by vise.yaml, but they are 
        overwritten by this command line options.""")
    parser_vasp_set.add_argument(
        "-d", "--prev_dir", type=Path,
        help="Parse the previous calculations for better input settings.")
    parser_vasp_set.add_argument(
        "--options", type=str, nargs="+",
        help="Manual options used at CategorizedInputOptions class.")
    parser_vasp_set.add_argument(
        "--uniform_kpt_mode", action="store_true",
        help="""Kpoints with uniform k-point sampling. The k-point sampling mesh 
        and centering are determined based on the given lattice. Note that only 
        when the angles are 90 degrees, the centering is shifted along the 
        perpendicular direction. This mode is useful when calculating the 
        supercells. """)
    parser_vasp_set.add_argument(
        "--file_transfer_type", type=str, nargs="+",
        help="File transfer which is used with the prev_dir argument and can "
             "be written with pairs of file names and transfer types of "
             "m (move), c (copy) or l (link). An example is: "
             "--file_transfer_type POSCAR c WAVECAR l ")

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
        help="KPOINTS file name.")
    parser_plot_band.add_argument(
        "-y", "--y_range", nargs="+", default=[-10, 10], type=float,
        help="Energy range, requiring two values.")
    parser_plot_band.add_argument(
        "-f", "--filename", type=str, default="band.pdf", help="Pdf file name.")

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
        help="Set max energies for each plot. AssertionError is raised when"
             "the number of ranges does not fit to the number of figures.")
    parser_plot_dos.add_argument(
        "--target", type=str, nargs="+",
        help="""Show specific PDOS. The input depends on AtomGroupingType.\n 
        AtomGroupingType.atoms:   1,2 3 \n
        AtomGroupingType.elements: Mg O""")
    parser_plot_dos.add_argument(
        "-f", "--filename", type=str, default="dos.pdf", help="Pdf file name.")
    parser_plot_dos.add_argument(
        "-b", "--base_energy", type=float,
        help="Set when showing the figure in the absolute energies scale.")
    parser_plot_dos.set_defaults(func=plot_dos)

    # -- plot_absorption -------------------------------------------------------
    parser_plot_absorption = subparsers.add_parser(
        name="plot_absorption",
        description="Tools for plotting absorption coefficient",
        parents=[vasprun_parser, outcar_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pa'])

    parser_plot_absorption.add_argument(
        "-f", "--filename", type=str, default="absorption.pdf",
        help="Pdf file name.")
    parser_plot_absorption.add_argument(
        "-y", "--y_ranges", nargs="+", default=[10**3, 10**8],
        type=lambda x: 10 ** float(x),
        help="Exponential parts of base-10 for energy range in cm-1, "
             "requiring two values.")
    parser_plot_absorption.add_argument(
        "-ckk", "--calc_kk", action="store_true",
        help="Set when real part of dielectric function is explicitly "
             "calculated using the Kramers-Kronig transformation.")
    parser_plot_absorption.add_argument(
        "-i", "--ita", type=float, default=0.01,
        help="Complex shift η in the Kramers-Kronig transformation.")

    parser_plot_absorption.set_defaults(func=plot_absorption)

    # -- effective_mass -------------------------------------------------------
    parser_effective_mass = subparsers.add_parser(
        name="effective_mass",
        description="Tools for calculating effective masses",
        parents=[vasprun_parser, outcar_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['em'])

    parser_effective_mass.add_argument(
        "-t", "--temperature", type=float, default=300,
        help="Temperature in K.")
    parser_effective_mass.add_argument(
        "-c", "--concentrations", type=lambda x: 10 ** float(x), nargs="+",
        help="Exponential parts of base-10 for concentrations in cm-3")

    parser_effective_mass.set_defaults(func=calc_effective_mass)

    # -- band_edge -------------------------------------------------------------
    parser_band_edge = subparsers.add_parser(
        name="band_edge",
        description="Calculate the band edge positions",
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

