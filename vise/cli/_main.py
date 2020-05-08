#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys

from vise import __version__
from vise.cli.main_function import get_poscar_from_mp, vasp_set
from vise.defaults import defaults
from vise.input_set.datasets.dataset_util import PotcarSet
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.tools import str2bool

logger = get_logger(__name__)

#def parse_args(args):



# def parse_args(args):
#     vise_yaml_files = '\n'.join(["* " + str(f) for f in defaults.yaml_files])
#     parser = argparse.ArgumentParser(
#         description=f"""
# Vise is a package that helps researchers to do first-principles calculations
# with the VASP code.
#
# Author: Yu Kumagai
# Version: {__version__}
#     """,
#         epilog=f"The parsed vise.yaml files are:\n{vise_yaml_files}",
#         formatter_class=argparse.RawTextHelpFormatter)
#
#     subparsers = parser.add_subparsers()
#
#     # -- get_poscars -----------------------------------------------------------
#     parser_get_poscar = subparsers.add_parser(
#         name="get_poscars",
#         description="Tools for generating a POSCAR file.",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         aliases=['gp'])
#
#     parser_get_poscar.add_argument(
#         "-p", "--poscar", type=str, default="POSCAR",
#         help="POSCAR-type file name.", metavar="FILE")
#     parser_get_poscar.add_argument(
#         "-m", "--mpid", type=int,
#         help="MP entry id without prefix, e.g., mp-1234.")
#
#     parser_get_poscar.set_defaults(func=get_poscar_from_mp)
#
#     # -- vasp_set ---------------------------------------------------------
#     parser_vasp_set = subparsers.add_parser(
#         name="vasp_set",
#         description="Tools for constructing vasp input set with vise",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         aliases=['vs'])
#
#     parser_vasp_set.add_argument(
#         "--print", action="store_true",
#         help="Whether to print vise.json.")
#     parser_vasp_set.add_argument(
#         "--length_tol", type=float, default=defaults.symmetry_length_tolerance,
#         help="Set length tolerance used for symmetry analysis [A].")
#     parser_vasp_set.add_argument(
#         "--angle_tol", type=float,
#         default=defaults.symmetry_angle_tolerance,
#         help="Set angle tolerance used for symmetry analysis [deg.].")
#     parser_vasp_set.add_argument(
#         "--json", type=str, default="vise.json",
#         help="Json file name for printing the ViseInputSet info.")
#     parser_vasp_set.add_argument(
#         "-p", "--poscar", default="POSCAR", type=str,
#         help="POSCAR-type file name.")
#     parser_vasp_set.add_argument(
#         "-k", "--kpt_density", default=defaults.kpoint_density, type=float,
#         help="K-point density in Angstrom along each direction .")
#     parser_vasp_set.add_argument(
#         "-s", "--standardize_structure", type=str2bool, default=True,
#         help="Store if one doesn't want the cell to be transformed to a "
#              "primitive cell.")
#     parser_vasp_set.add_argument(
#         "-pi", "--prior_info", type=str2bool, default=True,
#         help="Whether to use prior_info.json when it exists.")
#     parser_vasp_set.add_argument(
#         "--dirs", nargs="+", type=str, default=["."],
#         help="Make vasp set for the directories in the same condition.")
#     parser_vasp_set.add_argument(
#          "-d", "--prev_dir", type=str,
#         "file_transfer_type"
#         help="Parse the previous calculations for better input settings.")
#     parser_vasp_set.add_argument(
#         "--potcar", default=defaults.potcar_dict, type=str, nargs="+",
#         help="User specifying POTCAR set. E.g., Mg_pv O_h")
#     parser_vasp_set.add_argument(
#         "--potcar_set_name",
#         default=defaults.potcar_set_name,
#         type=str,
#         help=f"User specifying POTCAR set, from {PotcarSet.name_list()}.")
#     parser_vasp_set.add_argument(
#         "-x", "--xc",
#         default=defaults.xc,
#         type=str,
#         choices=[e.value for e in Xc],
#         help="Specify exchange-correlation (XC) interaction treatment "
#              f"from {Xc.name_list()}.")
#     parser_vasp_set.add_argument(
#         "-t", "--task",
#         default=defaults.task,
#         type=str,
#         choices=[t.value for t in Task],
#         help=f"Specify task name from {Task.name_list()}.")
#     parser_vasp_set.add_argument(
#         "--vise_opts",
#         type=str,
#         nargs="+",
#         default=defaults.vise_opts,
#         help="Keyword args for options in make_input classmethod of "
#              "ViseInputSet in vise. See document in vise for details.")
#     parser_vasp_set.add_argument(
#         "-uis", "--user_incar_settings",
#         type=str,
#         nargs="+",
#         default=defaults.user_incar_settings,
#         help="user_incar_settings in make_input classmethod of ViseInputSet "
#              "in vise. The default of this flag is set by the vise.yaml, "
#              "so if one does not want to override the default, use "
#              "additional_user_incar_setting instead. See also document in "
#              "vise input_set for details.")
#     parser_vasp_set.add_argument(
#         "-auis", "--additional_user_incar_settings",
#         type=str,
#         nargs="+",
#         help="Use this if one does not want to override "
#              "user_incar_settings written in the yaml file")
#     parser_vasp_set.add_argument(
#         "--ldauu",
#         type=str,
#         default=defaults.ldauu,
#         nargs="+",
#         help="LDAUU values, e.g., Mn 5 Ni 5")
#     parser_vasp_set.add_argument(
#         "--ldaul",
#         type=str,
#         default=defaults.ldaul,
#         nargs="+",
#         help="LDAUL values, e.g., Mn 2 Ni 2.")
#     parser_vasp_set.add_argument(
#         "-c", "--charge",
#         type=float,
#         default=0.0,
#         help="Charge state.")
#
#     parser_vasp_set.set_defaults(func=vasp_set)
#
#     # -- plot_band -----------------------------------------------------------
#     parser_plot_band = subparsers.add_parser(
#         name="plot_band",
#         description="Tools for plotting band structures",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         aliases=['pb'])
#
#     parser_plot_band.add_argument(
#         "-v", "--vasprun", type=str, default=defaults.vasprun,
#         help="vasprun.xml file name.")
#     parser_plot_band.add_argument(
#         "-v2", dest="vasprun2", type=str)
#     parser_plot_band.add_argument(
#         "-k", dest="kpoints", default="KPOINTS", type=str)
#     parser_plot_band.add_argument(
#         "-y", dest="y_range", nargs="+", type=float,
#         help="Energy range, requiring two values.")
#     parser_plot_band.add_argument(
#         "-f", dest="filename", type=str, help="pdf file name.")
#     parser_plot_band.add_argument(
#         "-a", dest="absolute", action="store_true",
#         help="Show in the absolute energies scale.")
#     parser_plot_band.add_argument(
#         "--legend", type=str2bool, default=True,
#         help="Not show the legend.")
#
# #    parser_plot_band.set_defaults(func=plot_band)
#
#     # -- plot_dos -----------------------------------------------------------
#     parser_plot_dos = subparsers.add_parser(
#         name="plot_dos",
#         description="Tools for plotting density of states",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         aliases=['pd'])
#
#     parser_plot_dos.add_argument(
#         "-v", "--vasprun", type=str, default=defaults.vasprun,
#         help="vasprun.xml file name.")
#     parser_plot_dos.add_argument(
#         "-cv", dest="cbm_vbm", type=float, nargs="+",
#         help="Set CBM and VBM.")
#     parser_plot_dos.add_argument(
#         "-t", dest="pdos_type", type=str, default="element",
#         choices=["element", "site", "none"],
#         help="How to group atoms for pdos.")
#     parser_plot_dos.add_argument(
#         "-s", dest="specific", type=str, nargs="+",
#         help="""
#         Show specific PDOS. If list elements are integers, PDOS at particular
#         sites are shown. If elements are shown, PDOS of particular elements are
#         shown. E.g.,
#         ["1", "2"] --> At site 1 and 2 compatible with pdos_type = "none"
#         ["Mg", "O"] --> Summed at Mg and O sites compatible with pdos_type
#                         = "element"
#         """)
#     parser_plot_dos.add_argument(
#         "-o", dest="orbital", type=str2bool, default=True,
#         help="Whether to decompose orbital components.")
#     parser_plot_dos.add_argument(
#         "-x", "--x_range", nargs="+", type=float,
#         help="Set energies minimum and maximum.")
#     parser_plot_dos.add_argument(
#         "-y", "--ymaxs", nargs="+", type=float,
#         help="Set max values of y ranges. Support two ways."
#              "1st: total_max, all_the_atoms"
#              "2nd: total_max, 1st_atom, 2nd_atom, ...")
#     parser_plot_dos.add_argument(
#         "-f", dest="filename", type=str, default="dos.pdf",
#         help="Pdf file name.")
#     parser_plot_dos.add_argument(
#         "-a", "--absolute", action="store_true",
#         help="Set when showing the figure in the absolute energies scale.")
#     parser_plot_dos.add_argument(
#         "-l", "--legend", type=str2bool, default=True,
#         help="Whether to show the figure legend.")
#     parser_plot_dos.add_argument(
#         "-c", "--crop_first_value", type=str2bool, default=True,
#         help="Whether to crop the first value in DOS.")
#
# #    parser_plot_dos.set_defaults(func=plot_dos)
#
#     # -- band_gap --------------------------------------------------------------
#     parser_band_gap = subparsers.add_parser(
#         name="band_gap",
#         description="Calculate the band gap from vasprun.xml",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         aliases=['bg'])
#
#     parser_band_gap.add_argument(
#         "-v", "--vasprun", type=str, default=defaults.vasprun,
#         help="vasprun.xml file name.")
# #    parser_band_gap.set_defaults(func=band_gap)
#
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
#     return parser.parse_args(args)
#
#

# def main():
#     args = parse_args(sys.argv[1:])
#     args.func(args)


# if __name__ == "__main__":
#     main()

