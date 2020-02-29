#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import argparse
import sys
from typing import Union

from vise.cli.main_function import (
    get_poscar_from_mp, vasp_set, chempotdiag, plot_band, plot_dos, vasp_run,
    band_gap)
from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL, KPT_DENSITY, TIMEOUT
from vise.custodian_extension.jobs import ViseVaspJob
from vise.input_set.input_set import ViseInputSet
from vise.util.logger import get_logger
from vise.util.main_tools import dict2list, get_user_settings, get_default_args
from vise.util.mp_tools import make_poscars_from_mp
from vise.util.tools import str2bool

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

__version__ = '0.0.1dev'
__date__ = 'will be inserted'

# The following keys are set by vise.yaml
setting_keys = ["vasp_cmd",
                "symprec",
                "angle_tolerance",
                "xc",
                "kpt_density",
                "initial_kpt_density",
                "vise_opts",
                "user_incar_settings",
                "ldauu",
                "ldaul",
                "potcar_set",
                "potcar_set_name",
                "relax_iter_num",
                "removed_files",
                "timeout"]


def simple_override(d: dict, keys: Union[list, str]) -> None:
    """Override dict if keys exist in vise.yaml.

    When the value in the user_settings is a dict, it will be changed to
    list using dict2list.
    """
    user_settings, yaml_path = get_user_settings(yaml_filename="vise.yaml",
                                                 setting_keys=setting_keys)
    if isinstance(keys, str):
        keys = [keys]

    for key in keys:
        if key in user_settings:
            v = user_settings[key]
            if isinstance(v, dict):
                v = dict2list(v)
            d[key] = v


def parse_args(args):

    parser = argparse.ArgumentParser(
        description="""                            
    Vise is a package that helps researchers to do first-principles calculations 
    with the VASP code.""",
        epilog=f"""                                 
    Author: Yu Kumagai
    Version: {__version__}                                                                 
    Last updated: {__date__}""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    # -- parent parser: prec
    prec = {"symprec": SYMMETRY_TOLERANCE,
            "angle_tolerance": ANGLE_TOL}
    simple_override(prec, ["symprec", "angle_tolerance"])

    prec_parser = argparse.ArgumentParser(description="Prec-related parser",
                                          add_help=False)
    prec_parser.add_argument(
        "--symprec", type=float, default=prec["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    prec_parser.add_argument(
        "--angle_tolerance", type=float, default=prec["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")

    # -- parent parser:  vasp set
    vasp_defaults = get_default_args(ViseInputSet.make_input)
    vasp_defaults.update(ViseInputSet.TASK_OPTIONS)
    vasp_defaults.update(ViseInputSet.XC_OPTIONS)
    vasp_defaults["potcar_set"] = None
    vasp_defaults["vise_opts"] = None
    vasp_defaults["user_incar_settings"] = None
    simple_override(vasp_defaults, ["xc",
                                    "task",
                                    "vise_opts",
                                    "user_incar_settings",
                                    "potcar_set",
                                    "potcar_set_name",
                                    "ldauu",
                                    "ldaul"])

    vasp_parser = argparse.ArgumentParser(description="Vasp set-related parser",
                                          add_help=False)

    vasp_parser.add_argument(
        "--potcar", dest="potcar_set", default=vasp_defaults["potcar_set"],
        type=str, nargs="+",
        help="User specifying POTCAR set. E.g., Mg_pv O_h")
    vasp_parser.add_argument(
        "--potcar_set_name", default=vasp_defaults["potcar_set_name"], type=str,
        help="User specifying POTCAR set, i ie. normal ,gw, or mp_relax_set.")
    vasp_parser.add_argument(
        "-x", "--xc", default=str(vasp_defaults["xc"]), type=str,
        help="Exchange-correlation (XC) interaction treatment.")
    vasp_parser.add_argument(
        "-t", "--task", default=str(vasp_defaults["task"]), type=str,
        help="The task name. See document of vise.")
    vasp_parser.add_argument(
        "--vise_opts", type=str, nargs="+", default=vasp_defaults["vise_opts"],
        help="Keyword arguments for options in make_input classmethod of "
             "ViseInputSet in vise. See document in vise for details.")
    vasp_parser.add_argument(
        "-uis", "--user_incar_settings", type=str, nargs="+",
        default=vasp_defaults["user_incar_settings"],
        help="user_incar_settings in make_input classmethod of ViseInputSet in "
             "vise. The default of this flag is set by the vise.yaml, so if one"
             "does not want to override the default, use "
             "additional_user_incar_setting instead. See also document in "
             "vise input_set for details.")
    vasp_parser.add_argument(
        "-auis", "--additional_user_incar_settings", type=str, nargs="+",
        default=None,
        help="Use this if one does not want to override user_incar_settings "
             "written in the yaml file")
    vasp_parser.add_argument(
        "--ldauu", type=str, default=vasp_defaults["ldauu"], nargs="+",
        help="Dict of LDAUU values")
    vasp_parser.add_argument(
        "--ldaul", type=str, default=vasp_defaults["ldaul"], nargs="+",
        help="Dict of LDAUL values.")
    vasp_parser.add_argument(
        "-c", "--charge", type=float, help="Charge state.")

    # -- get_poscars -----------------------------------------------------------
    parser_get_poscar = subparsers.add_parser(
        name="get_poscars",
        description="Tools for generating POSCAR file(s)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['gp'])

    gp_defaults = get_default_args(make_poscars_from_mp)

    parser_get_poscar.add_argument(
        "-p", "--poscar", type=str,
        help="POSCAR-type file name.", metavar="FILE")
    parser_get_poscar.add_argument(
        "-n", "--number", type=int,
        help="MP entry number without prefix 'mp-'")
    parser_get_poscar.add_argument(
        "-e", "--elements", type=str, nargs="+",
        help="Create directories with POSCARs containing the input elements.")
    parser_get_poscar.add_argument(
        "--e_above_hull", type=float, default=gp_defaults["e_above_hull"],
        help="Collect materials with this hull energy in eV/atom.")
    parser_get_poscar.add_argument(
        "--molecules", type=str2bool, default=gp_defaults["molecules"],
        help="Whether one doesn't want to replace pmg structures to molecules.")

    parser_get_poscar.set_defaults(func=get_poscar_from_mp)

    # -- vasp_set ---------------------------------------------------------
    parser_vasp_set = subparsers.add_parser(
        name="vasp_set",
        parents=[vasp_parser, prec_parser],
        description="Tools for constructing vasp input set with vise",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vs'])

    vs_defaults = {"kpt_density": KPT_DENSITY}
    simple_override(vs_defaults, ["kpt_density"])

    parser_vasp_set.add_argument(
        "--json", type=str,
        help="Json file name for printing the ViseInputSet info.")
    parser_vasp_set.add_argument(
        "-p", "--poscar", default="POSCAR", type=str,
        help="POSCAR-type file name.")
    parser_vasp_set.add_argument(
        "-k", "--kpt_density", default=vs_defaults["kpt_density"], type=float,
        help="K-point density in Angstrom along each direction .")
    parser_vasp_set.add_argument(
        "-s", "--standardize_structure", type=str2bool, default=True,
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive cell.")
    parser_vasp_set.add_argument(
        "-pi", "--prior_info", type=str2bool, default=True,
        help="Whether to use prior_info.json when it exists.")
    parser_vasp_set.add_argument(
        "--dirs", nargs="+", type=str, required=True,
        help="Make vasp set for the directories in the same condition.")
    parser_vasp_set.add_argument(
        "-d", "--prev_dir", type=str,
        help="Inherit input files from the previous directory.")

    del vs_defaults

    parser_vasp_set.set_defaults(func=vasp_set)

    # -- vasp_run --------------------------------------------------------------
    parser_vasp_run = subparsers.add_parser(
        name="vasp_run",
        parents=[vasp_parser, prec_parser],
        description="Tools for vasp run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vr'])

    vr_defaults = get_default_args(ViseVaspJob.kpt_converge)
    vr_defaults["vasp_cmd"] = None
    vr_defaults["timeout"] = TIMEOUT

    simple_override(vr_defaults, ["vasp_cmd",
                                  "initial_kpt_density"
                                  "timeout",
                                  "relax_iter_num",
                                  "convergence_criterion"])

    parser_vasp_run.add_argument(
        "-kc", "--kpoint_conv", action="store_true",
        help="Set if k-point convergence is checked.")
    parser_vasp_run.add_argument(
        "--json_file", default=None, type=str,
        help="Print structure_opt.json or kpt_conv.json (if -kc is set) info.")
    parser_vasp_run.add_argument(
        "-v", "--vasp_cmd", nargs="+", type=str,
        default=vr_defaults["vasp_cmd"],
        help="VASP command. If you are using mpirun, set this to something "
             "like \"mpirun pvasp\".",)
    parser_vasp_run.add_argument(
        "-ikd", "--initial_kpt_density", type=float,
        default=vr_defaults["initial_kpt_density"],
        help="Initial k-point density.")
    parser_vasp_run.add_argument(
        "--handler_name", type=str, default="default",
        help="Custodian error handler name listed in error_handlers.")
    parser_vasp_run.add_argument(
        "--timeout", type=int, default=vr_defaults["timeout"],
        help="Timeout used in TooLongTimeCalcErrorHandler.")
    parser_vasp_run.add_argument(
        "--remove_wavecar", action="store_true",
        help="Remove WAVECAR file after the calculation is finished.")
    parser_vasp_run.add_argument(
        "--max_relax_num", default=vr_defaults["max_relax_num"], type=int,
        help="Maximum number of relaxations.")
    parser_vasp_run.add_argument(
        "--criterion", dest="convergence_criterion",
        default=vr_defaults["convergence_criterion"], type=float,
        help="Convergence criterion of kpoints in eV / atom.")
    parser_vasp_run.add_argument(
        "--left_files", type=str, nargs="+", default=vr_defaults["left_files"],
        help="Filenames that are left at the calculation directory.")

    del vr_defaults
    parser_vasp_run.set_defaults(func=vasp_run)

    # -- chempotdiag -----------------------------------------------------------
    parser_cpd = subparsers.add_parser(
        name="chempotdiag",
        description="Tools for chemical potentials",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cpd'])

    # input
    # from file
    # parser_cpd.add_argument(
    #     "-e", "--energy", dest="energy_csv", type=str, default=None,
    #     help="Name of csv file of energies of compounds")
    parser_cpd.add_argument(
        "-dpd", "--draw_phase_diagram", action="store_true",
        help="set when drawing the phase diagram.")
    parser_cpd.add_argument(
        "-d", "--vasp_dirs", type=str, nargs='+',
        help="Drawing diagram from specified directories of vasp calculations")
    parser_cpd.add_argument(
        "-v", "--vasprun", default="vasprun.xml", type=str)
    parser_cpd.add_argument(
        "-e", "--elements", type=str, nargs="+",
        help="Element names. Obtain the total energies from the MP.")
    parser_cpd.add_argument(
        "-c", "--target_comp", type=str,
        help="Target compound focused for post process such as defect "
             "calculations.")
    parser_cpd.add_argument(
        "-f", dest="filename", type=str, help="Saved pdf file name.")
    parser_cpd.add_argument(
        "-pg", "--parse_gas", type=str2bool, default=True,
        help="Whether to parse results as gas for gas phases.")
    parser_cpd.add_argument(
        "-pp", "--partial_pressures", type=str, nargs='+',
        help="partial pressure of system in Pa. "
             "Example: -pp O2 1e+5 N2 20000 -> O2: 1e+5(Pa), N2: 20000(Pa)")
    parser_cpd.add_argument(
        "-t", "--temperature", type=float, default=0.0,
        help="temperature of system in K.")

    # thermodynamic status (P and T) input
    # # output
    # parser_cpd.add_argument("-y", "--yaml",
    #                     action="store_const", const=True,
    #                     default=False,
    #                     help="Dumps yaml of remarked_compound")

    parser_cpd.set_defaults(func=chempotdiag)

    # -- plot_band -----------------------------------------------------------
    parser_plot_band = subparsers.add_parser(
        name="plot_band",
        parents=[prec_parser],
        description="Tools for plotting band structures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pb'])

    parser_plot_band.add_argument(
        "-v", dest="vasprun", default="vasprun.xml", type=str)
    parser_plot_band.add_argument(
        "-v2", dest="vasprun2", type=str)
    parser_plot_band.add_argument(
        "-k", dest="kpoints", default="KPOINTS", type=str)
    parser_plot_band.add_argument(
        "-y", dest="y_range", nargs="+", type=float,
        help="Energy range, requiring two values.")
    parser_plot_band.add_argument(
        "-f", dest="filename", type=str, default=None, help="pdf file name.")
    parser_plot_band.add_argument(
        "-a", dest="absolute", action="store_true",
        help="Show in the absolute energy scale.")
    parser_plot_band.add_argument(
        "-l", dest="legend", action="store_false",
        help="Not show the legend.")

    parser_plot_band.set_defaults(func=plot_band)

    # -- plot_dos -----------------------------------------------------------
    parser_plot_dos = subparsers.add_parser(
        name="plot_dos",
        parents=[prec_parser],
        description="Tools for plotting density of states",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pd'])

    parser_plot_dos.add_argument(
        "-v", dest="vasprun", type=str, default="vasprun.xml")
    parser_plot_dos.add_argument(
        "-cv", dest="cbm_vbm", type=float, nargs="+",
        help="Set CBM and VBM.")
    parser_plot_dos.add_argument(
        "-t", dest="pdos_type", type=str, default="element",
        help=".")
    parser_plot_dos.add_argument(
        "-s", dest="specific", type=str, nargs="+", default=None, help=".")
    parser_plot_dos.add_argument(
        "-o", dest="orbital", type=str2bool, default=True,
        help="Whether to decompose orbital components.")
    parser_plot_dos.add_argument(
        "-x", "--x_range", nargs="+", type=float, default=None,
        help="Set energy minimum and maximum.")
    parser_plot_dos.add_argument(
        "-y", "--ymaxs", nargs="+", type=float, default=None,
        help="Set max values of y ranges. Support two ways."
             "1st: total_max, all_the_atoms" 
             "2nd: total_max, 1st_atom, 2nd_atom, ...")
    parser_plot_dos.add_argument(
        "-f", dest="filename", type=str, default="dos.pdf",
        help="Pdf file name.")
    parser_plot_dos.add_argument(
        "-a", "--absolute", action="store_true",
        help="Set when showing the figure in the absolute energy scale.")
    parser_plot_dos.add_argument(
        "-l", "--legend", type=str2bool, default=True,
        help="Whether to show the figure legend.")
    parser_plot_dos.add_argument(
        "-c", "--crop_first_value", type=str2bool, default=True,
        help="Whether to crop the first value in DOS.")

    parser_plot_dos.set_defaults(func=plot_dos)

    # -- band_gap --------------------------------------------------------------
    parser_band_gap = subparsers.add_parser(
        name="band_gap",
        description="Calculate the band gap from vasprun.xml",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['bg'])

    parser_band_gap.add_argument(
        "-v", "--vasprun", type=str, default="vasprun.xml", metavar="FILE")
    parser_band_gap.add_argument(
        "-o", "--outcar", type=str, default="OUTCAR", metavar="FILE")
    parser_band_gap.set_defaults(func=band_gap)

    # try:
    #     import argcomplete
    #     argcomplete.autocomplete(parser)
    #     # This supports bash autocompletion. To enable this, pip install
    #     # argcomplete, activate global completion, or add
    #     #      eval "$(register-python-argcomplete vise)"
    #     # into your .bash_profile or .bashrc
    # except ImportError:
    #     pass

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

