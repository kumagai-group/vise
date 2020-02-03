#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from typing import Union

from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL, KPT_DENSITY
from vise.main_function import vasp_set, plot_band, plot_dos, vasp_run, band_gap
from vise.util.logger import get_logger
from vise.util.main_tools import dict2list, get_user_settings

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

__version__ = '0.0.1dev'
__date__ = 'will be inserted'

# The following keys are set by vise.yaml
setting_keys = ["vasp_command",
                "symprec",
                "angle_tolerance",
                "xc",
                "kpt_density",
                "user_incar_setting",
                "ldauu",
                "ldaul",
                "potcar_set",
                "relax_iter_num",
                "kpoints_criteria"]

user_settings = get_user_settings(yaml_filename="vise.yaml",
                                  setting_keys=setting_keys)


def main():

    def simple_override(d: dict, keys: Union[list, str]) -> None:
        """Override dict if keys exist in vise.yaml.

        When the value in the user_settings is a dict, it will be changed to
        list using dict2list.
        """
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            if key in user_settings:
                v = user_settings[key]
                if isinstance(v, dict):
                    v = dict2list(v)
                d[key] = v

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

    # -- vasp_set ---------------------------------------------------------
    parser_vasp_set = subparsers.add_parser(
        name="vasp_set",
        description="Tools for constructing vasp input set with vise",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vs'])

    # all the defaults must be declared here.
    vs_defaults = {"symprec":            SYMMETRY_TOLERANCE,
                   "angle_tolerance":    ANGLE_TOL,
                   "xc":                 "pbesol",
                   "task":               "structure_opt",
                   "kpt_density":        KPT_DENSITY,
                   "vise_opts":          None,
                   "user_incar_setting": None,
                   "potcar_set":         None,
                   "ldauu":              None,
                   "ldaul":              None}

    simple_override(vs_defaults, list(vs_defaults.keys()))

    # write use potcar setting
    parser_vasp_set.add_argument(
        "--pj", dest="print_json", type=str,
        help="Print the ViseInputSet info from the given json file.")
    parser_vasp_set.add_argument(
        "-p", "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR-type file name.")
    parser_vasp_set.add_argument(
        "--potcar", dest="potcar_set",
        default=vs_defaults["potcar_set"],
        type=str, nargs="+",
        help="User specifying POTCAR set. E.g., Mg_pv O_h")
    parser_vasp_set.add_argument(
        "-x", "--xc", dest="xc", default=vs_defaults["xc"], type=str,
        help="Exchange-correlation (XC) interaction treatment.")
    parser_vasp_set.add_argument(
        "-t", "--task", dest="task", default=vs_defaults["task"], type=str,
        help="The task name. See document of vise.")
    parser_vasp_set.add_argument(
        "-k", "--kpt_density", dest="kpt_density",
        default=vs_defaults["kpt_density"], type=float,
        help="K-point density in Angstrom along each direction .")
    parser_vasp_set.add_argument(
        "-s", "--standardize", dest="standardize", action="store_false",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive cell.")
    parser_vasp_set.add_argument(
        "-d", "--prev_dir", dest="prev_dir", type=str,
        help="Inherit input files from the previous directory.")
    parser_vasp_set.add_argument(
        "-c", "--charge", dest="charge", type=int, default=0,
        help="Supercell charge state.")
    parser_vasp_set.add_argument(
        "-vise_opts", dest="vise_opts", type=str, nargs="+",
        default=vs_defaults["vise_opts"],
        help="Keyword arguments for options in make_input classmethod of "
             "ViseInputSet in vise. See document in vise for details.")
    # TODO: the vise.yaml is completely overridden. Should we partly override?
    parser_vasp_set.add_argument(
        "-uis", "--user_incar_setting", dest="user_incar_setting", type=str,
        nargs="+",
        default=vs_defaults["user_incar_setting"],
        help="user_incar_setting in make_input classmethod of ViseInputSet in "
             "vise. See document in vise for details.")
    parser_vasp_set.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str, default=None,
        help="Make vasp set for the directories in the same condition.")
    parser_vasp_set.add_argument(
        "-pi", "-prior_info", dest="prior_info", action="store_true",
        help="Set if prior_info.json is read.")
    parser_vasp_set.add_argument(
        "-ldauu", dest="ldauu", type=dict, default=vs_defaults["ldauu"],
        nargs="+", help="Dict of LDAUU values")
    parser_vasp_set.add_argument(
        "-ldaul", dest="ldaul", type=str, default=vs_defaults["ldaul"],
        nargs="+", help="Dict of LDAUL values.")
    parser_vasp_set.add_argument(
        "--symprec", dest="symprec", type=float,
        default=vs_defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_vasp_set.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=vs_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")

    del vs_defaults

    parser_vasp_set.set_defaults(func=vasp_set)

    # -- vasp_run --------------------------------------------------------------
    parser_vasp_run = subparsers.add_parser(
        name="vasp_run",
        description="Tools for vasp run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vr'])

    vr_defaults = {"symprec":             SYMMETRY_TOLERANCE,
                   "angle_tolerance":     ANGLE_TOL,
                   "left_files":          None,
                   "vasp_command":        None,
                   "user_incar_setting":  None,
                   "relax_iter_num":      10,
                   "kpoints_criteria":    0.003}

    simple_override(vr_defaults, list(vr_defaults.keys()))

    parser_vasp_run.add_argument(
        "-vc", "--vasp_command", dest="vasp_cmd", nargs="+", type=str,
        default=None,
        help="VASP command. If you are using mpirun, set this to something "
             "like \"mpirun pvasp\".",)
    parser_vasp_run.add_argument(
        "-handler_name", dest="handler_name", type=str, default="default",
        help="Custodian error handler name listed in error_handlers.")
    parser_vasp_run.add_argument(
        "-timeout", dest="timeout", type=int, default=None,
        help="Timeout used in TooLongTimeCalcErrorHandler.")
    parser_vasp_run.add_argument(
        "-rw", "--remove_wavecar", dest="rm_wavecar", action="store_true",
        help="Remove WAVECAR file after the calculation is finished.")
    parser_vasp_run.add_argument(
        "--max_relax_num", dest="max_relax_num",
        default=vr_defaults["relax_iter_num"], type=int,
        help="Maximum number of relaxations.")
    parser_vasp_run.add_argument(
        "-uis", "--user_incar_setting", dest="user_incar_setting", type=str,
        nargs="+",
        default=vr_defaults["user_incar_setting"],
        help="user_incar_setting in make_input classmethod of ViseInputSet in "
             "vise. See document in vise for details.")
    parser_vasp_run.add_argument(
        "-criteria", dest="kpoints_criteria",
        default=vr_defaults["kpoints_criteria"], type=float,
        help="Convergence criteria of kpoints in eV / atom.")
    parser_vasp_run.add_argument(
        "-files", dest="left_files", type=str, nargs="+",
        default=vr_defaults["left_files"],
        help="Filenames that are left at the calculation directory.")
    parser_vasp_run.add_argument(
        "-kc", "-kpoint_conv", dest="kpoint_conv", action="store_true",
        help="Set if k-point convergence is checked.")
    parser_vasp_run.add_argument(
        "-x", "--xc", dest="xc", default=None, type=str,
        help="Exchange-correlation (XC) interaction treatment.")

#    parser_vasp_run.add_argument(
#        "-kd", "-kpoint_density", dest="kpoint_density", type=float,
#        default=None, help="Initial k-point density.")

    del vr_defaults
    parser_vasp_run.set_defaults(func=vasp_run)

    # -- plot_band -----------------------------------------------------------
    parser_plot_band = subparsers.add_parser(
        name="plot_band",
        description="Tools for plotting band structures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pb'])

    pb_defaults = {"symprec":         SYMMETRY_TOLERANCE,
                   "angle_tolerance": ANGLE_TOL}

    simple_override(pb_defaults, list(pb_defaults.keys()))

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
    parser_plot_band.add_argument(
        "--symprec", dest="symprec", type=float,
        default=pb_defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_plot_band.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=pb_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")

    del pb_defaults
    parser_plot_band.set_defaults(func=plot_band)

    # -- plot_dos -----------------------------------------------------------
    parser_plot_dos = subparsers.add_parser(
        name="plot_dos",
        description="Tools for plotting density of states",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pd'])

    pd_defaults = {"symprec":         SYMMETRY_TOLERANCE,
                   "angle_tolerance": ANGLE_TOL}

    simple_override(pd_defaults, list(pd_defaults.keys()))

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
        "-o", dest="orbital", action="store_false",
        help="Switch off the orbital decomposition.")
    parser_plot_dos.add_argument(
        "-x", dest="x_range", nargs="+", type=float, default=None,
        help="Set energy minimum and maximum.")
    parser_plot_dos.add_argument(
        "-y", dest="ymaxs", nargs="+", type=float, default=None,
        help="Set max values of y ranges. Support two ways."
             "1st: total_max, all_the_atoms" 
             "2nd: total_max, 1st_atom, 2nd_atom, ...")
    parser_plot_dos.add_argument(
        "-f", dest="filename", type=str, default="dos.pdf",
        help="pdf file name.")
    parser_plot_dos.add_argument(
        "-a", dest="absolute", action="store_true",
        help="Show in the absolute energy scale.")
    parser_plot_dos.add_argument(
        "-l", dest="legend", action="store_false",
        help="Not show the legend.")
    parser_plot_dos.add_argument(
        "-c", dest="c", action="store_false",
        help="Not crop the first value.")
    parser_plot_dos.add_argument(
        "--symprec", dest="symprec", type=float,
        default=pd_defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_plot_dos.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=pd_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")

    del pd_defaults

    parser_plot_dos.set_defaults(func=plot_dos)

    # -- band_gap --------------------------------------------------------------
    parser_band_gap = subparsers.add_parser(
        name="band_gap",
        description="Calculate the band gap from vasprun.xml",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['bg'])

    vr_defaults = {"symprec":          SYMMETRY_TOLERANCE,
                   "angle_tolerance":  ANGLE_TOL,
                   "vasp_command":     None,
                   "relax_iter_num":   10,
                   "kpoints_criteria": 0.03}

    simple_override(vr_defaults, list(vr_defaults.keys()))

    parser_band_gap.add_argument(
        "-v", dest="vasprun", type=str, default="vasprun.xml", metavar="FILE")
    parser_band_gap.add_argument(
        "-o", dest="outcar", type=str, default="OUTCAR", metavar="FILE")
    parser_band_gap.set_defaults(func=band_gap)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

