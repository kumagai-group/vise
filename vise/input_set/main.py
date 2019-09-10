#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from copy import deepcopy
from inspect import signature
from itertools import chain

from pydefect.util.logger import get_logger
from pydefect.util.main_tools import list2dict, potcar_str2dict
from pymatgen import Structure
from pymatgen.core.periodic_table import Element
from vise.util.config import SYMMETRY_TOLERANCE, ANGLE_TOL, KPT_DENSITY
from vise.input_set.incar import incar_flags
from vise.input_set.new_input_set import InputSet
from vise.input_set.prior_info import PriorInfo

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

__version__ = '0.0.1'
__date__ = 'will be inserted'


def main():

    parser = argparse.ArgumentParser(
        description="""                            
    vise is a package that helps researchers to do first-principles calculations 
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
    vos_defaults = {"vos_kwargs": {"symprec": SYMMETRY_TOLERANCE,
                                   "angle_tolerance": ANGLE_TOL},
                    "xc":         "pbesol",
                    "kpt_density": KPT_DENSITY,
                    "perfect_incar_setting": None,
                    "potcar_set": None,
                    "ldauu":      None,
                    "ldaul":      None}

    # write use potcar setting
    parser_vasp_set.add_argument(
        "-p", "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR-type file name.")
    parser_vasp_set.add_argument(
        "--potcar", dest="potcar_set", default=vos_defaults["potcar_set"],
        type=str, nargs="+",
        help="User specifying POTCAR set. E.g., Mg_pv O_h")
    parser_vasp_set.add_argument(
        "-x", "--xc", dest="xc", default=vos_defaults["xc"], type=str,
        help="Exchange-correlation (XC) interaction treatment.")
    parser_vasp_set.add_argument(
        "-t", "--task", dest="task", default="structure_opt", type=str,
        help="The task name. See document of vise.")
    parser_vasp_set.add_argument(
        "-k", "--kpt_density", dest="kpt_density",
        default=vos_defaults["kpt_density"], type=float,
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
        "-vos_kw", "--vos_kwargs", dest="vos_kwargs", type=str, nargs="+",
        default=vos_defaults["vos_kwargs"],
        help="Keyword arguments in make_input classmethod of ObaSet in vise. "
             "See document in vise for details.")
    parser_vasp_set.add_argument(
        "-is", "--incar_setting", dest="incar_setting", type=str, nargs="+",
        default=vos_defaults["perfect_incar_setting"],
        help="user_incar_setting in make_input classmethod of ObaSet in vise. "
             "See document in vise for details.")
    parser_vasp_set.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Make vasp set for the directories in the same condition."
             "In pydefect, it is used especially for competing phases.")
    parser_vasp_set.add_argument(
        "-pi", "-prior_info", dest="prior_info", action="store_true",
        help="Set if prior_info.json is read for competing phase calculations.")
    parser_vasp_set.add_argument(
        "-ldauu", dest="ldauu", type=dict, default=vos_defaults["ldauu"],
        nargs="+", help=".")
    parser_vasp_set.add_argument(
        "-ldaul", dest="ldaul", type=str, default=vos_defaults["ldaul"],
        nargs="+", help=".")

    # delete the default dict to avoid bugs.
    del vos_defaults

    parser_vasp_set.set_defaults(func=vasp_set)


def vasp_set(args):

    flags = [str(s) for s in list(Element)]
    ldauu = list2dict(args.ldauu, flags)
    ldaul = list2dict(args.ldaul, flags)
    potcar_set = potcar_str2dict(args.potcar_set)
    base_kwargs = {"task": args.task,
                   "xc": args.xc,
                   "kpt_density": args.kpt_density,
                   "standardize_structure": args.standardize,
                   "ldauu": ldauu,
                   "ldaul": ldaul}

    flags = list(chain.from_iterable(incar_flags.values()))
    base_user_incar_settings = list2dict(args.incar_setting, flags)

    flags = list(signature(InputSet.make_input).parameters.keys())
    base_kwargs.update(list2dict(args.vos_kwargs, flags))

    original_dir = os.getcwd()
    dirs = args.dirs if args.dirs else ["."]

    for d in dirs:
        os.chdir(os.path.join(original_dir, d))
        logger.info(f"Constructing vasp set in {d}")
        user_incar_settings = deepcopy(base_user_incar_settings)
        kwargs = deepcopy(base_kwargs)

        if args.prior_info:
            if os.path.exists("prior_info.json"):
                prior_info = PriorInfo.load_json("prior_info.json")
                kwargs["band_gap"] = prior_info["band_gap"]
                kwargs["is_magnetization"] = \
                    abs(prior_info["total_magnetization"]) > 0.1

        if args.prev_dir:
            files = {"CHGCAR": "C", "WAVECAR": "M", "WAVEDER": "M"}
            vasp_set = InputSet.from_prev_calc(args.prev_dir,
                                               charge=args.charge,
                                               files_to_transfer=files,
                                               **kwargs)
        else:
            s = Structure.from_file(args.poscar)
            vasp_set = \
                InputSet.make_input(structure=s,
                                    charge=args.charge,
                                    user_incar_settings=user_incar_settings,
                                    override_potcar_set=potcar_set,
                                    **kwargs)

        vasp_set.write_input(".")

    os.chdir(original_dir)
