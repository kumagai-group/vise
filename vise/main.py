#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from glob import glob
import os
from os.path import join
import warnings

from pymatgen.core.structure import Structure

from obadb.vasp.input_set import ObaSet
from obadb.util.structure_handler import SYMPREC

__version__ = "0.0.1"
__date__ = "8.Aug.2018"

# Following defaults determine the condition of automatic defect calculations.
# electronegativity difference for antisites and substitutional impurities
_EN_DIFF = 1.0
# Maximum displacement distance
_DISTANCE = 0.2
# Cutoff radius in which atoms are perturbed.
_CUTOFF = 3.0
_SYMPREC = SYMPREC
_ANGLE_TOLERANCE = -1.0


def main():
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package for first-principles point defect calculations. It 
    allows us to construct input files, parse first-principles calculation 
    results, and analyze data.""",
        epilog="""                                 
    Author: Yu Kumagai, Akira Takahashi
    Version: {}                                                                 
    Last updated: {}""".format(__version__, __date__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #        allow_abbrev=False)

    subparsers = parser.add_subparsers()

    # -- vasp_input_maker ----------------------------------------------------
    parser_vasp_input_maker = subparsers.add_parser(
        name="vasp_input_maker",
        description="Tools for configuring vasp INCAR, KPOINTS, POTCAR files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vm'])
    parser_vasp_input_maker.add_argument(
        "-d", "--dirname", dest="dirname", type=str)
    parser_vasp_input_maker.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_vasp_input_maker.add_argument(
        "--standardize", dest="standardize", action="store_false", default=True,
        help="")
    parser_vasp_input_maker.add_argument(
        "--xc", "-f", dest="xc", type=str, help="Functional.")
    parser_vasp_input_maker.add_argument(
        "--task", "-t", dest="task", type=str, help="Task.")
    parser_vasp_input_maker.add_argument(
        "--cluster", dest="cluster", action="store_true", default=False,
        help="")
    parser_vasp_input_maker.add_argument(
        "--rough", dest="rough", action="store_true", default=False,
        help="")
    parser_vasp_input_maker.add_argument(
        "--factor", dest="factor", type=int, default=1,
        help=".")
    parser_vasp_input_maker.add_argument(
        "--charge", dest="charge", type=float, default=1,
        help=".")
    parser_vasp_input_maker.add_argument(
        "--encut", dest="encut", type=float, default=None,
        help=".")
    parser_vasp_input_maker.add_argument(
        "--magnetization", dest="magnetization", action="store_true",
        default=False, help="")
    parser_vasp_input_maker.add_argument(
        "--kpts_density", dest="kpts_density", type=float, default=None,
        help=".")
    parser_vasp_input_maker.add_argument(
        "--only_even", dest="only_even", action="store_true", default=False,
        help="")
    parser_vasp_input_maker.add_argument(
        "--hubbard_u", dest="hubbard_u", action="store_true", default=False,
        help="")
    parser_vasp_input_maker.add_argument(
        "--ionic_contribution", dest="ionic_contribution", action="store_true",
        default=False,
        help="")
    parser_vasp_input_maker.add_argument(
        "--kpt_mode", dest="kpt_mode", type=str,
        help="")
    parser_vasp_input_maker.add_argument(
        "--auto_kpar_npar", dest="auto_kpar_npar", action="store_false",
        default=True,
        help="")
    parser_vasp_input_maker.add_argument(
        "--band_ref_dist", dest="band_ref_dist", type=float,
        default=0.03,
        help="K-point density used for the band structure calculation.")

    # POTCAR is generated from POSCAR only.
    parser_vasp_input_maker.set_defaults(func=vasp_input_maker)

    args = parser.parse_args()
    args.func(args)


# def vasp_input_maker(args):
#     structure = Structure.from_file(args.poscar)
#     set = ObaSet(structure=structure,
#                  standardize_structure=args.
#                  task=args.task,
#                  poscar=args.poscar,
#                  num_split_kpoints=args.num_split_kpoints,
#                  is_metal=args.is_metal,
#                  kpts_shift=args.kpts_shift,
#                  kpts_density_opt=args.kpts_density_opt,
#                  kpts_density_defect=args.kpts_density_defect,
#                  factor_dos=args.factor_dos,
#                  factor_metal=args.factor_metal,
#                  prior_info=args.prior_info,
#                  is_magnetization=args.is_magnetization)

    # elements = Structure.from_file(args.poscar).symbol_set
    # make_potcar(elements)

    # MakeIncar(task=args.task,
    #           xc=args.xc,
    #           poscar=args.poscar,
    #           potcar="POTCAR",
    #           hfscreen=args.hfscreen,
    #           aexx=args.aexx,
    #           is_metal=args.is_metal,
    #           is_magnetization=args.is_magnetization,
    #           ldau=args.ldau,
    #           defect_in=args.defect_in,
    #           prior_info=args.prior_info,
    #           my_incar_setting=args.my_setting_file)


if __name__ == "__main__":
    main()
