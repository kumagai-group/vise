# -*- coding: utf-8 -*-

import argparse
import os

from pymatgen.core.structure import Structure
#from obadb.vasp.kpoints import make_kpoints
from obadb.vasp.input_set import ObaSet
from obadb.vasp.kpoints import make_kpoints
from obadb.util.structure_handler import find_hpkot_primitive


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p", "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR-type file name for the supercell.")
    parser.add_argument(
        "-x", "--xc", dest="xc", default="pbesol", type=str,
        help="XC interaction treatment.")
    parser.add_argument(
        "-t", "--task", dest="task", default="structure_opt", type=str,
        help="The task name.")
    parser.add_argument(
        "-k", "--kpt_density", dest="kpt_density", default=3, type=float,
        help="K-points density.")
    parser.add_argument(
        "-s", "--standardize", dest="standardize", action="store_false",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive one.")
    parser.add_argument(
        "-d", "--prev_dir", dest="prev_dir", type=str,
        help=".")
    parser.add_argument(
        "-gga", "--prev_gga", dest="prev_gga",
        action="store_true", help=".")
    parser.add_argument(
        "-gpre", "--prev_dir_gw_pre2", dest="prev_dir_gw_pre2", type=str,
        help=".")
    parser.add_argument(
        "-gw", "--prev_dir_gw0", dest="prev_dir_gw0", type=str,
        help=".")
    # parser.add_argument(
    #     "--defect", dest="defect", type=, help=".")
    parser.add_argument(
        "-vr", "-very_rough", dest="very_rough", action="store_true",
        help=".")
    parser.add_argument(
        "-is", "--incar_setting", dest="incar_setting", type=str, nargs="+",
        help="INCAR setting to be overwritten.")

    args = parser.parse_args()

    kwargs = {"task": args.task,
              "xc": args.xc,
              "kpt_density": args.kpt_density,
              "standardize_structure": args.standardize}

    if args.prev_dir_gw_pre2:
        obrs = ObaSet.from_prev_calc_gw_pre2(args.prev_dir_gw_pre2)
    elif args.prev_dir_gw0:
        obrs = ObaSet.from_prev_calc_gw(args.prev_dir_gw0)
    elif args.prev_dir:
        files = {"CHGCAR": "C", "WAVECAR": "L", "WAVEDER": "M"}
        obrs = ObaSet.from_prev_calc(args.prev_dir, copied_file_names=files)
            # , parse_calc_results=False,
            #                    parse_magnetization=False, parse_potcar=False,
            #                    parse_incar=False, parse_kpoints=False, **kwargs)
    elif args.prev_gga:
        s = Structure.from_file(args.poscar)
        user_incar_settings = {"LWAVE": True}
        default_potcar = "default_GW_POTCAR_list"

        obrs = ObaSet.make_input(s,
                                 standardize_structure=args.standardize,
                                 task="single_point",
                                 kpt_density=args.kpt_density,
                                 default_potcar=default_potcar,
                                 user_incar_settings=user_incar_settings,
                                 rough=args.rough)

    else:
        user_incar_settings = {}
        if args.incar_setting:
            if len(args.incar_setting) % 2 == 1:
                raise ValueError("INCAR settings must be even number.")
            for i, j in zip(args.incar_setting[0::2], args.incar_setting[1::2]):
                try:
                    j = float(j)
                except ValueError:
                    pass
                user_incar_settings[i] = j

        s = Structure.from_file(args.poscar)
        obrs = ObaSet.make_input(s,
                                 task=args.task,
                                 xc=args.xc,
                                 kpt_density=args.kpt_density,
                                 user_incar_settings=user_incar_settings,
                                 rough=args.very_rough)

    obrs.write_input(".")

    # structure = Structure.from_file(args.poscar)
    # obrs = ObaSet.make_input(structure, **kwargs)

