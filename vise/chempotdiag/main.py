#!/usr/bin/env python

#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from inspect import signature
import os
import warnings

from pymatgen.core.periodic_table import Element
from chempotdiag.chem_pot_diag import ChemPotDiag
from chempotdiag.make_inputs import make_vasp_inputs_from_mp

from atomate.utils.utils import get_logger

__author__ = "Akira Takahashi"
__maintainer__ = "Akira Takahashi"

logger = get_logger(__name__)


def overwrite_default_args(class_method, main_args):
    """ Use the defaults in class.classmethod

    Args:
        class_method (classmethod): classmethod. When using __init__, class is fine.
        main_args (dict): Args set by main

    Return:
        args (dict): Overwritten args by options
    """
    from inspect import _empty

    args_with_default = []
    sig = signature(class_method)

    for name, param in sig.parameters.items():
        if param.default != _empty:
            args_with_default.append(name)

    args = {}
    for a in args_with_default:
        if hasattr(main_args, a):
            if getattr(main_args, a) is not None:
                args[a] = getattr(main_args, a)

    return args


def main():
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package for first-principles point defect calculations. It 
    allows us to construct input files, parse first-principles calculation 
    results, and analyze data.""",
        epilog="""                                 
    Author: Akira Takahashi
    Version: {}                                                                 
    Last updated: {}""".format("0.0.1", "will be inserted"),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # get poscar from materials project
    parser.add_argument("-m", "--mat_proj_poscar",
                        help="",
                        action="store_true")

    parser.add_argument("-el", "--elements",
                        dest="elements", type=str, nargs='+',
                        default=None,
                        help="")

    parser.add_argument("-dp", "--dir_path",
                        dest="dir_path", type=str,
                        default=None,
                        help="")

    parser.add_argument("-ch", "--criterion_hull",
                        dest="criterion_hull", type=float,
                        default=None,
                        help="Collect materials only if energy above hull is "
                             "less than criterion_hull. Unit is meV/atom.")

    parser.add_argument("-k", "--mp_api_key",
                        help="",
                        action="store_true")

    parser.add_argument("-gp", "--gets_poly",
                        help="",
                        action="store_true")

    parser.add_argument("-no_mol", "--not_molecules",
                        help="",
                        action="store_true")

    # input

    # from file
    parser.add_argument("-e", "--energy", dest="energy_file",
                        type=str, default=None,
                        help="Name of text file of "
                             "energies of compounds")

    # from VASP
    parser.add_argument("-v", "--vasp_dirs",
                        dest="vasp_dirs", type=str, nargs='+',
                        default=None,
                        help="Drawing diagram from specified "
                             "directories of vasp calculations")
    # from VASP and MP
    parser.add_argument("-fmp_target", "--from_mp_target",
                        help="VASP result of target material,"
                             "when get competing phases from mp")

    parser.add_argument("-fmp_elem", "--from_mp_element", nargs="*",
                        help="VASP result of elements,"
                             "when get competing phases from mp")

    # VASP_option
    parser.add_argument("-p", "--poscar_name",
                        dest="poscar_name", type=str,
                        default="POSCAR",
                        help="Name of POSCAR, like CONTCAR, "
                             "POSCAR-finish,...")
    parser.add_argument("-o", "--outcar_name",
                        dest="outcar_name", type=str,
                        default="OUTCAR",
                        help="Name of OUTCAR, like OUTCAR-finish")

    parser.add_argument("-es", "--energy_shift", type=str,
                        dest="energy_shift",
                        nargs='+', default=None,
                        help="Energy shift, "
                             "e.g. -es N2/molecule 1 "
                             "-> make more unstable N2/molecule "
                             "by 1 eV")

    # thermodynamic status (P and T) input
    parser.add_argument("-pp", "--partial_pressures",
                        dest="partial_pressures", type=str,
                        nargs='+', default=None,
                        help="partial pressure of system. "
                             "e.g. -pp O2 1e+5 N2 20000 "
                             "-> O2: 1e+5(Pa), N2: 20000(Pa)")

    parser.add_argument("-t", "--temperature",
                        dest="temperature", type=float,
                        default=0,
                        help="temperature of system (unit: K)"
                             "e.g. -t 3000 -> 3000(K)")

    # drawing diagram
    parser.add_argument("-w", "--without_label",
                        help="Draw diagram without label.",
                        action="store_true")

    parser.add_argument("-c", "--remarked_compound",
                        dest="remarked_compound", type=str,
                        default=None,
                        help="Name of compound you are remarking."
                             "Outputted equilibrium_points are "
                             "limited to neighboring that "
                             "compounds, and those equilibrium "
                             "points are labeled in "
                             "chem_pot_diagram.")

    parser.add_argument("-d", "--draw_range",
                        dest="draw_range", type=float,
                        default=None,
                        help="Drawing range of diagram."
                             "If range is shallower than the "
                             "deepest vertex, "
                             "ValueError will occur")

    # output
    parser.add_argument("-s", "--save_file",
                        dest="save_file", type=str,
                        default=None,
                        help="File name to save the drawn diagram.")

    parser.add_argument("-y", "--yaml",
                        action="store_const", const=True,
                        default=False,
                        help="Dumps yaml of remarked_compound")

    parser.set_defaults(func=chempotdiag)

    args = parser.parse_args()
    args.func(args)


def chempotdiag(args):
    if args.mat_proj_poscar:
        kwargs_to_make_vasp_inputs = {}
        if args.dir_path:
            kwargs_to_make_vasp_inputs["dir_path"] = args.dir_path
        if args.criterion_hull is not None:
            kwargs_to_make_vasp_inputs["criterion_e_above_hull"] \
                = args.criterion_hull
        if args.mp_api_key:
            kwargs_to_make_vasp_inputs["api_key"] = args.mp_api_key
        if args.gets_poly:
            kwargs_to_make_vasp_inputs["gets_poly"] = True
        if args.not_molecules:
            kwargs_to_make_vasp_inputs["adds_molecule"] = False
        make_vasp_inputs_from_mp(args.elements, **kwargs_to_make_vasp_inputs)
    else:
        if len([f for f in
                [args.energy_file, args.vasp_dirs, args.from_mp_target]
                if bool(f)]) != 1:
            raise ValueError("Specify one of energy_file, vasp_dirs "
                             "and from_mp_target")

        # energy_shift
        energy_shift_dict = {}
        if args.energy_shift:
            if len(args.energy_shift) % 2 != 0:
                raise ValueError(f"Invalid energy shift "
                                 f"input {args.energy_shift}")
            for i in range(int(len(args.energy_shift) / 2)):
                output_name = \
                    args.energy_shift[2 * i] + "/" + args.outcar_name
                es = args.energy_shift[2 * i + 1]
                energy_shift_dict[output_name] = float(es)

        # pressure and temperature
        partial_pressure_dict = {}
        if args.partial_pressures:
            if len(args.partial_pressures) % 2 != 0:
                raise ValueError(f"Invalid partial pressures "
                                 f"input {args.partial_pressures}")
            for i in range(int(len(args.partial_pressures) / 2)):
                formula = args.partial_pressures[2 * i]
                pressure = args.partial_pressures[2 * i + 1]
                partial_pressure_dict[formula] = float(pressure)

        if args.energy_file:
            if args.temperature or args.partial_pressures:
                warnings.warn("Now temperature and pressures can not apply when"
                              " reading data from energy_file")
            cp = ChemPotDiag.from_file(args.energy_file)

        elif args.vasp_dirs:

            poscar_paths = [d + args.poscar_name for d in args.vasp_dirs]
            outcar_paths = [d + args.outcar_name for d in args.vasp_dirs]

            cp = ChemPotDiag. \
                from_vasp_calculations_files(poscar_paths,
                                             outcar_paths,
                                             temperature=args.temperature,
                                             pressure=partial_pressure_dict,
                                             energy_shift_dict=energy_shift_dict)
            if args.elements:
                cp.set_elements([Element(e) for e in args.elements])

        elif args.from_mp_target:
            elem_poscar_paths = [os.path.join(d, args.poscar_name) for d in args.from_mp_element]
            elem_outcar_paths = [os.path.join(d, args.outcar_name) for d in args.from_mp_element]
            cp = ChemPotDiag.from_vasp_and_materials_project(
                vasp_target_poscar=f"{args.from_mp_target}/{args.poscar_name}",
                vasp_target_output=f"{args.from_mp_target}/{args.outcar_name}",
                vasp_element_poscar=elem_poscar_paths,
                vasp_element_output=elem_outcar_paths,
                temperature=args.temperature,
                pressure=partial_pressure_dict,
                energy_shift_dict=energy_shift_dict
            )
            if args.elements:
                cp.set_elements([Element(e) for e in args.elements])

        print(f"Energies of elements ({cp.elements}) : {cp.element_energy}")
        #  Read args of drawing diagram from parser
        if args.remarked_compound:
            try:
                for vertex in cp.get_neighbor_vertices(args.remarked_compound):
                    print(vertex)
            except ValueError:
                print(f"{args.remarked_compound} is unstable."
                      f" No vertex is labeled.")

        kwargs_for_diagram = {}
        if args.remarked_compound:
            kwargs_for_diagram["remarked_compound"] = args.remarked_compound
        if args.save_file:
            kwargs_for_diagram["save_file_name"] = args.save_file
        if args.without_label:
            kwargs_for_diagram["with_label"] = False
        if args.draw_range:
            kwargs_for_diagram["draw_range"] = args.draw_range

        if cp.dim >= 4:
            print("Currently diagram is not available for quaternary or more.")
        else:
            cp.draw_diagram(**kwargs_for_diagram)
            # try:
            #     cp.draw_diagram(**kwargs_for_diagram)
            # except ValueError:
            #     kwargs_for_diagram.pop("remarked_compound")
            #     cp.draw_diagram(**kwargs_for_diagram)

        if args.yaml:
            if args.remarked_compound is None:
                raise ValueError("remarked_compound is needed to dump yaml")
            cp.dump_vertices_yaml(os.getcwd(), args.remarked_compound, args.elements)


if __name__ == "__main__":
    main()
