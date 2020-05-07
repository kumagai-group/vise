# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
from argparse import Namespace
from copy import deepcopy
from itertools import chain
from pathlib import Path
from typing import Tuple

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Vasprun, Outcar


from vise.analyzer.band_edge_properties import band_gap_properties
from vise.analyzer.band_plotter import PrettyBSPlotter
from vise.analyzer.dos_plotter import get_dos_plot
from vise.input_set.incar import all_incar_flags
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.input_set.prior_info import PriorInfo
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.cli.main_tools import potcar_str2dict, list2dict

logger = get_logger(__name__)


def vasp_settings_from_args(args: Namespace
                            ) -> Tuple[dict, dict]:
    """Generate vasp input settings from the given args. """

    flags = [str(s) for s in list(Element)]
    ldauu = list2dict(args.ldauu, flags)
    ldaul = list2dict(args.ldaul, flags)
    potcar_set = potcar_str2dict(args.potcar_set)
    key_candidates = list(VaspInputFiles.ALL_OPTIONS.keys())
    # Sanitize values of vis_base_kwargs and user_incar_settings with list2dict.
    vis_base_kwargs = list2dict(args.vise_opts, key_candidates)
    key_candidates = list(chain.from_iterable(all_incar_flags.values()))
    user_incar_settings = list2dict(args.user_incar_settings, key_candidates)

    if args.additional_user_incar_settings:
        d = list2dict(args.additional_user_incar_settings, key_candidates)
        user_incar_settings.update(d)

    vis_base_kwargs.update({"potcar_set_name": args.potcar_set_name,
                            "override_potcar_set": potcar_set,
                            "ldauu": ldauu,
                            "ldaul": ldaul})
    if hasattr(args, "charge"):
        vis_base_kwargs["charge"] = args.charge
    return user_incar_settings, vis_base_kwargs


def get_poscar_from_mp(args: Namespace) -> None:
    s = MPRester().get_structure_by_material_id(args.number)
    s.to(fmt="poscar", filename=args.poscar)


def vasp_set(args: Namespace) -> None:
    if args.print:
        vis = VaspInputFiles.load_json(args.json)
        print(vis)
        return

    # user_incar_settings and base_vis_kwargs can be modified by prior_info.json
    base_user_incar_settings, base_vis_kwargs = vasp_settings_from_args(args)
    task = Task.from_string(args.task)
    xc = Xc.from_string(args.xc)
    base_vis_kwargs.update(
        {"symprec": args.symprec,
         "angle_tolerance": args.angle_tolerance,
         "kpt_density": args.kpt_density,
         "standardize_structure": args.standardize_structure})

    started_dir = os.getcwd()
    cwd = Path.cwd()

    for d in args.dirs:
        d = cwd / d
        os.chdir(d)
        logger.info(f"Constructing vasp set in {d}")
        user_incar_settings = deepcopy(base_user_incar_settings)
        kwargs = deepcopy(base_vis_kwargs)

        if args.prior_info:
            yaml_file = d / "prior_info.yaml"
            if yaml_file.exists():
                prior_info = PriorInfo.load_yaml(str(yaml_file))
                if prior_info.is_magnetic:
                    kwargs["is_magnetization"] = True
                if prior_info.is_cluster:
                    if task != Task.cluster_opt:
                        logger.warning(f"task is changed from {task} to "
                                       f"{Task.cluster_opt}.")
                        task = Task.cluster_opt
                if prior_info.incar:
                    user_incar_settings.update(prior_info.incar)

        if args.prev_dir:
            logger.warning("user_incar_settings is not used when using "
                           "previous calculation.")
            files = {"CHGCAR": "C", "WAVECAR": "M", "WAVEDER": "M"}
            input_set = VaspInputFiles.from_prev_calc(args.prev_dir,
                                                      task=task,
                                                      xc=xc,
                                                      files_to_transfer_dict=files,
                                                      **kwargs)
        else:
            s = Structure.from_file(args.poscar)
            input_set = \
                VaspInputFiles.make_input(structure=s,
                                          task=task,
                                          xc=xc,
                                          user_incar_settings=user_incar_settings,
                                          **kwargs)

        input_set.create_input(".")

    os.chdir(started_dir)


# def plot_band(args) -> None:
#     p = PrettyBSPlotter.from_vasp_files(kpoints=args.kpoints,
#                                         vasprun=args.vasprun,
#                                         vasprun2=args.vasprun2,
#                                         absolute=args.absolute,
#                                         y_range=args.y_range,
#                                         legend=args.legend,
#                                         symprec=args.symprec,
#                                         angle_tolerance=args.angle_tolerance)
#     p.bsp.plot_brillouin()
# #    p.show(args.filename, format_type="pdf")


# def plot_dos(args) -> None:
#     if args.cbm_vbm:
#         if len(args.cbm_vbm) != 2 or args.cbm_vbm[0] < args.cbm_vbm[1]:
#             raise ValueError(f"cbm_vbm values {args.cbm_vbm} are not proper.")

    # dos = get_dos_plot(vasprun=args.vasprun,
    #                    cbm_vbm=args.cbm_vbm,
    #                    pdos_type=args.pdos_type,
    #                    specific=args.specific,
    #                    orbital=args.orbital,
    #                    xlim=args.x_range,
    #                    ymaxs=args.ymaxs,
    #                    zero_at_efermi=not args.absolute,
    #                    legend=args.legend,
    #                    crop_first_value=args.crop_first_value,
    #                    symprec=args.symprec,
    #                    angle_tolerance=args.angle_tolerance)
    # dos.savefig(args.filename, format="pdf")


# def band_gap(args) -> None:
#     try:
#         band_gap_value, vbm_info, cbm_info = \
#             band_gap_properties(vasprun=Vasprun(args.vasprun),
#                                 outcar=Outcar(args.outcar))
#         print(f"CBM info {cbm_info}")
#         print(f"VBM info {vbm_info}")
#         print(f"band gap info {band_gap_value}")
#     except TypeError:
#         print("Metallic system")


