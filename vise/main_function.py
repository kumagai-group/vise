# -*- coding: utf-8 -*-

import os
from copy import deepcopy
from itertools import chain

from custodian.custodian import Custodian

from pymatgen import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.outputs import BSVasprun

from vise.analyzer.band_gap import band_gap_properties
from vise.analyzer.band_plotter import PrettyBSPlotter
from vise.analyzer.dos_plotter import get_dos_plot
from vise.custodian_extension.jobs import ViseVaspJob
from vise.input_set.incar import incar_flags
from vise.input_set.input_set import ViseInputSet
from vise.input_set.prior_info import PriorInfo
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.error_classes import NoVaspCommandError
from vise.util.logger import get_logger
from vise.util.main_tools import potcar_str2dict, list2dict


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def vasp_set(args):
    if args.print_json:
        vis = ViseInputSet.load_json(args.print_json)
        print(vis)
        return

    flags = [str(s) for s in list(Element)]
    ldauu = list2dict(args.ldauu, flags)
    ldaul = list2dict(args.ldaul, flags)
    potcar_set = potcar_str2dict(args.potcar_set)
    task = Task.from_string(args.task)
    xc = Xc.from_string(args.xc)

    base_kwargs = {"kpt_density": args.kpt_density,
                   "standardize_structure": args.standardize,
                   "ldauu": ldauu,
                   "ldaul": ldaul}

    flags = list(ViseInputSet.ALL_OPTIONS.keys())
    base_kwargs.update(list2dict(args.vise_opts, flags))

    flags = list(chain.from_iterable(incar_flags.values()))
    base_user_incar_settings = list2dict(args.user_incar_setting, flags)

    original_dir = os.getcwd()
    dirs = args.dirs or ["."]

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
            input_set = ViseInputSet.from_prev_calc(args.prev_dir,
                                                    task=task,
                                                    xc=xc,
                                                    charge=args.charge,
                                                    files_to_transfer=files,
                                                    **kwargs)
        else:
            s = Structure.from_file(args.poscar)
            input_set = \
                ViseInputSet.make_input(structure=s,
                                        task=task,
                                        xc=xc,
                                        charge=args.charge,
                                        user_incar_settings=user_incar_settings,
                                        override_potcar_set=potcar_set,
                                        **kwargs)

        input_set.write_input(".")

    os.chdir(original_dir)


def vasp_run(args):
    if len(args.vasp_cmd) == 0:
        raise NoVaspCommandError
    elif len(args.vasp_cmd) == 1:
        vasp_cmd = args.vasp_cmd[0].split()
    else:
        vasp_cmd = args.vasp_cmd

    flags = list(chain.from_iterable(incar_flags.values()))
    user_incar_settings = list2dict(args.user_incar_setting, flags)
    print(user_incar_settings)
    # handlers
    handlers = []
    # handlers = [ObaVaspErrorHandler(), MeshSymmetryErrorHandler(),
    #             ObaUnconvergedErrorHandler(), NonConvergingErrorHandler(),
    #             PotimErrorHandler()]
    #    kpoints_criteria = args.kpoints_criteria

    # if args.converge_kpt:
    #     c = Custodian(handlers,
    #                   ViseVaspJob.kpt_converge(
    #                       cmd, max_relax_number, kpoints_criteria,
    #                   removes_wavecar=removes_wavecar),
    #                   polling_time_step=5, monitor_freq=1,
    #                   max_errors=10, gzipped_output=False)
    # else:
    # st = Poscar.from_file("POSCAR").structure
    # oba_vis = ObaSet.make_input(st)
    # # oba_vis = \
    # #     ObaSet.from_prev_calc(".",
    # #                           parse_calc_results=False,
    # #                           standardize_structure=True,
    # #                           parse_potcar=True,
    # #                           parse_incar=True,
    # #                           parse_kpoints=True)
    # oba_vis.write_input(".")
    # job = ViseVaspJob.structure_optimization_run(vasp_cmd=args.vasp_cmd,
    #                                              max_relax_num=args.max_relax_num,
    #                                              removes_wavecar=args.rm_wavecar)
    #    job =
    c = Custodian(handlers=handlers,
                  jobs=ViseVaspJob.kpt_converge(
                      vasp_cmd=vasp_cmd,
                      user_incar_settings=user_incar_settings,
                      left_files=args.left_files),
                  polling_time_step=5,
                  monitor_freq=1,
                  max_errors=10,
                  gzipped_output=False)
    c.run()


def plot_band(args):

    p = PrettyBSPlotter.from_vasp_files(kpoints=args.kpoints,
                                        vasprun=args.vasprun,
                                        vasprun2=args.vasprun2,
                                        absolute=args.absolute,
                                        y_range=args.y_range,
                                        legend=args.legend,
                                        symprec=args.symprec,
                                        angle_tolerance=args.angle_tolerance)

    p.show(args.filename, format_type="pdf")


def plot_dos(args):
    dos = get_dos_plot(vasprun_file=args.vasprun,
                       cbm_vbm=args.cbm_vbm,
                       pdos_type=args.pdos_type,
                       specific=args.specific,
                       orbital=args.orbital,
                       xlim=args.x_range,
                       ymaxs=args.ymaxs,
                       zero_at_efermi=args.absolute,
                       legend=args.legend,
                       crop_first_value=args.cfv,
                       symprec=args.symprec,
                       angle_tolerance=args.angle_tolerance)

    if args.filename:
        dos.savefig(args.filename, format="pdf")
    else:
        dos.show()


def band_gap(args):
    v = BSVasprun(args.vasprun)
    try:
        band_gap_value, vbm_info, cbm_info = band_gap_properties(v)
        print(f"CBM info {cbm_info}")
        print(f"VBM info {vbm_info}")
        print(f"band gap info {band_gap_value}")
    except TypeError:
        print("Metallic system")