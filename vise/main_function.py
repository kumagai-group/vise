# -*- coding: utf-8 -*-

from copy import deepcopy
from itertools import chain
import os

from custodian.custodian import Custodian
from pathlib import Path

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import BSVasprun, Outcar

from vise.analyzer.band_gap import band_gap_properties
from vise.analyzer.band_plotter import PrettyBSPlotter
from vise.analyzer.dos_plotter import get_dos_plot
from vise.chempotdiag.chem_pot_diag import ChemPotDiag
from vise.chempotdiag.free_energy_entries import FreeEnergyEntrySet
from vise.custodian_extension.handler_groups import handler_group
from vise.custodian_extension.jobs import (
    ViseVaspJob, KptConvResult, StructureOptResult)
from vise.input_set.incar import incar_flags
from vise.input_set.input_set import ViseInputSet
from vise.input_set.prior_info import PriorInfo
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.error_classes import NoVaspCommandError
from vise.util.logger import get_logger
from vise.util.main_tools import potcar_str2dict, list2dict
from vise.util.mp_tools import make_poscars_from_mp


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_poscar_from_mp(args):
    if args.number:
        s = MPRester().get_structure_by_material_id(f"mp-{args.number}")
        s.to(filename=args.poscar)
        return True

    if args.elements:
        make_poscars_from_mp(elements=args.elements,
                             e_above_hull=args.e_above_hull,
                             molecules=args.molecules)
        return True

    logger.warning("Set mp number or elements")
    return False


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

    vise_base_kwargs = {"kpt_density": args.kpt_density,
                        "standardize_structure": args.standardize,
                        "potcar_set_name": args.potcar_set_name,
                        "ldauu": ldauu,
                        "ldaul": ldaul}

    key_candidates = list(ViseInputSet.ALL_OPTIONS.keys())
    vise_base_kwargs.update(list2dict(args.vise_opts, key_candidates))

    key_candidates = list(chain.from_iterable(incar_flags.values()))
    base_user_incar_settings = list2dict(args.user_incar_setting,
                                         key_candidates)
    if args.additional_user_incar_setting:
        d = list2dict(args.additional_user_incar_setting, key_candidates)
        base_user_incar_settings.update(d)

    original_dir = os.getcwd()
    dirs = args.dirs or ["."]

    for d in dirs:
        d = Path.cwd() / d
        os.chdir(d)
        logger.info(f"Constructing vasp set in {d}")
        user_incar_settings = deepcopy(base_user_incar_settings)
        kwargs = deepcopy(vise_base_kwargs)

        if args.prior_info:
            yaml_file = d / "prior_info.yaml"
            if yaml_file.exists():
                print(yaml_file)
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

    if args.print:
        if args.kpoint_conv:
            filename = args.json_file or "kpt_conv.json"
            kpt_conv = KptConvResult.load_json(filename)
            print(kpt_conv)
        else:
            filename = args.json_file or "structure_opt.json"
            str_opt = StructureOptResult.load_json(filename)
            print(str_opt)
        return

    # TODO: implement a way to use different type of vasp (vasp_ncl, e.g.)
    if isinstance(args.vasp_cmd, str):
        vasp_cmd = args.vasp_cmd.split()
    elif isinstance(args.vasp_cmd, list):
        if len(args.vasp_cmd) == 1:
            vasp_cmd = args.vasp_cmd[0].split()
        else:
            vasp_cmd = args.vasp_cmd
    else:
        raise NoVaspCommandError("Vasp command must be specified properly.")

    flags = list(chain.from_iterable(incar_flags.values()))
    user_incar_settings = list2dict(args.user_incar_setting, flags)
    if args.additional_user_incar_setting:
        key_candidates = list(chain.from_iterable(incar_flags.values()))
        d = list2dict(args.additional_user_incar_setting, key_candidates)
        user_incar_settings.update(d)

    handlers = handler_group(args.handler_name, timeout=args.timeout)

    optimization_args = {"vasp_cmd": vasp_cmd,
                         "removes_wavecar": args.rm_wavecar,
                         "max_relax_num": args.max_relax_num,
                         "left_files": args.left_files,
                         "removed_files": ["PCDAT", "vasprun.xml"],
                         "symprec": args.symprec,
                         "angle_tolerance": args.angle_tolerance}

    custodian_args = {"handlers": handlers,
                      "polling_time_step": 5,
                      "monitor_freq": 1,
                      "max_errors": 10,
                      "gzipped_output": False}

    if args.kpoint_conv:
        custodian_args["jobs"] = ViseVaspJob.kpt_converge(
            xc=args.xc,
            convergence_criterion=args.convergence_criterion,
            initial_kpt_density=args.initial_kpt_density,
            user_incar_settings=user_incar_settings,
            **optimization_args)
    else:
        custodian_args["jobs"] = ViseVaspJob.structure_optimization_run(
            **optimization_args)

    c = Custodian(**custodian_args)
    c.run()


def chempotdiag(args):
    if args.elements:
        entry_set = FreeEnergyEntrySet.from_mp(args.elements)
    else:
        entry_set = FreeEnergyEntrySet.from_vasp_files(
            directory_paths=args.dirs,
            vasprun=args.vasprun,
            parse_gas=args.parse_gas,
            temperature=args.temperature,
            partial_pressures=args.partial_pressures)

    entry_set.to_json()
    pd = PhaseDiagram(entries=entry_set.entries)
    if args.draw_phase_diagram:
        pd_plotter = PDPlotter(pd)
        if args.pd_filename:
            image_format = args.pd_filename.split(".")[-1]
            pd_plotter.write_image(args.pd_filename, image_format)
        else:
            pd_plotter.show()

    cpd = ChemPotDiag.from_phase_diagram(pd, target_comp=args.target_comp)
    cpd.to_json_file()

    cpd.draw_diagram(filename=args.filename)


def plot_band(args):
    p = PrettyBSPlotter.from_vasp_files(kpoints_filenames=args.kpoints,
                                        vasprun_filenames=args.vasprun,
                                        vasprun2_filenames=args.vasprun2,
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
                       zero_at_efermi=not args.absolute,
                       legend=args.legend,
                       crop_first_value=args.c,
                       symprec=args.symprec,
                       angle_tolerance=args.angle_tolerance)

    dos.savefig(args.filename, format="pdf")


def band_gap(args):
    v = BSVasprun(args.vasprun)
    o = Outcar(args.outcar)
    try:
        band_gap_value, vbm_info, cbm_info = band_gap_properties(v, o)
        print(f"CBM info {cbm_info}")
        print(f"VBM info {vbm_info}")
        print(f"band gap info {band_gap_value}")
    except TypeError:
        print("Metallic system")
