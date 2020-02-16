# -*- coding: utf-8 -*-

from copy import deepcopy
from itertools import chain
import os

from custodian.custodian import Custodian

from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import BSVasprun, Outcar

from vise.analyzer.band_gap import band_gap_properties
from vise.analyzer.band_plotter import PrettyBSPlotter
from vise.analyzer.dos_plotter import get_dos_plot
from vise.chempotdiag.chem_pot_diag import ChemPotDiag
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


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_poscar_from_mp(args):
    s = MPRester().get_structure_by_material_id(f"mp-{args.number}")
    s.to(filename=args.poscar)


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
        os.chdir(os.path.join(original_dir, d))
        logger.info(f"Constructing vasp set in {d}")
        user_incar_settings = deepcopy(base_user_incar_settings)
        kwargs = deepcopy(vise_base_kwargs)

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
                logger.warning("Now temperature and pressures can not apply "
                               "when reading data from energy_file")
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
