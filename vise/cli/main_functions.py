# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import sys
from argparse import Namespace
from copy import deepcopy
from pathlib import Path

import yaml
from pymatgen.core import Structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import DieleFuncData
from vise.analyzer.plot_band_dos import BandDosPlotlyPlotter
from vise.analyzer.plot_diele_func_data import DieleFuncMplPlotter, \
    DieleFuncPlotType
from vise.analyzer.plot_band import BandMplPlotter
from vise.analyzer.plot_dos import DosPlotter
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.analyzer.vasp.dos_data import DosDataFromVasp
from vise.analyzer.vasp.make_diele_func import make_diele_func
from vise.analyzer.vasp.make_effective_mass import make_effective_mass
from vise.analyzer.vasp.make_irreps import special_points_from_kpoints, \
    make_irreps_from_wavecar
from vise.analyzer.vasp.plot_band import BandPlotInfoFromVasp
from vise.atom_energies.make_atom_vasp_set import make_atom_poscar_dirs
from vise.cli.main_tools import potcar_str2dict, list2dict
from vise.defaults import defaults
from vise.input_set.datasets.dataset_util import all_incar_flags
from vise.input_set.input_options import CategorizedInputOptions, \
    assignable_option_set
from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.prior_info import prior_info_from_calc_dir, PriorInfo
from vise.input_set.task import Task
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.util.file_transfer import FileTransfers
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


def structure_info(args: Namespace) -> None:
    input_structure = Structure.from_file(args.poscar)
    symmetrizer = StructureSymmetrizer(input_structure,
                                       symprec=args.symprec,
                                       angle_tolerance=args.angle_tolerance)
    if args.show_primitive:
        primitive = symmetrizer.primitive
        if primitive != input_structure:
            print(primitive.to(fmt="poscar"))
        else:
            logger.info("Input structure is a primitive cell.")
    elif args.show_conventional:
        conventional = symmetrizer.conventional
        if conventional != input_structure:
            print(conventional.to(fmt="poscar"))
        else:
            logger.info("Input structure is a conventional cell.")
    else:
        print(symmetrizer)


def get_poscar_from_mp(args: Namespace) -> None:
    s = MPRester().get_structure_by_material_id(args.mpid)
    s.to(fmt="poscar", filename=args.poscar)
    data = MPRester().get_data(args.mpid)[0]
    d = {"total_magnetization": data["total_magnetization"],
         "band_gap": data["band_gap"],
         "data_source": args.mpid,
         "icsd_ids": data["icsd_ids"]}
    args.prior_info.write_text(yaml.dump(d), None)


class VaspSet:
    def __init__(self, args: Namespace):
        self.args = args

        try:
            self._prior_info = PriorInfo.load_yaml()
        except FileNotFoundError:
            self._prior_info = PriorInfo()
        self.task = Task.cluster_opt if self._prior_info.is_cluster \
            else args.task

        options = CategorizedInputOptions(
            structure=self._structure(),
            task=self.task,
            xc=args.xc,
            kpt_density=args.kpt_density,
            overridden_potcar=self._overridden_potcar(),
            **self._option_kwargs())

        vif = VaspInputFiles(options, self._overridden_incar_settings())
        vif.create_input_files(Path.cwd())
        if hasattr(self, "_file_transfers"):
            self._file_transfers.transfer()
        vif.vise_log.to_yaml_file()

    def _structure(self):
        # avoid overlapping structure for e.g., phonon forces.
        if self.args.prev_dir and self.task.fine_to_inherit_structure_from_prev:
            return Structure.from_file(self.args.prev_dir / defaults.contcar)
        return Structure.from_file(self.args.poscar)

    def _overridden_incar_settings(self):
        result = deepcopy(defaults.user_incar_settings)
        result.update(self._prior_info.incar)

        if self.args.user_incar_settings:
            args = list2dict(self.args.user_incar_settings, all_incar_flags)
            result.update(args)

        return result

    def _overridden_potcar(self):
        result = {}
        if self.args.overridden_potcar:
            result.update(potcar_str2dict(self.args.overridden_potcar))
        return result

    def _option_kwargs(self):
        result = deepcopy(defaults.options)
        if self._prior_info:
            result.update(self._prior_info.input_options_kwargs)

        if self.args.prev_dir:
            pi = prior_info_from_calc_dir(prev_dir_path=self.args.prev_dir,
                                          vasprun=self.args.vasprun,
                                          outcar=self.args.outcar)
            result.update(pi.input_options_kwargs)

            self._file_transfers = FileTransfers(self._file_transfer(),
                                                 path=self.args.prev_dir)

        if self.args.options:
            args = list2dict(self.args.options, assignable_option_set)
            if "kpt_density" in args and self.args.kpt_density:
                logger.warning(
                    f"Because kpt_density is set at kpt_density argument and "
                    f"optional argument, kpt_density argument value "
                    f"{self.args.kpt_density} is used.")
                args.pop("kpt_density")
            result.update(args)
        if self.args.uniform_kpt_mode:
            result["kpt_mode"] = KpointsMode.uniform

        return result

    def _file_transfer(self):
        if not self.args.file_transfer_type:
            return {}
        else:
            result = {}
            for filename, transfer in zip(self.args.file_transfer_type[0::2],
                                          self.args.file_transfer_type[1::2]):
                result[filename] = transfer
            return result


def plot_band(args: Namespace):
    band_plot_info_from_vasp = BandPlotInfoFromVasp(
        vasprun=Vasprun(args.vasprun), kpoints_filename=args.kpoints_filename)
    plot_info = band_plot_info_from_vasp.make_band_plot_info()

    if args.wavecar_filename:
        try:
            with open(args.kpoints_filename) as f:
                first_line = f.readline()
            sg_num = int(first_line.split()[-1])
            logger.info(f"Space group number is set to {sg_num} read from "
                        f"{args.kpoints_filename} header.")
        except ValueError:
            logger.info(f"Failed to get space group number from "
                        f"{args.kpoints_filename} header.")
            sg_num = None

        special_points, kpt_indices = \
            special_points_from_kpoints(args.kpoints_filename)
        irreps = make_irreps_from_wavecar(
            special_points, kpt_indices, sg_num, args.wavecar_filename,
            args.poscar)
        plot_info.band_energy_infos["1"].irreps = irreps

    plot_info.to_json_file()
    if args.plotly:
        plotter = BandDosPlotlyPlotter(band_plot_info=plot_info)

        plotter.fig.show()
        return
    plotter = BandMplPlotter(plot_info, energy_range=args.y_range)
    plotter.construct_plot()
    plotter.plt.savefig(args.filename, format="pdf")


def plot_dos(args: Namespace):
    vasprun = Vasprun(args.vasprun)
    outcar = Outcar(args.outcar)
    band_edge = VaspBandEdgeProperties(vasprun, outcar)

    if band_edge.band_gap:
        vertical_lines = [band_edge.vbm_info.energy, band_edge.cbm_info.energy]
    else:
        vertical_lines = [vasprun.efermi]

    if args.base_energy is None:
        base = vertical_lines[0]
    else:
        base = args.base_energy

    dos_data_from_vasp = \
        DosDataFromVasp(vasprun, vertical_lines, base, args.crop_first_value)
    dos_data = dos_data_from_vasp.make_dos_data()

    ylim_set = None
    if args.y_max_ranges:
        if dos_data.spin:
            ylim_set = [[-y_max, y_max] for y_max in args.y_max_ranges]
        else:
            ylim_set = [[0, y_max] for y_max in args.y_max_ranges]

    structure = vasprun.final_structure
    grouped_atom_indices = args.type.grouped_atom_indices(structure,
                                                          args.target)
    logger.info(f"Grouped atom indices: {grouped_atom_indices}")
    plot_data = dos_data.dos_plot_data(grouped_atom_indices,
                                       energy_range=args.x_range,
                                       dos_ranges=ylim_set)
    plot_data.to_json_file()
    plotter = DosPlotter(plot_data, args.legend)
    plotter.construct_plot()
    plotter.plt.savefig(args.filename, format="pdf")


def plot_diele_func(args: Namespace):
    if args.input_csv_name:
        diele_func_data = DieleFuncData.from_csv_file(args.input_csv_name)
    else:
        diele_func_data = make_diele_func(Vasprun(args.vasprun),
                                          Outcar(args.outcar),
                                          use_vasp_real=not args.calc_kk,
                                          ita=args.ita)
    diele_func_data.to_json_file()
    if args.to_csv:
        diele_func_data.to_csv_file()
    if args.y_range:
        if args.plot_type is DieleFuncPlotType.absorption_coeff:
            y_range = [10 ** args.y_range[0], 10 ** args.y_range[1]]
        else:
            y_range = args.y_range
    else:
        y_range = None

    plotter = DieleFuncMplPlotter(diele_func_data)
    plotter.construct_plot(directions=args.directions,
                           plot_type=args.plot_type, y_range=y_range)
    filename = args.filename or str(args.plot_type) + ".pdf"
    plotter.plt.savefig(filename, format="pdf")


def calc_effective_mass(args: Namespace):
    vasprun, outcar = Vasprun(args.vasprun), Outcar(args.outcar)
    band_edge_prop = VaspBandEdgeProperties(vasprun, outcar)
    try:
        vbm, cbm = band_edge_prop.vbm_cbm
    except TypeError:
        logger.warning("Band gap does not exist, so not suited for effective"
                       "mass calculation.")
        return
    effective_mass = make_effective_mass(vasprun,
                                         args.temperature,
                                         args.concentrations,
                                         vbm, cbm)
    print(effective_mass)


def band_edge_properties(args: Namespace):
    vasprun = Vasprun(args.vasprun)
    outcar = Outcar(args.outcar)
    print(VaspBandEdgeProperties(vasprun, outcar))
