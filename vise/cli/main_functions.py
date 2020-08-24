# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from copy import deepcopy
from pathlib import Path

import yaml
from pymatgen import Structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.plot_band import BandPlotter
from vise.analyzer.plot_dos import DosPlotter
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.analyzer.vasp.dos_data import DosDataFromVasp
from vise.analyzer.vasp.plot_band import BandPlotInfoFromVasp
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


def get_poscar_from_mp(args: Namespace) -> None:
    s = MPRester().get_structure_by_material_id(args.mpid)
    s.to(fmt="poscar", filename=args.poscar)
    data = MPRester().get_data(args.mpid)[0]
    d = {"total_magnetization": data["total_magnetization"],
         "band_gap": data["band_gap"],
         "data_source": args.mpid}
    args.prior_info.write_text(yaml.dump(d), None)


class VaspSet:
    def __init__(self, args: Namespace):
        self.args = args
        self._file_transfers = None

        try:
            self._prior_info = PriorInfo.load_yaml()
        except FileNotFoundError:
            self._prior_info = PriorInfo()
        task = Task.cluster_opt if self._prior_info.is_cluster else args.task

        options = CategorizedInputOptions(
            structure=self._structure(),
            task=task,
            xc=args.xc,
            kpt_density=args.kpt_density,
            overridden_potcar=self._overridden_potcar(),
            **self._option_kwargs())

        vif = VaspInputFiles(options, self._overridden_incar_settings())
        vif.create_input_files(Path.cwd())
        if self._file_transfers:
            self._file_transfers.transfer()

    def _structure(self):
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
    plot_info.to_json_file()
    plotter = BandPlotter(plot_info, energy_range=args.y_range)
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

    dos_data_from_vasp = DosDataFromVasp(vasprun, vertical_lines, base,
                                         args.crop_first_value)
    dos_data = dos_data_from_vasp.make_dos_data()

    ylim_set = None
    if args.y_max_ranges:
        if dos_data.spin:
            ylim_set = [[-y_max, y_max] for y_max in args.y_max_ranges]
        else:
            ylim_set = [[0, y_max] for y_max in args.y_max_ranges]

    structure = vasprun.final_structure
    grouped_atom_indices = args.type.grouped_atom_indices(structure, args.target)
    plot_data = dos_data.dos_plot_data(grouped_atom_indices,
                                       xlim=args.x_range,
                                       ylim_set=ylim_set)
    plot_data.to_json_file()
    plotter = DosPlotter(plot_data, args.legend)
    plotter.construct_plot()
    plotter.plt.savefig(args.filename, format="pdf")


def band_edge_properties(args: Namespace):
    vasprun = Vasprun(args.vasprun)
    outcar = Outcar(args.outcar)
    print(VaspBandEdgeProperties(vasprun, outcar))
