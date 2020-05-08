# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from pathlib import Path
from copy import deepcopy

from pymatgen.ext.matproj import MPRester
from pymatgen import Structure

from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.defaults import defaults
from vise.input_set.prior_info import PriorInfoFromCalcDir
from vise.input_set.kpoints_mode import KpointsMode
from vise.cli.main_tools import potcar_str2dict


def get_poscar_from_mp(args: Namespace) -> None:
    s = MPRester().get_structure_by_material_id(args.mpid)
    s.to(fmt="poscar", filename=args.poscar)


class VaspSet:
    def __init__(self, args: Namespace):
        self.args = args
        self._file_transfers = None

        options = CategorizedInputOptions(
            structure=self._structure(),
            task=args.task,
            xc=args.xc,
            kpt_density=args.kpt_density,
            overridden_potcar=self._overridden_potcar(),
            charge=args.charge,
            **self._option_kwargs())

        vif = VaspInputFiles(options, self._overridden_incar_settings())
        vif.create_input_files(Path.cwd())
        if self._file_transfers:
            self._file_transfers.transfer()

    def _structure(self):
        return Structure.from_file(self.args.poscar)

    def _overridden_incar_settings(self):
        result = deepcopy(defaults.user_incar_settings)
        if self.args.user_incar_settings:
            result.update(self.args.user_incar_settings)

        return result

    def _overridden_potcar(self):
        result = deepcopy(defaults.overridden_potcar)
        if self.args.overridden_potcar:
            result.update(potcar_str2dict(self.args.overridden_potcar))
        return result

    def _option_kwargs(self):
        result = deepcopy(defaults.options)
        if self.args.prev_dir:
            pi = PriorInfoFromCalcDir(self.args.prev_dir,
                                      self.args.file_transfer_type)
            result.update(pi.input_options_kwargs)
            self._file_transfers = pi.file_transfers

        if self.args.options:
            result.update(self.args.options)
        if self.args.uniform_kpt_mode:
            result["kpt_mode"] = KpointsMode.uniform

        return result
