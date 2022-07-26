# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy
from pathlib import Path
from typing import Optional, Dict, Any

from pymatgen.io.vasp.sets import Poscar, Kpoints, Potcar

from vise import __version__
from vise.input_set.incar import ViseIncar
from vise.input_set.incar_settings_generator import IncarSettingsGenerator
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.kpoints import ViseKpoints
from vise.input_set.potcar_generator import generate_potcar
from vise.input_set.structure_kpoints_generator import StructureKpointsGenerator
from vise.input_set.vise_log import ViseLog
from vise.util.logger import get_logger
from vise.util.structure_handler import create_symbol_list

logger = get_logger(__name__)


class VaspInputFiles:
    def __init__(self,
                 input_options: CategorizedInputOptions,
                 overridden_incar_settings: Optional[dict] = None):

        self._version = __version__
        self._input_options = input_options
        self._initial_structure = input_options.initial_structure
        self._generate_structure_kpoints(input_options)
        self._generate_potcar_incar_settings(input_options)
        self._overridden_incar_settings = overridden_incar_settings
        self._incar_settings.update(overridden_incar_settings or {})

    def _generate_structure_kpoints(self, input_options):
        str_kpoints_generator = StructureKpointsGenerator(
            **input_options.structure_kpoints_options)
        str_kpoints_generator.generate_input()

        self._structure = str_kpoints_generator.structure
        self._kpoints = str_kpoints_generator.kpoints
        self._num_kpts = str_kpoints_generator.num_kpts
        self._num_kpt_factor = str_kpoints_generator.num_kpt_factor

    def _generate_potcar_incar_settings(self, input_options):
        symbol_list = create_symbol_list(self._structure)
        self._potcar = generate_potcar(symbol_list=symbol_list,
                                       **input_options.potcar_options)
        incar_settings_generator = \
            IncarSettingsGenerator(self._structure,
                                   symbol_list,
                                   self._num_kpts,
                                   self._num_kpt_factor,
                                   self._potcar,
                                   **input_options.incar_settings_options)
        self._incar_settings = incar_settings_generator.incar_settings

    @property
    def vise_log(self):
        result = copy(self._input_options.vise_log)
        result.user_incar_settings = self._overridden_incar_settings
        return result

    @property
    def initial_structure(self):
        return self._initial_structure

    @property
    def version(self) -> str:
        return self._version

    @property
    def incar(self) -> ViseIncar:
        incar = ViseIncar.from_dict(self._incar_settings)
        return incar

    @property
    def poscar(self) -> Poscar:
        return Poscar(self._structure)

    @property
    def kpoints(self) -> Kpoints:
        return self._kpoints

    @property
    def potcar(self) -> Potcar:
        return self._potcar

    @property
    def input_files(self) -> Dict[str, Any]:
        return {"INCAR": self.incar, "KPOINTS": self.kpoints,
                "POSCAR": self.poscar, "POTCAR": self.potcar}

    def create_input_files(self,
                           dirname: Path,
                           poscar_significant_figures=10) -> None:
        dirname.mkdir(exist_ok=True)
        for filename, obj in self.input_files.items():
            if isinstance(obj, Poscar):
                # Need to show more digit than pymatgen one.
                obj.write_file(dirname / filename,
                               significant_figures=poscar_significant_figures)
            elif isinstance(obj, Kpoints):
                kpoints = ViseKpoints.from_dict(obj.as_dict())
                kpoints.write_file(dirname / filename)
            else:
                obj.write_file(dirname / filename)

