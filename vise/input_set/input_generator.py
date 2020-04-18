# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
from typing import Optional, Dict, Any

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import Poscar, Kpoints, Potcar, VaspInput

from vise import __version__
from vise.input_set.incar import ViseIncar
from vise.input_set.structure_kpoints_generator import StructureKpointsGenerator
from vise.input_set.potcar_generator import generate_potcar
from vise.input_set.incar_settings_generator import IncarSettingsGenerator
from vise.util.logger import get_logger
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.structure_handler import create_symbol_list

logger = get_logger(__name__)


class ViseInputSet:
    def __init__(self,
                 input_structure: Structure,
                 task: Task,
                 xc: Xc,
                 sort_structure: bool = True,
                 input_options: dict = None,
                 user_incar_settings: Optional[dict] = None):

        self._version = __version__
        self._initial_structure = input_structure.copy()

        if sort_structure:
            input_structure = input_structure.get_sorted_structure()
        else:
            input_structure = input_structure.copy()

        skj = StructureKpointsGenerator(input_structure, task=task)
        skj.generate_input()
        self._structure = skj.structure
        self._kpoints = skj.kpoints
        symbol_list = create_symbol_list(skj.structure)
        self._potcar = generate_potcar(symbol_list=symbol_list, xc=xc)
        incar_settings_generator = \
            IncarSettingsGenerator(skj.structure.composition,
                                   symbol_list,
                                   skj.num_kpts,
                                   skj.num_kpt_factor,
                                   self._potcar,
                                   task,
                                   xc)
        self._incar_settings = incar_settings_generator.incar_settings
        self._incar_settings.update(user_incar_settings or {})

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
        return {"INCAR": self.incar,
                "KPOINTS": self.kpoints,
                "POSCAR": self.poscar,
                "POTCAR": self.potcar}

    def create_input_files(self,
                           dirname: Path,
                           poscar_significant_figures=10) -> None:
        dirname.mkdir(exist_ok=True)
        for filename, obj in self.input_files.items():
            if isinstance(obj, Poscar):
                obj.write_file(dirname / filename,
                               significant_figures=poscar_significant_figures)
            else:
                obj.write_file(dirname / filename)

