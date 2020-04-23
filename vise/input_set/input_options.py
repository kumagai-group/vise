# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from inspect import getfullargspec
from typing import Dict, Any, List
from copy import deepcopy

from pymatgen.core.structure import Structure

from vise.input_set.incar_settings_generator import IncarSettingsGenerator
from vise.input_set.potcar_generator import generate_potcar
from vise.input_set.structure_kpoints_generator import StructureKpointsGenerator
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)


class DividedInputOptions:
    def __init__(self,
                 input_structure: Structure,
                 task: Task,
                 xc: Xc,
                 **input_options):

        self._input_options = deepcopy(input_options)
        self._input_options.update(
            {"initial_structure": input_structure.copy(),
             "task": task,
             "xc": xc})
        self._raise_error_when_unknown_options_exist()

    def _raise_error_when_unknown_options_exist(self) -> None:
        unknown_args_set = self.input_option_set - self.assignable_option_set
        if unknown_args_set:
            raise ViseInputOptionsError(
                f"Options {unknown_args_set} are invalid")

    @property
    def input_option_set(self) -> set:
        return set(self._input_options.keys())

    @property
    def assignable_option_set(self) -> set:
        return (set(self.structure_kpoints_args)
                | set(self.potcar_args)
                | set(self.incar_settings_args))

    @property
    def structure_kpoints_args(self) -> List[str]:
        return self.constructor_args(StructureKpointsGenerator)

    @property
    def potcar_args(self) -> List[str]:
        return self.function_args(generate_potcar)

    @property
    def incar_settings_args(self) -> List[str]:
        return self.constructor_args(IncarSettingsGenerator)

    @staticmethod
    def function_args(func) -> List[str]:
        return getfullargspec(func).args

    def constructor_args(self, target_class) -> List[str]:
        result = self.function_args(target_class.__init__)
        result.remove("self")
        return result

    def pick_target_options(self, target_args: List[str]) -> Dict[str, Any]:
        result = {}
        for target_arg in target_args:
            if target_arg in self._input_options:
                result[target_arg] = self._input_options[target_arg]
        return result

    @property
    def structure_kpoints_options(self) -> Dict[str, Any]:
        return self.pick_target_options(target_args=self.structure_kpoints_args)

    @property
    def potcar_options(self) -> Dict[str, Any]:
        return self.pick_target_options(target_args=self.potcar_args)

    @property
    def incar_settings_options(self) -> Dict[str, Any]:
        return self.pick_target_options(target_args=self.incar_settings_args)

    @property
    def initial_structure(self) -> Structure:
        return self._input_options["initial_structure"]


class ViseInputOptionsError(KeyError):
    pass

