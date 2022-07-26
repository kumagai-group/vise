# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
import re
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path
from typing import List, Dict, Any

from monty.io import zopen
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp import Incar
from tabulate import tabulate

from vise.input_set.datasets.dataset_util import incar_categories
from vise.util.logger import get_logger

MODULE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))


logger = get_logger(__name__)


class ViseIncar(Incar):
    """Incar class modified for pretty writing of INCAR file. """
    @classmethod
    def from_file(cls, filename: str) -> "ViseIncar":
        """
        Since from_file and from_string methods in Incar class use 'Incar' class
        constructor, we need to override them.
        """
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ViseIncar":
        kwargs = deepcopy(d)
        if kwargs.get("MAGMOM") and isinstance(kwargs["MAGMOM"][0], dict):
            kwargs["MAGMOM"] = [Magmom.from_dict(m) for m in kwargs["MAGMOM"]]

        return cls({k: v for k, v in d.items()
                    if k not in ("@module", "@class")})

    @classmethod
    def from_string(cls, string: str) -> "ViseIncar":
        """Reads an Incar object from a string.

        This method is different from that of Incar superclass as it supports
        ";" semantic, which splits the incar flags.
        """
        params = {}
        for line in string.splitlines():
            # Support ";" semantic for split of flags in INCAR file
            for sub_line in line.split(";"):
                m = re.match(r"\s*(\w+)\s*=\s*(.*)\s*", sub_line)
                if m:
                    key = m.group(1).strip()
                    val = m.group(2).strip()
                    params[key] = ViseIncar.proc_val(key, val)
        return cls(params)

    def __add__(self, other: Incar) -> "ViseIncar":
        """Add all the values of another INCAR object to this object. """
        original_settings = {k: v for k, v in self.items()}
        for key, value in other.items():
            if key in self and value != self[key]:
                raise ValueError("INCARs have conflicting values!")
            else:
                original_settings[key] = value
        return ViseIncar(original_settings)

    def get_string(self, **kwargs) -> str:
        lines = OrderedDict()
        input_incar_flags = list(self.keys())

        for category, flags_by_category in incar_categories.items():
            settings_by_category = []

            for flag in flags_by_category:
                if flag in input_incar_flags:
                    settings_by_category.append(
                        [flag, self.setting_to_str(flag)])
                    input_incar_flags.remove(flag)

            if settings_by_category:
                lines[f"# {category}"] = tabulated_string(settings_by_category)

        if input_incar_flags:
            logger.error(f"{input_incar_flags} may be invalid in INCAR.")
            lines[f"# others"] = tabulated_string(input_incar_flags)

        return "\n\n".join(["\n".join([k, v]) for k, v in lines.items()])

    def setting_to_str(self, incar_flag: str) -> str:
        if isinstance(self[incar_flag], list):
            return " ".join([str(i) for i in self[incar_flag]])
        else:
            return str(self[incar_flag])

    @property
    def is_ncl_calc(self) -> bool:
        return any([self.get("LNONCOLLINEAR", False),
                    self.get("LSORBIT", False)])


def tabulated_string(target_list: List[list]) -> str:
    """Tabulate pairs of keys and values with sandwiching "=".

    Examples: target_list = [["NSW", 1], "ISMEAR", -5]]
     ->
NSWW   =  1
ISMER  = -5"""
    # if not disable_numparse, e.g., LOPTICS = True becomes 1.
    equal_added = [[key, "=", value] for (key, value) in target_list]
    return str(tabulate(equal_added, tablefmt="plain", disable_numparse=True))

