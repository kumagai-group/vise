# -*- coding: utf-8 -*-

import itertools
import os
import re
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path

from monty.io import zopen
from monty.serialization import loadfn

from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp import Incar
from pymatgen.util.io_utils import clean_lines

from vise.util.logger import get_logger

from tabulate import tabulate


MODULE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
# This incar_flags should be OrderedDict, but from python 3.6, dict uses
# order-preserving semantics. Besides, it does not affect vasp result.
incar_tags = loadfn(MODULE_DIR / "datasets" / "incar_flags.yaml")


logger = get_logger(__name__)


class ViseIncar(Incar):
    """Incar class modified for pretty writing of INCAR file.

    """
    @classmethod
    def from_file(cls, filename: str) -> "ViseIncar":
        """Reads an Incar object from a file.

        Since from_file and from_string methods in Incar class use Incar class
        constructor, we need to override them.
        """
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())

    @classmethod
    def from_dict(cls, d: dict) -> "ViseIncar":
        kwargs = deepcopy(d)
        if kwargs.get("MAGMOM") and isinstance(kwargs["MAGMOM"][0], dict):
            kwargs["MAGMOM"] = [Magmom.from_dict(m) for m in kwargs["MAGMOM"]]

        return cls({k: v for k, v in d.items()
                    if k not in ("@module", "@class")})

    @classmethod
    def from_string(cls, string: str):
        """Reads an Incar object from a string.

        This method is different from that of Incar superclass as it supports
        ";" semantic, which splits the incar flags.

        Args:
            string (str): Incar string

        Returns:
            ViseIncar object
        """
        params = {}
        for line in string.splitlines():
            # Support ";" semantic for split of tags in INCAR file
            for sub_line in line.split(";"):
                m = re.match(r"\s*(\w+)\s*=\s*(.*)\s*", sub_line)
                if m:
                    key = m.group(1).strip()
                    val = m.group(2).strip()
                    params[key] = ViseIncar.proc_val(key, val)
        return cls(params)

    def __add__(self, other: Incar) -> "ViseIncar":
        """Add all the values of another INCAR object to this object. """
        params = {k: v for k, v in self.items()}
        for k, v in other.items():
            if k in self and v != self[k]:
                raise ValueError("INCARs have conflicting values!")
            else:
                params[k] = v
        return ViseIncar(params)

    def get_string(self, **kwargs) -> str:
        """Override for the pretty printing. """
        lines = OrderedDict()
        check_incar_keys = list(self.keys())

        for category, flags in incar_tags.items():
            tag_list_by_category = []

            for incar_tag in flags.keys():
                if incar_tag in check_incar_keys:
                    tag_list_by_category.append(
                        [incar_tag, "=", self.str_value(incar_tag)])
                    check_incar_keys.remove(incar_tag)

            if tag_list_by_category:
                # if not disable_numparse, e.g., LOPTICS = True becomes 1.
                tabulated_str = str(tabulate(tag_list_by_category,
                                             tablefmt="plain",
                                             disable_numparse=True))
                lines[f"# {category}"] = tabulated_str

        if check_incar_keys:
            logger.error(f"{check_incar_keys} may be invalid in INCAR.")

        return "\n\n".join(["\n".join([k, v]) for k, v in lines.items()])

    def str_value(self, incar_tag: str) -> str:
        if isinstance(self[incar_tag], list):
            return " ".join([str(i) for i in self[incar_tag]])
        else:
            return self[incar_tag]

    @property
    def is_ncl_calc(self) -> bool:
        return any([self.get("LNONCOLLINEAR", False),
                    self.get("LSORBIT", False)])
