# -*- coding: utf-8 -*-

import itertools
import os
import re
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
incar_flags = loadfn(MODULE_DIR / "datasets" / "incar_flags.yaml")


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

    def get_string(self, raise_error: bool = False, **kwargs) -> str:
        """Override for the pretty printing. """
        lines = []
        check_incar_keys = list(self.keys())
        for key, val in incar_flags.items():
            comment = False
            blank_line = False
            ll = []

            for v in val.keys():
                if v in check_incar_keys:
                    if comment is False:
                        lines.append(f"# {key} \n")
                        comment = True
                    if v == "MAGMOM" and isinstance(self[v], list):
                        value = []
                        if type(self[v][0]) in list or Magmom \
                                and self.is_ncl_calc:
                            value.append(
                                " ".join(str(i) for j in self[v] for i in j))
                        elif self.get("LSORBIT") or self.get("LNONCOLLINEAR"):
                            for m, g in itertools.groupby(self[v]):
                                value.append("3*{}*{}".format(len(tuple(g)), m))
                        # YK: Add this
                        elif len(self[v]) == 1:
                            value.append(str(self[v][0]))
                        else:
                            # float() to ensure backwards compatibility between
                            # float magmoms and Magmom objects
                            try:
                                for m, g in itertools.groupby(self[v],
                                                              lambda x: float(
                                                                  x)):
                                    value.append(
                                        "{}*{}".format(len(tuple(g)), m))
                            except ValueError:
                                raise ValueError("MAGMOM could be improper.")

                        ll.append([v, " ".join(value)])
                    elif isinstance(self[v], list):
                        ll.append([v, " ".join([str(i) for i in self[v]])])
                    else:
                        ll.append([v, self[v]])
                    blank_line = True
                    check_incar_keys.remove(v)
            if blank_line:
                # if not disable_numparse, LOPTICS = True becomes 1.
                lines.append(str(tabulate([[l[0], "=", l[1]] for l in ll],
                                          tablefmt="plain",
                                          disable_numparse=True)) + "\n")
                lines.append("\n")  # blank space

        for mson_key in ["@module", "@class"]:
            try:
                check_incar_keys.remove(mson_key)
            except ValueError:
                pass

        if check_incar_keys:
            logger.error(f"{check_incar_keys} are invalid in INCAR.")
            if raise_error:
                raise ValueError

        return "".join(lines)

    @property
    def is_ncl_calc(self) -> bool:
        return any([self.get("LNONCOLLINEAR", False),
                    self.get("LSORBIT", False)])
