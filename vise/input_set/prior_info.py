# -*- coding: utf-8 -*-

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen.io.vasp import Vasprun, Outcar

from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.defaults import defaults
from vise.input_set.input_options import CategorizedInputOptions


@dataclass()
class PriorInfo(MSONable):
    structure: Structure = None
    energy_per_atom: float = None
    band_gap: float = None
    total_magnetization: float = None
    data_source: str = None
    is_cluster: bool = None
    magnetization_criterion: float = defaults.integer_criterion
    band_gap_criterion: float = defaults.band_gap_criterion
    incar: dict = None

    def dump_yaml(self, filename: str = "prior_info.yaml") -> None:
        with open(filename, "w") as f:
            f.write(yaml.dump(self.as_dict()))

    @classmethod
    def load_yaml(cls, filename: str = "prior_info.yaml"):
        with open(filename, "r") as f:
            d = yaml.load(f, Loader=yaml.SafeLoader)
        d["incar"] = d.get("incar", None)
        return cls.from_dict(d)

    def dump_json(self, filename: str = "prior_info.json") -> None:
        with open(filename, "w") as fw:
            json.dump(self.as_dict(), fw, indent=2)

    @classmethod
    def load_json(cls, filename: str = "prior_info.json"):
        return loadfn(filename)

    @property
    def is_magnetic(self) -> Optional[bool]:
        try:
            return self.total_magnetization > self.magnetization_criterion
        except TypeError:
            return

    @property
    def has_band_gap(self) -> bool:
        return self.band_gap > self.band_gap_criterion

    @property
    def is_metal(self) -> bool:
        return not self.has_band_gap

#
# class PriorInfoFromCalcDir:
#     def __init__(self,
#                  prev_dir_path: Path,
#                  contcar: str = "CONTCAR",
#                  vasprun: str = "vasprun.xml",
#                  outcar: str = "OUTCAR",
#                  integer_criterion: float = defaults.integer_criterion):
#         self._vasprun = Vasprun(str(prev_dir_path / vasprun))
#         self._outcar = Outcar(str(prev_dir_path / outcar))
#         self._integer_criterion = integer_criterion
#
#         self.structure = Structure.from_file(str(prev_dir_path / contcar))
#
#     def vbm_cbm(self):
#         band_edges = VaspBandEdgeProperties(self._vasprun, self._outcar)
#         return [band_edges.vbm_info.energy, band_edges.cbm_info.energy]
#
#     def is_magnetic(self):
#         return self._outcar.total_mag > self._integer_criterion
#
#     def generate_input_options(self, **input_options_w_task_xc):
#         input_options_w_task_xc.update({"input_structure": self.structure,
#                                         "vbm_cbm": self.vbm_cbm,
#                                         "is_magnetization": self.is_magnetic})
#         return CategorizedInputOptions(**input_options_w_task_xc)