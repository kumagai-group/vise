# -*- coding: utf-8 -*-

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen.io.vasp import Vasprun, Outcar

from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.defaults import defaults


@dataclass()
class PriorInfo(MSONable):
    structure: Structure = None
    energy_per_atom: float = None
    band_gap: float = None
    vbm_cbm: list = field(default_factory=list)
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

    @property
    def input_options_kwargs(self):
        return {"vbm_cbm": self.vbm_cbm,
                "is_magnetization": self.is_magnetic}


def prior_info_from_calc_dir(prev_dir_path: Path,
                             vasprun: str = "vasprun.xml",
                             outcar: str = "OUTCAR"):

    vasprun = Vasprun(str(prev_dir_path / vasprun))
    outcar = Outcar(str(prev_dir_path / outcar))

    structure = vasprun.final_structure.copy()
    energy_per_atom = vasprun.final_energy / len(structure)
    band_edge_property = VaspBandEdgeProperties(vasprun, outcar)
    total_magnetization = outcar.total_mag

    return PriorInfo(structure=structure,
                     energy_per_atom=energy_per_atom,
                     band_gap=band_edge_property.band_gap,
                     vbm_cbm=band_edge_property.vbm_cbm,
                     total_magnetization=total_magnetization)



