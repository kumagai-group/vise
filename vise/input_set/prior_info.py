# -*- coding: utf-8 -*-

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict

import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen.io.vasp import Vasprun, Outcar

from vise.defaults import defaults
from vise.input_set.input_options import CategorizedInputOptions
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.util.file_transfer import FileTransfers


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

    def input_options_kwargs(self):
        return {"input_structure": self.structure,
                "band_gap": self.band_gap,
                "vbm_cbm": self.vbm_cbm,
                "is_magnetization": self.is_magnetic}


class PriorInfoFromCalcDir(PriorInfo):
    def __init__(self,
                 prev_dir_path: Path,
                 vasprun: str = "vasprun.xml",
                 outcar: str = "OUTCAR",
                 file_transfer_type: Optional[Dict[str, str]] = None):

        vasprun = Vasprun(str(prev_dir_path / vasprun))
        outcar = Outcar(str(prev_dir_path / outcar))

        structure = vasprun.final_structure.copy()
        energy_per_atom = vasprun.final_energy / len(structure)
        band_edge_property = VaspBandEdgeProperties(vasprun, outcar)
        total_magnetization = outcar.total_mag

        super().__init__(structure=structure,
                         energy_per_atom=energy_per_atom,
                         band_gap=band_edge_property.band_gap,
                         vbm_cbm=band_edge_property.vbm_cbm,
                         total_magnetization=total_magnetization)

        self.file_transfer_type = FileTransfers(file_transfer_type,
                                                path=prev_dir_path)

