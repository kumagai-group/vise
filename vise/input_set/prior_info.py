# -*- coding: utf-8 -*-

import json
from typing import Optional
import yaml
from dataclasses import dataclass

from monty.json import MSONable
from monty.serialization import loadfn

from vise.defaults import defaults


@dataclass()
class PriorInfo(MSONable):
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
            d = yaml.load(f, Loader=yaml.FullLoader)
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
