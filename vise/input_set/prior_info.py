# -*- coding: utf-8 -*-
import json
from typing import Optional

from monty.json import MSONable
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class PriorInfo(MSONable):
    """ Prior information for controlling parameters in DFT calculations """
    def __init__(self,
                 energy_per_atom: float = None,
                 band_gap: float = None,
                 total_magnetization: float = None,
                 data_source: str = None,
                 is_cluster: bool = None,
                 mag_threshold: float = 0.05,
                 band_gap_threshold: float = 0.1):
        """
            Args:
                energy_per_atom (float):
                    Energy per atom calculated in the data_source.
                band_gap (float):
                    Band gap calculated in the data_source.
                total_magnetization (float):
                    Total total_magnetization in the data_source.
                data_source (str):
                    The data source
                is_cluster (bool):
                    Whether the system is molecule or not.
                mag_threshold (float):
                    Threshold to judge if the system is magnetic.
                band_gap_threshold (float):
                    Threshold to judge if the system is metal.
        """
        self.energy_per_atom = energy_per_atom
        self.band_gap = band_gap
        self.total_magnetization = total_magnetization
        self.data_source = data_source
        self.is_cluster = is_cluster
        self.mag_threshold = mag_threshold
        self.band_gap_threshold = band_gap_threshold

    def __repr__(self):
        outs = ["energy_per_atom: " + str(self.energy_per_atom),
                "band gap: " + str(self.band_gap),
                "total magnetization: " + str(self.total_magnetization),
                "data source: " + str(self.data_source),
                "is cluster: " + str(self.is_cluster)]
        return "\n".join(outs)

    def dump_json(self, filename: str = "prior_info.json") -> None:
        with open(filename, "w") as fw:
            json.dump(self.as_dict(), fw, indent=2)

    @classmethod
    def load_json(cls, filename: str = "prior_info.json"):
        return loadfn(filename)

    @property
    def is_magnetic(self) -> Optional[bool]:
        try:
            return self.total_magnetization > self.mag_threshold
        except TypeError:
            return

    @property
    def has_band_gap(self) -> bool:
        return self.band_gap > self.band_gap_threshold

    @property
    def is_metal(self) -> bool:
        return self.band_gap < self.band_gap_threshold
