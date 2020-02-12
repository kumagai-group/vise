# -*- coding: utf-8 -*-
#  Copyright (c) Oba-group
#  Distributed under the terms of the MIT License.

from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union, Any

from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core.composition import Composition
from pymatgen.io.vasp import Vasprun

from vise.util.logger import get_logger
from vise.util.tools import parse_file
from vise.chempotdiag.gas import Gas

logger = get_logger(__name__)


class EnergyEntry(PDEntry):

    def __init__(self,
                 composition: Union[Composition, str],
                 total_energy: float,
                 zero_point_vib: float = None,
                 free_e_shift: float = None,
                 name: str = None,
                 attribute: object = None):
        """

        Args:
            composition:
            total_energy:
            zero_point_vib:
            free_e_shift:
            name:
            attribute:
        """
        super().__init__(composition, 0.0, name, attribute)
        self.total_energy = total_energy
        self.zero_point_vib = zero_point_vib
        self.free_e_shift = free_e_shift

        self.energy = total_energy
        if zero_point_vib:
            self.energy += zero_point_vib
        if free_e_shift:
            self.energy += free_e_shift

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "composition": self.composition.as_dict(),
                "total_energy": self.total_energy,
                "zero_point_vib": self.zero_point_vib,
                "free_e_shift": self.free_e_shift,
                "name": self.name,
                "attribute": self.attribute}

    @classmethod
    def from_dict(cls, d):
        d["composition"] = Composition(d["composition"])
        return cls(**d)


class EnergyEntries:
    def __init__(self,
                 energy_entries: List[EnergyEntry],
                 temperature: float = 0,
                 pressure: float = 1e+5):
        self.energy_entries = energy_entries
        self.temperature = temperature
        self.pressure = pressure

    @classmethod
    def from_vasp_files(cls,
                        directory_paths: List[Path],
                        vasprun: str = "vasprun.xml",
                        parse_gas: bool = True,
                        temperature: float = 0,
                        pressure: float = 1e+5,
                        ) -> "EnergyEntries":

        energy_entries = []
        for d in directory_paths:
            logger.info(f"Parsing data in {d} ...")
            vasprun = parse_file(Vasprun, Path(d) / vasprun)
            composition = vasprun.unit_cell_formula
            if parse_gas and composition in Gas.name_list():
                gas = Gas[composition]
                d = {"zero_point_vib": gas.zero_point_vibrational_energy,
                     "free_e_shift": gas.free_e_shift(temperature, pressure)}
                energy_entries.append(
                    EnergyEntry(composition=composition,
                                total_energy=vasprun.final_energy, **d))

        return cls(energy_entries, temperature, pressure)

    @classmethod
    def from_mp(cls,
                elements: List[str],
                align_substances: bool,
                substance_energies: dict) -> "EnergyEntries":
        """ """


    @classmethod
    def from_yaml(cls, filename: str) -> "EnergyEntries":
        pass

