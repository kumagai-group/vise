# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
from pathlib import Path
from typing import List, Dict, Union
import re

from monty.json import MSONable

from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core.composition import Composition
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.io.vasp import Vasprun

from vise.util.mp_tools import get_mp_materials
from vise.chempotdiag.gas import Gas
from vise.config import REFERENCE_PRESSURE
from vise.util.logger import get_logger
from vise.util.tools import parse_file

logger = get_logger(__name__)


class FreeEnergyEntry(PDEntry):

    def __init__(self,
                 composition: Union[Composition, str],
                 total_energy: float,
                 zero_point_vib: float = None,
                 free_e_shift: float = None,
                 name: str = None,
                 data: dict = None,
                 attribute: object = None):
        """Free energy version of PDEntry with zero_point_vib and free_e_shift

        Name is used for plotting chemical potential and phase diagrams. It is
        especially important for CompoundPhaseDiagram, as the composition is
        described with DummySpecie.

        Args:
            composition (Composition/str):
                Composition of the entry.
            total_energy (float):
                Energy of the entry. Usually final calculated energy from VASP.
            zero_point_vib (float):
                Zero point vibration energy.
            free_e_shift (float):
                Entropic contribution to the free energy
            name (str):
                Name of the Entry. Mostly used for drawing the figure.
            data (dict):
                An optional dict of any additional data associated
                with the entry. Defaults to None.
            attribute:
                Optional attribute of the entry. This can be used to specify
                that the entry is a newly found compound, or to specify a
                particular label for the entry, or else ... Used for further
                analysis and plotting purposes. An attribute can be anything
                but must be MSONable.
        """
        super().__init__(composition, 0.0, name, attribute)
        self.total_energy = total_energy
        self.zero_point_vib = zero_point_vib
        self.free_e_shift = free_e_shift
        self.data = data
        self.energy = total_energy

        if zero_point_vib:
            self.energy += zero_point_vib
        if free_e_shift:
            self.energy += free_e_shift

    def as_dict(self):
        return {"composition": self.composition.as_dict(),
                "total_energy": self.total_energy,
                "zero_point_vib": self.zero_point_vib,
                "free_e_shift": self.free_e_shift,
                "name": self.name,
                "data": self.data,
                "attribute": self.attribute}

    @classmethod
    def from_dict(cls, d: dict):
        kwargs = deepcopy(d)
        kwargs["composition"] = Composition(d["composition"])
        kwargs.pop("@module", None)
        kwargs.pop("@class", None)
        return cls(**kwargs)

    def __repr__(self):
        zero_point_vib = \
            f"{self.zero_point_vib:.3f}" if self.zero_point_vib else "none"
        free_e_shift = \
            f"{self.free_e_shift:.3f}" if self.free_e_shift else "none"
        return f"PDEntry : {self.name} with composition {self.composition}. " \
               f"Energy: {self.energy:.3f}. " \
               f"Zero point vib energy: {zero_point_vib}. " \
               f"Free energy contribution: {free_e_shift}. " \
               f"Data: {self.data}"


class FreeEnergyEntrySet(EntrySet, MSONable):
    """List of FreeEnergyEntries with temperature and pressure."""
    def __init__(self,
                 entries: List[FreeEnergyEntry],
                 temperature: float = None,
                 pressure: float = None):
        super().__init__(entries)
        self.temperature = temperature
        self.pressure = pressure

    @classmethod
    def from_vasp_files(cls,
                        directory_paths: List[Path],
                        vasprun: str = "vasprun.xml",
                        parse_gas: bool = True,
                        temperature: float = 0.0,
                        partial_pressures: Dict[str, float] = None,
                        ignore_file_not_found: bool = True,
                        ) -> "FreeEnergyEntrySet":
        """ Constructs class object from vasp output files.

        Args:
            directory_paths (list):
                List of relative directory paths, containig vasp results.
            vasprun (str):
                Name of the vasprun.xml file.
            parse_gas (bool):
                Whether to parse gases as gas phases if the composition exists
                in the Gas.name_list().
            temperature (float):
                Temperature in K.
            partial_pressures (dict):
                Dict of species as keys (str) and pressures in Pa as values.
                Example: {"O2": 2e5, "N2": 70000}
            ignore_file_not_found (bool):

        Returns:
            FreeEnergyEntrySet class object.
        """
        energy_entries = []
        for d in directory_paths:
            logger.info(f"Parsing data in {d} ...")
            try:
                v: Vasprun = parse_file(Vasprun, Path(d) / vasprun)
            except FileNotFoundError:
                if ignore_file_not_found:
                    logger.critical(f"{d} is not parsed as vasprun.xml does "
                                    f"not exist in it.")
                    continue
                else:
                    raise
            composition = v.final_structure.composition.formula

            kwargs = {}
            mol_dir = re.match(r"^mol_", str(d))
            if parse_gas and mol_dir:
                if composition not in Gas.name_list():
                    logger.critical(f"{d} does not exist in Gas.name_list, so "
                                    f"skipp to generate the zero point vib "
                                    f"and free energies.")

                if partial_pressures is None:
                    logger.info(f"Pressure of {composition} is set to "
                                f"{REFERENCE_PRESSURE} (Pa).")
                    pressure = REFERENCE_PRESSURE
                else:
                    try:
                        pressure = partial_pressures[composition]
                    except KeyError as e:
                        msg = f"Partial pressure for {composition} is not set."
                        raise Exception(msg) from e
                        
                gas = Gas[composition]
                zpve = gas.zero_point_vibrational_energy
                free_e = gas.free_e_shift(temperature, pressure)
                logger.warning(
                    f"{composition} is parsed as gas. "
                    f"Pressure: {pressure}, "
                    f"temperature: {temperature}, "
                    f"zero point vib energy: {zpve}, "
                    f"free energy: {free_e}")

                kwargs = {"zero_point_vib": zpve, "free_e_shift": free_e}

            energy_entries.append(
                FreeEnergyEntry(composition=composition,
                                total_energy=v.final_energy, **kwargs))

        return cls(energy_entries, temperature, partial_pressures)

    @classmethod
    def from_mp(cls, elements: List[str]) -> "FreeEnergyEntrySet":
        """Obtain the energies from Materials Project."""
        properties = ["task_id", "full_formula", "final_energy"]
        materials = get_mp_materials(elements, properties=properties)

        energy_entries = []
        for m in materials:
            energy_entries.append(
                FreeEnergyEntry(composition=m["full_formula"],
                                total_energy=m["final_energy"],
                                data={"mp_id": m["task_id"]}))

        return cls(energy_entries)

    def __repr__(self):
        lines = [f"pressure {self.pressure}",
                 f"temperature {self.temperature}"]
        lines.extend([str(i) for i in self.entries])

        return "\n".join(lines)


class ConstrainedFreeEnergyEntrySet(FreeEnergyEntrySet):
    def __init__(self,
                 entries: List[FreeEnergyEntry],
                 original_entries:  List[FreeEnergyEntry],
                 constr_chempot: Dict[str, float],
                 temperature: float = None,
                 pressure: float = None):
        """FreeEnergyEntrySet obtained at the constrained chemical potential.

       Args:
           entries (List[FreeEnergyEntry]):
               Entries for the
           original_entries (List[FreeEnergyEntry]):
           constr_chempot (Dict[str, float]):
           temperature (float):

           pressure (float):

        """
        super().__init__(entries, temperature, pressure)
        self.original_entries = original_entries
        self.constr_chempot = constr_chempot

    @classmethod
    def from_entry_set(cls,
                       entry_set: FreeEnergyEntrySet,
                       constr_chempot: Dict[str, float]
                       ) -> "ConstrainedFreeEnergyEntrySet":
        """

        Args:
            entry_set:
                Original FreeEnergyEntrySet.
            constr_chempot:
                Constrained chemical potentials in absolute scale.

        Returns:
            ConstrainedFreeEnergyEntrySet object.
        """
        new_entries = []
        for e in entry_set.entries:
            new_entry = deepcopy(e)
            comp = e.composition
            sub_comp = Composition({e: comp[e] for e in constr_chempot})
            if comp == sub_comp:
                continue
            new_entry.composition = comp - sub_comp
            sub_energy = sum([comp[e] * c for e, c in constr_chempot.items()])
            new_entry.energy = e.energy - sub_energy
            new_entries.append(new_entry)

        return cls(entries=new_entries,
                   original_entries=entry_set.entries,
                   constr_chempot=constr_chempot,
                   temperature=entry_set.temperature,
                   pressure=entry_set.pressure)


