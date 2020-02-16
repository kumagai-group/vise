# -*- coding: utf-8 -*-
import json
import shutil
from copy import deepcopy
from pathlib import Path
from typing import List, Dict, Union

from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import Vasprun
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
        """

        Name is used for the plot of the chemical potential.

        Args:
            composition (Composition):
                Composition of the entry.
            total_energy (float):
                Energy of the entry. Usually final calculated energy from VASP.
            zero_point_vib (float):
                Zero point vibration energy.
            free_e_shift (float):
                Entropic contribution to the free energy
            name (str):
                Name of the Entry.
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
        return f"PDEntry : {self.name} with composition {self.composition}" \
               f" and energy = {self.energy:.4f}"


class FreeEnergyEntrySet(EntrySet):
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
                        temperature: float = 0,
                        pressure: float = REFERENCE_PRESSURE
                        ) -> "FreeEnergyEntrySet":

        energy_entries = []
        for d in directory_paths:
            logger.info(f"Parsing data in {d} ...")
            v: Vasprun = parse_file(Vasprun, Path(d) / vasprun)
            composition = v.final_structure.composition.reduced_formula
            kwargs = {}
            if parse_gas and composition in Gas.name_list():
                gas = Gas[composition]
                kwargs = \
                    {"zero_point_vib": gas.zero_point_vibrational_energy,
                     "free_e_shift": gas.free_e_shift(temperature, pressure)}

            energy_entries.append(
                FreeEnergyEntry(composition=composition,
                                total_energy=v.final_energy, **kwargs))

        return cls(energy_entries, temperature, pressure)

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


class ConstrainedFreeEnergyEntrySet(FreeEnergyEntrySet):
    """FreeEnergyEntrySet obtained at the constrained chemical potential."""
    def __init__(self,
                 entries: List[FreeEnergyEntry],
                 original_entries:  List[FreeEnergyEntry],
                 constr_chempot: Dict[str, float],
                 temperature: float = None,
                 pressure: float = None):
        """ """
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


def get_mp_materials(elements: List[str],
                     properties: List[str],
                     e_above_hull: float = 1e-4,
                     api_key=None):
    exclude_z = list(i for i in range(1, 100))  # from H to Cn
    excluded_elements = [str(Element.from_Z(z)) for z in exclude_z]
    for e in elements:
        excluded_elements.remove(e)
    with MPRester(api_key) as m:
        materials = \
            m.query(criteria={"elements": {"$in": elements,
                                           "$nin": excluded_elements},
                              "e_above_hull": {"$lte": e_above_hull}},
                    properties=properties)

    return materials


def make_poscars_from_mp(elements,
                         path: Path = Path.cwd(),
                         e_above_hull=0.01,
                         api_key=None,
                         add_molecules=True,
                         only_molecules=True) -> None:
    """

    Args:
        elements(list): like ["Cu", "O"]
        path(str):
        e_above_hull(float):
        api_key(str):
        add_molecules(bool):
        only_molecules:

    Returns:
        None
    """
    if not path.is_dir:
        raise NotADirectoryError(f"{path} is not directory.")

    mol_dir = Path(__file__).parent / "molecules"

    molecules_formula_list = []
    if add_molecules:
        for g in Gas:
            comp = Composition(str(g))
            if set([str(e) for e in comp.elements]) < set(elements):
                molecules_formula_list.append(comp.reduced_formula)
                dirname = path / f"mol_{str(comp)}"
                if dirname.exists():
                    logger.critical(f"{dirname} exists! So, skip creating it.")
                else:
                    dirname.mkdir()
                    shutil.copyfile(mol_dir / str(comp) / "POSCAR",
                                    dirname / "POSCAR")

    properties = ["task_id",
                  "full_formula",
                  "final_energy",
                  "structure",
                  "spacegroup",
                  "band_gap",
                  "total_magnetization",
                  "magnetic_type"]
    materials = get_mp_materials(elements, properties, e_above_hull, api_key)

    for m in materials:
        comp = Composition(m.pop("full_formula")).reduced_formula
        if only_molecules and comp in molecules_formula_list:
            continue
        m_path = path / f"{m['task_id']}_{comp}"
        m_path.mkdir()
        m.pop("structure").to(filename=m_path / "POSCAR")
        json_path = m_path / "prior_info.json"
        with open(str(json_path), "w") as fw:
            json.dump(m, fw)

