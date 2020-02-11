#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
import copy
from collections import defaultdict
import os
from typing import Sequence, Union, Optional, List, Dict, Tuple
from pathlib import Path
from collections import UserList

import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.ext.matproj import MPRester
from pymatgen.entries.computed_entries import ComputedEntry

from vise.chempotdiag.gas import Gas
from vise.config import MOLECULE_SUFFIX, REFERENCE_PRESSURE

# TODO: Once __eq__ is implemented, __hash__ has to be implemented?


ElemOrderType = Sequence[Union[str, Element]]


class Compound:
    def __init__(self,
                 name: Optional[str],
                 composition: Composition,
                 energy: float,
                 standard_energies: dict = None,
                 gas: Optional[Gas] = None):
        """ Create a Compound object.

        Args:
            name (None/str):
                Name of compound. Usually reduced formula.
            composition (Composition):
                Composition of compound, like {"Ba": 1, "Ti": 1, "O": 3}
            energy (float):
                Energy of compound.
            standard_energies (dict):
                Referenced standard energies for each element.
            gas (Gas):
                If compound is gas, related gas data (like Gas.O2) will be used
                when temperature and pressure are considered.
        """
        self.name = name
        self.composition = composition.fractional_composition
        self.energy = energy / composition.num_atoms
        self.standard_energies = standard_energies
        self.gas = gas

    @classmethod
    def from_vasp_files(cls,
                        path: Path,
                        vasprun: str = "vasprun.xml",
                        is_molecule_gas: bool = True,
                        name: Optional[str] = None,
                        energy_shift: float = 0) -> "Compound":
        """Create Compound object from vasp calculation.

        Args:
            path (Path):
                Path to the calculation directory.
            vasprun (str):
                Name of the vasprun.xml file.
            is_molecule_gas (bool):
                Read molecule data as gas. If directory name contains
                "molecule", zero-point energy and free energy are considered,
                via database.
            name (None/str):
                Name of the compound.
            energy_shift (float):
                Energy shift value in eV.

        Returns:
            (Compound) Compound object from vasp files.
        """
        v = Vasprun(str(path / vasprun))
        composition = v.structures[-1].composition

        gas = None
        if is_molecule_gas and "molecule" in path.name:
            gas_formula = composition.get_reduced_formula_and_factor()[0]
            # search gas
            for g in Gas:
                if str(g) == gas_formula:
                    gas = g
                    break
            if gas is None:
                raise ValueError(f"{gas_formula} not found in molecule db")

        energy = v.final_energy + energy_shift
        name = name or composition.reduced_formula

        return cls(name, composition, energy, gas=gas)

    @classmethod
    def from_mp_database(cls, mp_id: str) -> "Compound":
        """
        Args:
            mp_id (str): material_id of materials project database

        Returns:
            Compound object.
        """
        mpr = MPRester()
        mp_dat: ComputedEntry = mpr.get_entry_by_material_id(mp_id)
        energy = mp_dat.energy
        composition = mp_dat.composition

        return cls(f"{mp_id}-{composition}", composition, energy)

    @property
    def elements(self) -> List[Element]:
        """(list of Element) List of elements of the composition vector. """
        return self.composition.elements

    @property
    def dim(self) -> int:
        """(int) Number of considered atoms. """
        return len(self.elements)

    def composition_vector(self, element_list: ElemOrderType) -> np.ndarray:
        """Create composition as numpy array

        The composition order is assumed to be the same as given elements.

        x = Composition("N2").fractional_composition
        x["N"] = 1.0
        x["O"] = 0

        Args:
            element_list:

        Returns:
            np.ndarray for the composition
        """
        return np.array([self.composition[e] for e in element_list])

    def gas_energy_shift(self, temperature: float, pressure: float) -> float:
        """Energy shift caused by zero point vibration and entropy.

        Args:
            temperature(float):
                Temperature in K.
            pressure(float):
                Pressure in Pa.

        Returns (float):
            Value of gas energy shift in eV/atom.
        """
        if self.gas is None:
            return 0.0
        else:
            return self.gas.energy_shift(temperature, pressure)

    def free_energy(self,
                    temperature: Optional[float] = None,
                    pressure: Optional[float] = None) -> float:
        """Free energy

        Args:
            temperature(float or None): (K)
            pressure(float or None): (Pa)

        Returns (float):
            Value of free energy shift in eV/atom.
        """
        return self.energy + self.gas_energy_shift(temperature, pressure)

    @property
    def standardized_energy(self):
        if self.standard_energies is None:
            raise ValueError(f"standard_energies is not set yet.")

        sum_standard_energy = 0.0
        for e, c in zip(self.elements, self.composition_vector(self.elements)):
            sum_standard_energy += self.standard_energies[e] * c

        return self.energy - sum_standard_energy

    def __repr__(self):
        return (f"Name: {self.name}, Composition: {self.composition}, "
                f"Energy: {self.energy}, Gas: {self.gas}")

    def __eq__(self, other: "Compound"):
        """Used for sorting. """
        if not isinstance(other, Compound):
            raise TypeError("Compared obj is not Compound object.")
        if self.name != other.name:
            return False
        elif self.energy != other.energy:
            return False
        elif any([a != b for a, b in zip(self.composition, other.composition)]):
            return False

        return True

    def almost_equal(self, other: "Compound", atol: float = 1e-5) -> bool:
        """
        Args:
            other (Compound): Compared compound.
            atol: Absolute tolerance.

        Returns (bool):
            If self and other almost_equals.
        """
        return self.composition.almost_equals(other.composition, atol=atol)

    def __ne__(self, other: "Compound"):
        """Used for sorting. """
        return not self == other

    def __lt__(self, other: "Compound"):
        """Used for sorting. """
        if not isinstance(other, Compound):
            raise TypeError("Compound can't be compared with another class.")

        if self.name != other.name:
            return self.name < other.name
        elif self.energy != other.energy:
            return self.energy < other.energy
        return False


class DummyCompoundForDiagram(Compound):
    """ Object for dummy compound for drawing boundary of diagram.

    This class is intended to be used for scipy.spatial.halfspaces.
    It will mean halfspace expressed by (x > boundary_energy),
    but halfspace has to be notated by (Ax + b < 0), then composition is -1.
    """

    # HACK: This __init__ is needed to allow composition -1,
    # without standardization.
    # def __init__(self, name: str, composition: Composition, energy: float):
    #     """ Constructor of dummy boundary. """
    #     super().__init__(name, composition, energy)

    @classmethod
    def construct_boundary(cls,
                           element: Union[Element, str],
                           energy: float) -> "Compound":
        name = f"{element} boundary"
        composition = Composition({Element(element): -1}, allow_negative=True)
        return cls(name, composition, energy)

    @classmethod
    def from_vasp_files(cls, *args, **kwargs):
        raise NotImplementedError("DummyCompound can't be made from files.")


class CompoundsList(UserList):

    def __init__(self):
        super().__init__(self)
        for c in self:
            if not isinstance(c, Compound):
                raise TypeError("c must not contain objects except Compound.")

    def __str__(self):
        return super(CompoundsList, self).__str__()

    def append(self, compound: Compound):
        super().append(compound)

    def extend(self, compounds: "CompoundsList"):
        super().extend(compounds)

    def __add__(self, compounds: "CompoundsList"):
        super().__add__(compounds)

    @property
    def elements(self) -> List[str]:
        elements = set()
        for c in self:
            elements |= c.elements
        return list(sorted(elements))

    @property
    def dim(self) -> int:
        return len(self.elements)

    def standard_energy_array(self,
                              temperature: Optional[float],
                              pressure: Optional[Dict[str, float]]
                              ) -> Tuple["CompoundsList", Dict[str, float]]:
        """Energies of stable simple substance are set to 0

        Args:
            temperature (float):
            pressure (Dict[str, int]):

        Returns:
            Tuple of new CompoundsList and dict of element energy.

        """
        copied_obj = copy.deepcopy(self)
        comp_list: List[Composition] = [c.composition for c in copied_obj]

        pressure = pressure or defaultdict(lambda: REFERENCE_PRESSURE)
        energy_list = np.array(copied_obj.free_energies(temperature, pressure))

        element_energy = {e: float("inf") for e in self.elements}
        # search minimum energies of substances
        for c, e in zip(comp_list, energy_list):
            # index = np.where(abs(c - 1) < 1e-8)[0]
            if c.is_element:
                elem = c.elements[0]
                e_per_atom = e / c.num_atoms
                if e_per_atom < element_energy[elem]:
                    element_energy[elem] = e_per_atom

        if any(en == float("inf") for en in element_energy.values()):
            not_found = \
                [e for e, en in element_energy.items() if en == float('inf')]
            raise ValueError(f"Standardization of energies of compounds failed "
                             f"because calculation of simple element "
                             f"{not_found} was not found.")

        for c in self:
            for e in c.elements:
                c.standard_energy = c.composition[e] * element_energy[e]

        return copied_obj, element_energy

    def gas_energy_shifts(self,
                          temperature: Optional[float],
                          pressure: Optional[Dict[str, float]]
                          ) -> List[float]:
        """Energy_shift of compounds for the gases

        Except for Gas, zero is inserted.

        Args:
            temperature (float or None)
            pressure (Dict[str, float] or None)

        Returns (list of float):
            List of gas energy shifts.
        """
        if len([o is None for o in [temperature, pressure]]) == 1:
            raise ValueError(f"Only one of temperature ({temperature}) "
                             f"and pressure ({pressure}) is specified."
                             f"Specify two or zero of these option.")
        energy_shifts = []
        for c in self:
            if not c.gas:
                energy_shift = 0
            else:
                energy_shift = c.gas_energy_shift(temperature,
                                                  pressure[c.gas.formula])
            energy_shifts.append(energy_shift)

        return energy_shifts

    def free_energies(self,
                      temperature: Optional[float],
                      pressure: Optional[Dict[str, float]]
                      ) -> List[float]:
        """Shifted free energy by temperature and pressure for gas phases,

        Args:
            temperature (float):
            pressure (Dict[str, float]):

        Returns (list of float):
            List of free energies.
        """
        if pressure is None:
            pressure = defaultdict(lambda: REFERENCE_PRESSURE, {})
        comp_energy = zip(self, self.gas_energy_shifts(temperature, pressure))
        return [c.energy + es for c, es in comp_energy]

    def get_indices_and_compounds(self, compound_name: str):
        """
        Find object of Compound from self by name(str) of compound.

        Args:
            compound_name (str): name of compound.

        Returns:
            (None/list) Matched compound data. If no compounds match, return None.
        """
        compound_name = Composition(compound_name).reduced_formula
        return [(i, c) for i, c in enumerate(self)
                if Composition(c.name).reduced_formula == compound_name]

    @classmethod
    def from_file(cls, file_name: str):
        """
        Create a object of CompDat from file.

        Args:
            file_name (str): File name of information of energies.
            The style of file is like below
            Mg-O
            Mg    1  0   -2.11327587
            Mg    1  0   -1.2342344
            O     0  1   -2.46256238
            Mg2O  2  1  -14.46256238 ...
            Note that the energies of element substances are essential.
        """
        with open(file_name) as f:
            lines = f.readlines()
        elements = lines[0].strip().split('-')
        compounds_list = []
        num_elements = len(elements)
        for l in lines[1:]:
            line = l.split()
            name = line[0]
            num_atoms = np.array([float(i) for i in line[1:num_elements + 1]])
            num_total_atoms = np.sum(num_atoms)
            composition = \
                Composition({Element(e): amt
                             for e, amt in
                             zip(elements, num_atoms / num_total_atoms)})
            energy = float(line[num_elements + 1]) / num_total_atoms
            compounds_list.append(Compound(name, composition, energy))
        compounds_list = CompoundsList(compounds_list)
        return cls(compounds_list)

    # TODO: make unittest
    def to(self, file_name: Optional[str] = None):
        """
        CompoundsList to str. If filename is provided, this method outputs file.

        Args:
            file_name (str): Name of file to output.
                             If it is not provided, string will be returned.
        Returns:
            (None/str): If file_name is provided, returns None.
            Otherwise, returns str.
            The style of file or str must be consistent with from_file.
        """
        return_str = "-".join(str(e) for e in self.elements)
        for compound in self:
            return_str += "\n"
            return_str += f" {compound.name}"
            for c in compound.composition:
                return_str += f" {c}"
            return_str += f" {compound.energy}"
        if file_name:
            with open(file_name, 'w') as f:
                f.write(return_str)
            return None
        else:
            return return_str

    @classmethod
    def from_vasp_files(cls,
                        paths: List[str],
                        vasprun: str = "vasprun.xml",
                        energy_shift_dict: Optional[Dict[str, float]] = None):
        """
        Args:
            paths (Path):
                Paths to the calculation difectories.
            vasprun (str):
                Name of the vasprun.xml file.
            energy_shift_dict (dict):
                Energy shift per cell in eV. Keys are directory names,
                e.g. {N2_molecule: 1}
        Returns:
            (CompoundsList) CompoundsList object from vasp files.
        """
        energy_shift_dict = energy_shift_dict or {}
        energy_shift_dict = defaultdict(float, energy_shift_dict)

        compounds_list = cls([])
        for path in paths:
            compound = Compound.from_vasp_files(
                path=Path(path),
                vasprun=vasprun,
                energy_shift=energy_shift_dict[str(path)])
            compounds_list.append(compound)

        return compounds_list

    @classmethod
    def from_vasp_and_materials_project(cls,
                                        vasp_target_poscar: str,
                                        vasp_target_output: str,
                                        vasp_element_poscar: List[str],
                                        vasp_element_output: List[str],
                                        fmt: str = "outcar",
                                        energy_shift_dict:
                                        Optional[Dict[str, float]] = None):
        """

        Args:
            vasp_target_poscar (str): path
            vasp_target_output (str): path
            vasp_element_poscar (list): path
            vasp_element_output (list): path
            fmt (str):
            energy_shift_dict (dict):

        Returns:
            (CompoundsList) CompoundsList object from materials project.

        """
        if energy_shift_dict is None:
            energy_shift_dict = defaultdict(float, {})
        # to match energy_shift_dict key,
        # output_paths is transformed to abs_path
        # (for example, it is preferable to match "OUTCAR" and "./OUTCAR")
        energy_shift_dict = {os.path.abspath(k): v
                             for k, v in energy_shift_dict.items()}
        energy_shift_dict = defaultdict(float, energy_shift_dict)

        # compounds_list = cls([], temperature=temperature, pressure=pressure)
        compounds_list = cls([])

        # target compound
        target_s = Poscar.from_file(vasp_target_poscar).structure
        elements = target_s.composition.elements
        compound = Compound.\
            from_vasp_files(vasp_target_poscar,
                            vasp_target_output,
                            fmt=fmt,
                            elements=elements,
                            energy_shift=
                                        energy_shift_dict[vasp_target_output])
        compounds_list.append(compound)

        vasp_element_energy = defaultdict(float, {})
        element_set = set()
        # vasp elements
        for p, o in zip(vasp_element_poscar, vasp_element_output):
            compound = Compound.\
                from_vasp_files(p,
                                o,
                                fmt=fmt,
                                elements=elements,
                                energy_shift=energy_shift_dict[o],
                                is_molecule_gas=False)
            p = Poscar.from_file(p)
            if not p.structure.composition.is_element:
                raise ValueError(f"Unexpected case: VASP calculation of "
                                 f"not simple element, "
                                 f"{Poscar.from_file(p).structure.composition}",
                                 f"{p}")
            poscar_element = p.structure.composition.elements[0]
            if poscar_element in element_set:
                raise ValueError(f"Unexpected case:"
                                 f"2 or more results of same species, "
                                 f"{poscar_element}")
            vasp_element_energy[str(poscar_element)] = compound.energy
            element_set.add(poscar_element)

        # from mp data
        # MPRester elements and make alignment
        mpr = MPRester()
        mp_element_energy = {}
        for e in elements:
            q = {"elements": [str(e)],
                 "e_above_hull": 0}
            data = mpr.query(q, ["energy_per_atom"])[0]
            mp_element_energy[str(e)] = data["energy_per_atom"]

        mp_to_vasp = \
            {e: (vasp_element_energy[str(e)] - mp_element_energy[str(e)])
             for e in elements}

        def align_mp_to_vasp(mp_energy: float, composition_: Composition):
            aligned = mp_energy
            for elem_, amount in composition_.items():
                aligned += mp_to_vasp[elem_] * amount
            return aligned

        q = {"elements":
             {"$in": [str(e) for e in elements],
              "$nin": [str(e) for e in Element if e not in elements]
              },
             "e_above_hull": 0
             }
        for mp_dat in mpr.query(q, ["energy_per_atom", "pretty_formula",
                                    "e_above_hull", "material_id"]):
            energy: float = mp_dat["energy_per_atom"]
            composition = Composition(mp_dat["pretty_formula"])
            composition = composition / composition.num_atoms
            composition_vector = []
            for elem in elements:
                if elem in elements:
                    composition_vector.append(composition[elem])
                else:
                    composition_vector.append(0)
            compound = Compound(mp_dat["pretty_formula"],
                                elements,
                                composition_vector,
                                align_mp_to_vasp(energy, composition))
            compounds_list.append(compound)
        return compounds_list

    # TODO: unittest
    def almost_equal(self, other: "CompoundsList", tol: float = 1e-5):
        """
        Check whether two CompoundsList are almost same.

        Args:
            other (CompoundsList): Compared compounds_list.
            tol (float): Tolerance for numeric error of energy and composition.

        Returns (bool): If self almost equals to other.
        """
        if len(self) != len(other):
            return False
        else:
            for c1, c2 in zip(sorted(self), sorted(other)):
                if not c1.almost_equal(c2, atol=tol):
                    return False
        return True

