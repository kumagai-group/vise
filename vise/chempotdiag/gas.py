#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
from abc import ABCMeta, abstractmethod
from typing import Union, List, Tuple, Dict, Any
from enum import Enum
import json
import math
import os
import re
from collections import defaultdict
from typing import Optional

from scipy import constants
from ruamel import yaml

from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.inputs import Poscar
from chempotdiag.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE

"""
Data is from NIST Chemistry WebBook, SRD 69
"""


class FundamentalFrequency:

    def __init__(self, wave_number: float, degeneration: int):
        """

        Args:
            wave_number (float): frequency (cm^-1)
            degeneration (int): number of degeneration
        """
        self._wave_number = wave_number
        self._degeneration = degeneration

    def __str__(self):
        return f"wave_number: {self.wave_number} /cm, " \
               f"degeneration: {self.degeneration}"

    @property
    def wave_number(self) -> float:
        """

        Returns (float): cm^(-1)

        """
        return self._wave_number

    @property
    def degeneration(self) -> int:
        """

        Returns (int): number of degeneration

        """
        return self._degeneration

    def as_dict(self) -> Dict[str, Union[float, int]]:
        return {"wave_number": self._wave_number,
                "degeneration": self._degeneration}

    @property
    def zero_point_vibrational_energy(self) -> float:
        """

        Returns (float): (eV)

        """
        h = constants.physical_constants["Planck constant in eV s"][0]
        # h/2 * degneration *
        # wave_number(cm^-1) * 100(cm^-1 -> m^-1) * c(speed_of_light)
        return 1 / 2 * h * self._degeneration * \
            self._wave_number * 1e+2 * constants.speed_of_light


class FundamentalFrequencies:

    def __init__(self,
                 frequencies: List[FundamentalFrequency],
                 reference: str):
        """
        Args:
            frequencies (list of FundamentalFrequency): Fundamental frequencies.
            reference (str): Citation of the data.
        """
        self._frequencies = frequencies
        self._reference = reference

    def __iter__(self):
        return self._frequencies.__iter__()

    @classmethod
    def from_yaml(cls, file_path: str) -> "FundamentalFrequencies":
        """

        Args:
            file_path (str):

        Returns (FundamentalFrequencies):

        """

        with open(file_path, 'r') as fr:
            read_data = yaml.safe_load(fr)

        frequencies = \
            [FundamentalFrequency(f["wave_number"], f["degeneration"])
             for f in read_data["frequencies"]]
        try:
            reference = read_data["reference"]
        except KeyError:
            reference = None

        return cls(frequencies, reference)

    @property
    def zero_point_vibrational_energy(self) -> float:
        """

        Returns (float): E_ZPVE (eV)

        """
        return sum(f.zero_point_vibrational_energy for f in self._frequencies)


class AbstractThermodynamicsFunction(metaclass=ABCMeta):
    # TODO: implement abstract factory

    @abstractmethod
    def heat_capacity(self, temperature: float) -> float:
        pass

    @abstractmethod
    def standard_enthalpy(self, temperature: float) -> float:
        pass

    @abstractmethod
    def standard_entropy(self, temperature: float) -> float:
        pass

    @abstractmethod
    def r_ln_p_p0(self, pressure: float) -> float:
        pass

    @property
    @abstractmethod
    def enthalpy_zero(self) -> float:
        """

        Returns: Energy at pressure = self._reference_pressure, T=0K

        """


class ShomateParameters:
    """
        Parameters of Shomate Equations.
        Contains range of temperature like (100, 700).
    """

    def __init__(self,
                 temperature_range: Tuple[float, float],
                 a: float, b: float, c: float, d: float,
                 e: float, f: float, g: float, h: float):
        """

        Args:
            temperature_range (tuple): Unit is (K).
            a (float): parameter A.
            b (float): parameter B.
            c (float): parameter C.
            d (float): parameter D.
            e (float): parameter E.
            f (float): parameter F.
            g (float): parameter G.
            h (float): parameter H.
        """
        self._temperature_range = temperature_range
        self._a = a
        self._b = b
        self._c = c
        self._d = d
        self._e = e
        self._f = f
        self._g = g
        self._h = h

    @property
    def temperature_range(self) -> Tuple[float, float]:
        return self._temperature_range

    @property
    def min_temperature(self) -> float:
        return self._temperature_range[0]

    @property
    def max_temperature(self) -> float:
        return self._temperature_range[1]

    @property
    def a(self) -> float:
        return self._a

    @property
    def b(self) -> float:
        return self._b

    @property
    def c(self) -> float:
        return self._c

    @property
    def d(self) -> float:
        return self._d

    @property
    def e(self) -> float:
        return self._e

    @property
    def f(self) -> float:
        return self._f

    @property
    def g(self) -> float:
        return self._g

    @property
    def h(self) -> float:
        return self._h

    def can_apply(self, temp: float) -> bool:
        if self.min_temperature <= temp <= self.max_temperature:
            return True
        else:
            return False

    def __str__(self):
        return f"Temperature_range: {self.temperature_range}\n"\
               f"A: {self.a}\n"\
               f"B: {self.b}\n"\
               f"C: {self.c}\n"\
               f"D: {self.d}\n"\
               f"E: {self.e}\n"\
               f"F: {self.f}\n"\
               f"G: {self.g}\n"\
               f"H: {self.h}\n"


class ShomateThermodynamicsFunction(AbstractThermodynamicsFunction):
    """
    Shomate thermodynamics Function.
    """

    def __init__(self,
                 func_param_list: List[ShomateParameters],
                 enthalpy_zero: float,
                 reference_pressure: float = REFERENCE_PRESSURE):
        """
        Args:
            func_param_list (list): parameters of functions.
                                    Must contain _Shomate_parameters.
            enthalpy_zero (float): enthalpy (J/mol) at T=0K,
                                   p=reference_pressure
            reference_pressure (float): reference_pressure (Pa)

        """
        self._params = func_param_list
        self._enthalpy_zero = enthalpy_zero
        self._reference_pressure = reference_pressure

    @classmethod
    def from_nist_table(cls, file_path: str) -> "ShomateThermodynamicsFunction":
        """
        shomate_nist.dat is written e.g.
        Enthalpy_zero   -8.825
        Temperature (K)	298. - 6000.
        A	31.44510
        B	8.413831
        C	-2.778850
        D	0.218104
        E	-0.211175
        F	-10.43260
        G	237.2770
        H	0.000000
        Reference	Chase, 1998
        Comment	Data last reviewed in June, 1982

        where Enthalpy_zero is kJ/mol

        Args:
            file_path (str):

        Returns:

        """
        # TODO: unit of file is confusing
        # Is there better expression with regular expression?
        temperature_ranges = []
        param_dicts = []
        enthalpy_zero = None
        with open(file_path, "r") as fr:
            for line in fr:
                if enthalpy_zero is None:
                    is_tmp_line = \
                        bool(re.match(r"Enthalpy_zero", line))
                    if is_tmp_line:
                        matched = re.findall(r"[\-\d\.]+", line)
                        enthalpy_zero = float(matched[0]) / 1000
                    continue

                if not temperature_ranges:
                    is_tmp_line = \
                        bool(re.match(r"Temperature\s+\(K\)\s+", line))
                    if is_tmp_line:
                        matched = re.findall(r"[\d\.]+", line)
                        for i in range(int(len(matched)/2)):
                            temperature_ranges.append(
                                (float(matched[2*i]), float(matched[2*i+1])))
                        param_dicts = \
                            [dict() for _ in range(len(temperature_ranges))]
                else:
                    symbol_match = re.search(r"\s*([A-H])\s+", line)
                    if symbol_match:
                        symbol = symbol_match.groups()[0].lower()
                        values = re.findall(r"[\-\d\.]+", line)
                        if len(values) != len(temperature_ranges):
                            raise ValueError(f"Number of temperature_ranges"
                                             f"({len(temperature_ranges)})"
                                             f"and parameters ({line})"
                                             f"is not consistent.")
                        for i in range(len(values)):
                            param_dicts[i][symbol] = float(values[i])
        return cls([ShomateParameters(tr, **pd)
                    for tr, pd in zip(temperature_ranges, param_dicts)],
                   enthalpy_zero)

    def params(self, temperature: float) -> ShomateParameters:
        """
        Return parameters of function at temperature like {a:1.1, b:2.2, ...}
        Args:
            temperature:
        Returns:

        """
        params = None
        for p in self._params:
            if p.can_apply(temperature):
                params = p
        if params is None:
            raise ValueError(f"Temperature {temperature} is out of temperature "
                             f"range to apply")
        return params

    @property
    def temperature_range(self) -> Tuple[float, float]:
        return self.min_temperature, self.max_temperature

    @property
    def max_temperature(self) -> float:
        max_data = max(self._params, key=lambda x: x.max_temperature)
        return max_data.max_temperature

    @property
    def min_temperature(self) -> float:
        min_data = min(self._params, key=lambda x: x.min_temperature)
        return min_data.min_temperature

    def heat_capacity(self, temperature: float) -> float:
        """

        Args:
            temperature (float): (K)

        Returns (float): J/(K mol)

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e = \
            params.a, params.b, params.c, params.d, params.e
        return a + b * t + c * t ** 2 + d * t ** 3 + e / (t ** 2)

    def standard_enthalpy(self, temperature: float) -> float:
        """

        Args:
            temperature (float): (K)

        Returns (float): J/mol

        """
        params = self.params(temperature) or None
        t = temperature / 1000
        a, b, c, d, e, f, h = \
            params.a, params.b, params.c, params.d, params.e, \
            params.f, params.h
        enthalpy_kj = a * t + b * t ** 2 / 2 + c * t ** 3 / 3 + \
            d * t ** 4 / 4 - e / t + f - h
        return enthalpy_kj / 1000

    def standard_entropy(self, temperature: float) -> float:
        """

        Args:
            temperature (float): (K)

        Returns (float): J/(mol K)

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e, g = \
            params.a, params.b, params.c, params.d, \
            params.e, params.g
        return a * math.log(t) + b * t + c * t ** 2 / 2 + d * t ** 3 / 3 - \
            e / (2 * t ** 2) + g

    def r_ln_p_p0(self, pressure: float) -> float:
        """
        R ln(p/p0)
        Args:
            pressure (float): Pa

        Returns (float): Rln(p/p0), (J / (mol K))

        """
        return constants.R * math.log(pressure / self._reference_pressure)

    @property
    def enthalpy_zero(self) -> float:
        """

        Returns: Energy at pressure = self._reference_pressure, T=0K (J/mol)

        """
        return self._enthalpy_zero


class InvalidFileError(Exception):
    pass


class Gas(Enum):
    H2 = "H2"
    N2 = "N2"
    O2 = "O2"
    F2 = "F2"
    P2 = "P2"
    P4 = "P4"
    H2O = "H2O"
    NH3 = "NH3"
    NO2 = "NO2"

    def __init__(self, formula: Union[str, Composition]):
        self._formula = str(formula)

        dirname = \
            os.path.join(os.path.abspath(os.path.dirname(__file__)),
                         "molecules",
                         formula)

        # Read prior_info.json
        try:
            path = os.path.join(dirname, "prior_info.json")
            with open(path, "r") as fr:
                self._properties = json.load(fr)
        except FileNotFoundError:
            raise ValueError(f"Data of {formula} molecule is not found.")

        # Read data of thermodynamics function
        try:
            path = os.path.join(dirname, "shomate_nist.dat")
            # TODO: Implement flexibly with abstract class
            self._thermodynamics_function = \
                ShomateThermodynamicsFunction.from_nist_table(path)
        except FileNotFoundError:
            raise ValueError(f"shomate_nist.dat of {formula} "
                             f"molecule is not found.")
        except:
            raise InvalidFileError(f"shomate_nist.dat of {formula} molecule")

        # Read POSCAR
        try:
            path = os.path.join(dirname, "POSCAR")
            self._structure = Poscar.from_file(path).structure
        except FileNotFoundError:
            raise ValueError(f"POSCAR of {formula} molecule is not found.")
        except:
            raise InvalidFileError(f"Failed to read POSCAR of "
                                   f"{formula} molecule")

        # Read Zero point vibrational frequency
        try:
            path = os.path.join(dirname, "fundamental_frequencies.yaml")
            self._fundamental_frequencies = \
                FundamentalFrequencies.from_yaml(path)
        except FileNotFoundError:
            raise ValueError(f"fundamental_frequencies.yaml of "
                             f"{formula} molecule is not found.")
        except:
            raise InvalidFileError(f"Failed to read "
                                   f"fundamental_frequencies.yaml "
                                   f"of {formula} molecule")

    def __str__(self):
        return self.value

    @property
    def composition(self) -> Composition:
        return Composition(str(self))

    @property
    def properties(self) -> Dict[str, Any]:
        """
        Returns:

        """
        return self._properties

    @property
    def structure(self) -> Structure:
        """
        Returns (Structure): Molecule's structure.
        """
        return self._structure

    @property
    def temperature_range(self) -> Tuple[float, float]:
        """
        Returns (Tuple[float, float]): Temperature range that can be applied.
        """
        return self._thermodynamics_function.temperature_range

    @property
    def min_temperature(self) -> float:
        """
        Returns (float): Minimum temperature that can be applied.
        """
        return self._thermodynamics_function.min_temperature

    @property
    def max_temperature(self) -> float:
        """
        Returns (float): Maximum temperature that can be applied.
        """
        return self._thermodynamics_function.max_temperature

    def r_ln_p_p0_joule_per_mol(self, pressure: float) -> float:
        """
        
        Args:
            pressure (float): 

        Returns (float): Rln(p/p0), (J / (mol K))

        """
        return self._thermodynamics_function.r_ln_p_p0(pressure)
    
    def r_ln_p_p0_ev_per_atom(self, pressure: float) -> float:
        """

        Args:
            pressure (float): 

        Returns (float): Rln(p/p0), (eV / atom)

        """
        joule_to_ev = \
            constants.physical_constants["joule-electron volt relationship"][0]
        return joule_to_ev * \
            self.r_ln_p_p0_joule_per_mol(pressure) / \
            (self.n_atom * constants.Avogadro)

    @property
    def formula(self) -> str:
        return self._formula

    @property
    def n_atom(self) -> int:
        return Composition(self._formula).num_atoms

    @property
    def fundamental_frequencies(self) -> FundamentalFrequencies:
        """

        Returns: Frequencies of vibration at zero point (cm^-1)

        """
        return self._fundamental_frequencies

    @property
    def zero_point_vibrational_energy(self) -> float:
        """

        Returns: Energy of vibrational energy at zero point (eV/atom)

        """
        return self._fundamental_frequencies.zero_point_vibrational_energy /\
            self.n_atom

    def heat_capacity(self, temperature) -> float:
        """
        Args:
            temperature (float): Unit is (K)

        Returns (float): Heat capacity C_p (J/mol*K).

        """
        return self._thermodynamics_function.heat_capacity(temperature)

    def standard_enthalpy(self, temperature: float) -> float:
        """
        Standard enthalpy from room temperature.
        Args:
            temperature (float): Unit is (K)

        Returns (float): Standard enthalpy (kJ/mol)
                        from room temperature (T=298.15(K))

        """
        return self._thermodynamics_function.standard_enthalpy(temperature)

    def standard_entropy(self, temperature: float) -> float:
        """
        Args:
            temperature(float): Unit is (K)

        Returns(float): Standard entropy (J/mol*K).

        """
        return self._thermodynamics_function.standard_entropy(temperature)

    def free_energy_shift(self, temperature: float) -> float:
        """
        G(T, P0) - H(T=0K, P0) (eV/mol)

        Returns:

        """
        val_joule =\
            self.standard_enthalpy(temperature) - \
            temperature * self.standard_entropy(temperature) - \
            self._thermodynamics_function.enthalpy_zero
        joule_to_ev = \
            constants.physical_constants["joule-electron volt relationship"][0]
        return val_joule * joule_to_ev / (self.n_atom * constants.Avogadro)

    def energy_shift(self,
                     temperature: float = 0,  
                     pressure: float = 1e+5) -> float:
        """

        Args:
            temperature (float): (K) When 0K, only zero_point_vibrational_energy
                                 is considered.
            pressure (float): (Pa)
                              When specified value is None, 
                              default value is also 1e+5 Pa.

        Returns (float): Free energy shift (eV/atom)

        """
        if abs(temperature) < 1e-12:
            return self.zero_point_vibrational_energy

        return self.zero_point_vibrational_energy + \
            self.free_energy_shift(temperature) + \
            temperature * self.r_ln_p_p0_ev_per_atom(pressure)
