# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from typing import Union, List, Tuple, Dict, Any
from enum import Enum
import json
import math
import os
from pathlib import Path
import re
import yaml

from scipy.constants import speed_of_light, physical_constants, R, Avogadro

from monty.json import MSONable

from pymatgen.core.composition import Composition
from pymatgen.io.vasp.inputs import Poscar

from vise.config import REFERENCE_PRESSURE
from vise.util.error_classes import InvalidFileError
from vise.input_set.prior_info import PriorInfo

""" Data is from NIST Chemistry WebBook, SRD 69 """


class FundamentalFrequency(MSONable):

    def __init__(self, wave_number: float, degeneration: int):
        """

        Args:
            wave_number (float): frequency (cm^-1)
            degeneration (int): number of degeneration
        """
        self.wave_number = wave_number
        self.degeneration = degeneration

    def __str__(self):
        return f"wave_number: {self.wave_number} /cm, " \
               f"degeneration: {self.degeneration}"

    @property
    def zero_point_vibrational_energy(self) -> float:
        h = physical_constants["Planck constant in eV s"][0]
        sol = speed_of_light
        # h/2 * degeneration * wave_num(cm^-1) * 100(cm^-1->m^-1) * light speed
        return 0.5 * h * self.degeneration * self.wave_number * 100 * sol


class FundamentalFrequencies:

    def __init__(self,
                 frequencies: List[FundamentalFrequency],
                 reference: str):
        """
        Args:
            frequencies (list of FundamentalFrequency):
                Fundamental frequencies.
            reference (str):
                Citation of the data.
        """
        self.frequencies = frequencies
        self.reference = reference

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
        return sum(f.zero_point_vibrational_energy for f in self.frequencies)


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
        self.temperature_range = temperature_range
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.g = g
        self.h = h

    def can_apply(self, temp: float) -> bool:
        return self.temperature_range[0] <= temp <= self.temperature_range[1]

    @property
    def min_temperature(self) -> float:
        return self.temperature_range[0]

    @property
    def max_temperature(self) -> float:
        return self.temperature_range[1]

    def __str__(self):
        outs = [f"Temperature_range: {self.temperature_range}\n",
                f"A: {self.a}\n", f"B: {self.b}\n", f"C: {self.c}\n",
                f"D: {self.d}\n", f"E: {self.e}\n", f"F: {self.f}\n",
                f"G: {self.g}\n", f"H: {self.h}\n"]
        return "\n".join(outs)


class ShomateThermodynamicsFunction(AbstractThermodynamicsFunction):
    """Shomate thermodynamics Function. """

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
        self.enthalpy_zero = enthalpy_zero
        self.reference_pressure = reference_pressure

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
            ShomateThermodynamicsFunction object
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
        """Return parameters of function at temperature like {a:1.1, b:2.2, ...}

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
        return max_data.temperature_range[1]

    @property
    def min_temperature(self) -> float:
        min_data = min(self._params, key=lambda x: x.min_temperature)
        return min_data.temperature_range[0]

    def heat_capacity(self, temperature: float) -> float:
        """

        Args:
            temperature (float): (K)

        Returns (float): J/(K mol)

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e = params.a, params.b, params.c, params.d, params.e
        return a + b * t + c * t ** 2 + d * t ** 3 + e / (t ** 2)

    def standard_enthalpy(self, temperature: float) -> float:
        """

        Args:
            temperature (float): (K)

        Returns (float): J/mol

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e, f, h = \
            params.a, params.b, params.c, params.d, params.e, params.f, params.h
        enthalpy_kj = (a * t + b * t ** 2 / 2 + c * t ** 3 / 3
                       + d * t ** 4 / 4 - e / t + f - h)
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
            params.a, params.b, params.c, params.d, params.e, params.g
        return (a * math.log(t) + b * t + c * t ** 2 / 2 + d * t ** 3 / 3
                - e / (2 * t ** 2) + g)

    def r_ln_p_p0(self, pressure: float) -> float:
        """ R ln(p/p0)
        Args:
            pressure (float): Pa

        Returns (float): Rln(p/p0), (J / (mol K))

        """
        return R * math.log(pressure / self.reference_pressure)


class Gas(Enum):
    F2 = "F2"
    H2 = "H2"
    H2O = "H2O"
    N2 = "N2"
    NH3 = "NH3"
    NO2 = "NO2"
    O2 = "O2"
    P4 = "P4"

    def __str__(self):
        return self.name

    def __init__(self, formula: Union[str, Composition]):
        self.formula = str(formula)
        dirname = Path(__file__).parent / "molecules" / self.formula

        try:
            self.properties = PriorInfo.load_json(dirname / "prior_info.yaml")
        except:
            pass
        self.thermodynamics_function = ShomateThermodynamicsFunction.\
            from_nist_table(dirname / "shomate_nist.dat")
        self.fundamental_frequencies = FundamentalFrequencies.\
            from_yaml(dirname / "fundamental_frequencies.yaml")

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])

    @property
    def temperature_range(self) -> Tuple[float, float]:
        """Temperature range that can be applied. """
        return self.thermodynamics_function.temperature_range

    @property
    def min_temperature(self) -> float:
        """Minimum temperature that can be applied. """
        return self.thermodynamics_function.min_temperature

    @property
    def max_temperature(self) -> float:
        """Maximum temperature that can be applied. """
        return self.thermodynamics_function.max_temperature

    def r_ln_p_p0_j_mol(self, pressure: float) -> float:
        """Rln(p/p0), (J / (mol K)) """
        return self.thermodynamics_function.r_ln_p_p0(pressure)
    
    def r_ln_p_p0_ev_per_atom(self, pressure: float) -> float:
        """Rln(p/p0), (eV / atom) """
        j2ev = physical_constants["joule-electron volt relationship"][0]
        return j2ev * self.r_ln_p_p0_j_mol(pressure) / (self.n_atom * Avogadro)

    @property
    def n_atom(self) -> int:
        return Composition(self.formula).num_atoms

    @property
    def zero_point_vibrational_energy(self) -> float:
        """Energy of vibrational energy at zero point (eV/atom) """
        zpve = self.fundamental_frequencies.zero_point_vibrational_energy
        return zpve / self.n_atom

    def heat_capacity(self, temperature: float) -> float:
        """Returns (float): Heat capacity C_p (J/mol*K). """
        return self.thermodynamics_function.heat_capacity(temperature)

    def standard_enthalpy(self, temperature: float) -> float:
        """Standard enthalpy (kJ/mol) from room temperature (T=298.15(K)) """
        return self.thermodynamics_function.standard_enthalpy(temperature)

    def standard_entropy(self, temperature: float) -> float:
        """Standard entropy (J/mol*K). """
        return self.thermodynamics_function.standard_entropy(temperature)

    def vib_rot_term(self, temperature: float) -> float:
        """G(T, P0) - H(T=0K, P0) (eV/mol) """
        h = self.standard_enthalpy(temperature)
        s = self.standard_entropy(temperature)
        h_0 = self.thermodynamics_function.enthalpy_zero
        val_joule = h - temperature * s - h_0
        joule_to_ev = physical_constants["joule-electron volt relationship"][0]
        return val_joule * joule_to_ev / (self.n_atom * Avogadro)

    def free_e_shift(self,
                     temperature: float = 0,
                     pressure: float = 1e+5) -> float:
        """Free energy shift (eV/atom) """
        if abs(temperature) < 1e-12:
            return 0.0

        vib_rot_term = self.vib_rot_term(temperature)
        trans_term = temperature * self.r_ln_p_p0_ev_per_atom(pressure)
        return vib_rot_term + trans_term

    def energy_shift(self,
                     temperature: float = 0,  
                     pressure: float = 1e+5) -> float:
        """Free energy shift including zero point vibration (eV/atom) """
        if abs(temperature) < 1e-12:
            return self.zero_point_vibrational_energy

        zero_point_vib = self.zero_point_vibrational_energy
        free_e_shift = self.free_e_shift(temperature, pressure)
        return zero_point_vib + free_e_shift
