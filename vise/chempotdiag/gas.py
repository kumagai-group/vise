# -*- coding: utf-8 -*-
from functools import reduce
from math import pi, log, exp, sqrt
from pathlib import Path
from typing import Union

import yaml
from pymatgen.core.composition import Composition
from scipy.constants import hbar, eV, k, m_u

from vise.util.logger import get_logger
from vise.config import REFERENCE_PRESSURE, ROOM_TEMPERATURE
""" 
Characteristic data is from
McQuarrie and Simon, Physical Chemistry A molecular approach,
"""

logger = get_logger(__name__)
parent = Path(__file__).parent

with open(parent / "molecules" / "molecule_data.yaml", 'r') as fr:
    MOLECULE_DATA = yaml.safe_load(fr)


class Gas:
    def __init__(self,
                 formula: Union[str, Composition],
                 temperature: float = ROOM_TEMPERATURE,
                 pressure: float = REFERENCE_PRESSURE):
        """
        Args:
            formula (list):
            temperature (float): enthalpy (J/mol) at T=0K,
                                   p=reference_pressure
            pressure (float): reference_pressure (Pa)
        """
        self.formula = str(formula)
        if formula not in MOLECULE_DATA:
            raise NotRegisteredError(
                f"{self.formula} is not registered as a molecule in vise."  
                f"Candidates are {', '.join(MOLECULE_DATA.keys())}.")

        self.char_vib_temp = MOLECULE_DATA[formula]["vib"]
        self.char_rot_temp = MOLECULE_DATA[formula]["rot"]
        self.mag_degeneracy = MOLECULE_DATA[formula]["mag_degeneracy"]
        self.sym_num = MOLECULE_DATA[formula]["sym_num"]
        if len(self.char_rot_temp) == 1:
            self.is_linear = True
        else:
            self.is_linear = False
        self.t = temperature
        self.pressure = pressure

    @property
    def zero_point_vibrational_energy(self) -> float:
        return sum(t for t in self.char_vib_temp) * k / eV / 2  # eV

    @property
    def kt_in_ev(self) -> float:
        return self.t * k / eV  # [eV]

    @property
    def kt(self) -> float:
        return self.t * k  # [J]

    @property
    def trans_free_energy(self) -> float:
        weight = Composition(self.formula).weight * m_u  # kg
        quantum_volume = (2 * pi * hbar ** 2 / weight / self.kt) ** 1.5  # m^3
        dist_func = self.kt / (self.pressure * quantum_volume)  # -
        return - self.kt_in_ev * log(dist_func)  # eV

    def rot_free_energy(self, high_t_limit_factor: float = 20.0) -> float:
        if self.is_linear:
            # high-temperature limit
            if self.char_rot_temp[0] * high_t_limit_factor < self.t:
                q = self.t / self.char_rot_temp[0]
            else:
                factor = self.char_rot_temp[0] / self.t
                q = sum([(2 * j + 1) * exp(-j * (j + 1) * factor)
                         for j in range(0, int(high_t_limit_factor * 5))])
        else:
            if max(self.char_rot_temp) * high_t_limit_factor > self.t:
                raise ValueError(
                    "Now, rotational free energy for non-linear molecule is "
                    "estimated only at the high-temperature limit. Set more"
                    f"than {max(self.char_rot_temp) * high_t_limit_factor }K")
            char_temp_mul = reduce(lambda x, y: x*y, self.char_rot_temp)
            q = sqrt(pi * self.t ** 3) / sqrt(char_temp_mul)

        return - self.kt_in_ev * log(q / self.sym_num)

    @property
    def vib_free_energy(self) -> float:
        logs = [log(1 / (1 - exp(-t / self.t))) for t in self.char_vib_temp]
        return - self.kt_in_ev * sum(logs)

    @property
    def spin_free_energy(self) -> float:
        return - self.kt_in_ev * log(self.mag_degeneracy)

    @property
    def n_atoms(self) -> int:
        return Composition(self.formula).num_atoms

    @property
    def chem_pot_shift(self) -> float:
        """Free energy shift (eV/atom) """
        return (self.trans_free_energy + self.rot_free_energy() +
                self.vib_free_energy + self.spin_free_energy) / self.n_atoms

    @property
    def energy_shift(self) -> float:
        """Free energy shift including zero point vibration (eV/atom) """
        return (self.chem_pot_shift +
                self.zero_point_vibrational_energy / self.n_atoms)


class NotRegisteredError(Exception):
    pass

