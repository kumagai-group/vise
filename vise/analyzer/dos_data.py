# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy, deepcopy
from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List

import numpy as np
from monty.json import MSONable
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn


logger = get_logger(__name__)


@dataclass
class PDos(MSONable):
    s:   np.ndarray  # [by spin][by energy]
    px:  np.ndarray
    py:  np.ndarray
    pz:  np.ndarray
    dxy: np.ndarray
    dyz: np.ndarray
    dxz: np.ndarray
    dx2: np.ndarray
    dz2: np.ndarray
    f_3: Optional[np.ndarray] = None
    f_2: Optional[np.ndarray] = None
    f_1: Optional[np.ndarray] = None
    f0:  Optional[np.ndarray] = None
    f1:  Optional[np.ndarray] = None
    f2:  Optional[np.ndarray] = None
    f3:  Optional[np.ndarray] = None

    @property
    def p(self):
        return self.px + self.py + self.pz

    @property
    def d(self):
        return self.dxy + self.dyz + self.dxz + self.dx2 + self.dz2

    @classmethod
    def from_dict(cls, d):
        if "px" in d:
            return super().from_dict(d)

        logger.warning("Because the PDOS are not decomposed to each orbital,"
                       "e.g., px, py, pz, they are equally divided.")
        x = dict(s=d["s"], px=d["p"] / 3, py=d["p"] / 3, pz=d["p"] / 3,
                 dxy=d["d"] / 5, dyz=d["d"] / 5, dxz=d["d"] / 5, dx2=d["d"] / 5,
                 dz2=d["d"] / 5)
        if "f" in d:
            x.update({
                "f_3": d["f"] / 7,
                "f_2": d["f"] / 7,
                "f_1": d["f"] / 7,
                "f0": d["f"] / 7,
                "f1": d["f"] / 7,
                "f2": d["f"] / 7,
                "f3": d["f"] / 7})
        return cls(**x)

    @property
    def f(self):
        try:
            return (self.f_3 + self.f_2 + self.f_1 +
                    self.f0 + self.f1 + self.f2 + self.f3)
        except TypeError:
            return None

    def __add__(self, other: "PDos"):
        args = {}
        for k, v in self.__dict__.items():
            if v is None:
                args[k] = None
            else:
                args[k] = v + getattr(other, k)
        return PDos(**args)


@dataclass
class DosData(MSONable):
    energies: List[float]
    total: np.ndarray
    pdos: List[PDos]
    vertical_lines: List[float]
    base_energy: Optional[float] = 0.0

    @property
    def spin(self):
        return False if len(self.total) == 1 else True

    def dos_plot_data(self,
                      grouped_atom_indices: Dict[str, List[int]],
                      energy_range: Optional[List[float]] = None,
                      dos_ranges: Optional[List[List[float]]] = None,
                      ) -> "DosPlotData":
        """
        Args:
            grouped_atom_indices:
                key: name of the grouped
            energy_range:
            dos_ranges:

        Returns:

        """

        if dos_ranges is not None:
            assert len(grouped_atom_indices) + 1 == len(dos_ranges)  # total+pdos

        # Total dos does not have spin decomposition.
        doses = [[DosBySpinEnergy("", self.total.tolist())]]
        names = ["total"]

        for name, atom_indices_by_group in grouped_atom_indices.items():
            pdos_list = [self.pdos[idx] for idx in atom_indices_by_group]
            pdos = reduce(lambda x, y: x + y, pdos_list)
            pdos_by_ax = [DosBySpinEnergy("s", pdos.s.tolist()),
                          DosBySpinEnergy("p", pdos.p.tolist()),
                          DosBySpinEnergy("d", pdos.d.tolist())]
            if pdos.f is not None:
                pdos_by_ax.append(DosBySpinEnergy("f", pdos.f.tolist()))

            doses.append(pdos_by_ax)
            names.append(name)

        energy_range = energy_range or [-5, 10]
        abs_xlim = [x + self.base_energy for x in energy_range]

        if dos_ranges is None:
            dos_ranges = default_dos_ranges(abs_xlim, doses, self.energies,
                                            self.spin)

        energies = [e - self.base_energy for e in self.energies]
        shifted_energy_lines = [e - self.base_energy
                                for e in self.vertical_lines]

        return DosPlotData(energies, doses, names, energy_range, dos_ranges,
                           shifted_energy_lines)


def default_dos_ranges(energy_range, doses, energies, spin_polarized,
                       multi: float = 1.1, round_digit: int = 2):
    mask = np.ma.masked_outside(energies, energy_range[0], energy_range[1]).mask
    max_dos_by_ax = []
    for dos_by_ax in doses:
        max_dos_by_ax.append(np.max([dos_by_spin.max_dos(mask)
                                     for dos_by_spin in dos_by_ax]))

    total_dos_max = max_dos_by_ax[0]
    plot_maxes = [total_dos_max]
    if len(max_dos_by_ax) > 1:
        pdos_max = max(max_dos_by_ax[1:])
        plot_maxes.extend([pdos_max] * (len(doses) - 1))

    max_y_ranges = [round(i * multi, round_digit) for i in plot_maxes]

    if spin_polarized:
        return [[-y, y] for y in max_y_ranges]
    return [[0, y] for y in max_y_ranges]


@dataclass
class DosBySpinEnergy(MSONable):
    name: str  # e.g., "s", "p", ...
    dos: List[List[float]]  # [by spin][by energy]

    def max_dos(self, mask: List[bool] = None):
        return max([np.max(np.ma.masked_array(d, mask).compressed())
                    for d in self.dos])


@dataclass
class DosPlotData(MSONable, ToJsonFileMixIn):
    relative_energies: List[float]
    doses: List[List[DosBySpinEnergy]]  # [by ax][by orbital, i.e., s, p, d, f]
    names: List[str]  # For each ax, e.g., ["total", "H"]
    energy_range: List[float]  # e.g., [-5.0, 6.0]
    dos_ranges: List[List[float]]  # e.g., [[0.0, 6.0], [-1.0, 5.0], ...]
    energy_lines: List[float]  # e.g., VBM and CBM or Fermi level

    @classmethod
    def from_dict(cls, d):
        for k in copy(d):
            if k[0] == "@":
                d.pop(k)
        # For backward compatibility
        if "xlim" in d:
            d["energy_range"] = d.pop("xlim")
        if "ylim_set" in d:
            d["dos_ranges"] = d.pop("ylim_set")
        if "vertical_lines" in d:
            d["energy_lines"] = d.pop("vertical_lines")
        for i, x in enumerate(d["doses"]):
            for j, y in enumerate(x):
                d["doses"][i][j] = DosBySpinEnergy.from_dict(y)

        if "energy_range" not in d:
            d["energy_range"] = [-5, 10]
        if "dos_ranges" not in d:
            spin_polarized = len(d["doses"][0][0].dos) == 2
            d["dos_ranges"] = default_dos_ranges(
                d["energy_range"], d["doses"], d["relative_energies"],
                spin_polarized, multi=1.1, round_digit=2)

        return cls(**d)


def scissor_energy(dos_plot_data: DosPlotData,
                   energy_shift: float) -> DosPlotData:
    result = deepcopy(dos_plot_data)
    gap_middle = np.mean(dos_plot_data.energy_lines)
    energy_over_gap_middle = np.array(result.relative_energies) >= gap_middle
    shifted_energy_start_idx = int(np.argmax(energy_over_gap_middle))
    for i in range(shifted_energy_start_idx, len(result.relative_energies)):
        result.relative_energies[i] += energy_shift

    result.energy_lines[1] += energy_shift
    return result
