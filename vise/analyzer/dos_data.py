# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List

import numpy as np


class PDos:
    def __init__(self,
                 s:   np.ndarray,  # [by spin][by energy]
                 px:  np.ndarray,
                 py:  np.ndarray,
                 pz:  np.ndarray,
                 dxy: np.ndarray,
                 dyz: np.ndarray,
                 dxz: np.ndarray,
                 dx2: np.ndarray,
                 dz2: np.ndarray,
                 f_3: Optional[np.ndarray] = None,
                 f_2: Optional[np.ndarray] = None,
                 f_1: Optional[np.ndarray] = None,
                 f0:  Optional[np.ndarray] = None,
                 f1:  Optional[np.ndarray] = None,
                 f2:  Optional[np.ndarray] = None,
                 f3:  Optional[np.ndarray] = None,
                 ):
        self.s = s
        self.px = px
        self.py = py
        self.pz = pz
        self.dxy = dxy
        self.dyz = dyz
        self.dxz = dxz
        self.dx2 = dx2
        self.dz2 = dz2
        self.f_3 = f_3
        self.f_2 = f_2
        self.f_1 = f_1
        self.f0 = f0
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3

    @property
    def p(self):
        return self.px + self.py + self.pz

    @property
    def d(self):
        return self.dxy + self.dyz + self.dxz + self.dx2 + self.dz2

    @property
    def f(self):
        try:
            return (self.f_3 + self.f_2 + self.f_1 +
                    self.f0 + self.f1 + self.f2 + self.f3)
        except TypeError:
            return None

    def __add__(self, other: "PDos"):
        args = {k: v + getattr(other, k) for k, v in self.__dict__.items()}
        return PDos(**args)


class DosData:
    def __init__(self,
                 energies: List[float],
                 total: np.ndarray,
                 pdos: List[PDos]):

        self.energies = energies
        self.total = total
        self.pdos = pdos

    def dos_plot_data(self,
                      grouped_atom_indices: Dict[str, List[int]],
                      vertical_lines: List[float],
                      base_energy: Optional[float] = 0.0,
                      xlim: Optional[List[float]] = None,
                      ylim_set: Optional[List[List[float]]] = None,
                      ) -> "DosPlotData":

        if ylim_set is not None:
            assert len(grouped_atom_indices) + 1 == len(ylim_set)  # total+pdos

        doses = [[DosBySpinEnergy("total", self.total)]]
        for name, atom_indices_by_group in grouped_atom_indices.items():
            pdos_list = [self.pdos[idx] for idx in atom_indices_by_group]
            pdos = reduce(lambda x, y: x + y, pdos_list)
            pdos_by_ax = [DosBySpinEnergy(f"{name}-s", pdos.s),
                          DosBySpinEnergy(f"{name}-p", pdos.p),
                          DosBySpinEnergy(f"{name}-d", pdos.d)]
            if pdos.f is not None:
                pdos_by_ax.append(DosBySpinEnergy(f"{name}-f", pdos.f))

            doses.append(pdos_by_ax)

        xlim = xlim or [-10, 10]
        ylim_set = ylim_set or [[-y, y] for y in self.max_y_ranges(doses)]

        energies = [e - base_energy for e in self.energies]
        shifted_vertical_lines = [e - base_energy for e in vertical_lines]

        return DosPlotData(energies, doses, xlim, ylim_set, shifted_vertical_lines)

    @staticmethod
    def max_y_ranges(doses, multi=1.1, round_digit=2):
        max_modulus_dos_by_ax = []
        for dos_by_ax in doses:
            max_dos_by_ax = np.max([dos_by_spin.max_dos()
                                    for dos_by_spin in dos_by_ax])
            max_modulus_dos_by_ax.append(max_dos_by_ax)
        return [round(i * multi, round_digit) for i in max_modulus_dos_by_ax]


@dataclass
class DosBySpinEnergy:
    name: str
    dos: np.array  # [by spin][by energy]

    def max_dos(self):
        return np.max(self.dos)


@dataclass()
class DosPlotData:
    relative_energies: List[float]
    doses: List[List[DosBySpinEnergy]]  # [by ax][by orbital]
    xlim: List[float]
    ylim_set: List[List[float]]
    vertical_lines: List[float]
