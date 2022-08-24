# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List, Optional

import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Vasprun

from vise.analyzer.dos_data import DosData
from vise.analyzer.dos_data import PDos
from vise.util.logger import get_logger
logger = get_logger(__name__)


class DosDataFromVasp:

    def __init__(self, vasprun: Vasprun,
                 vertical_lines: Optional[List[float]] = None,
                 base_energy: float = 0.0,
                 crop_first_value=False,
                 energy_window: List[float] = None):
        self.complete_dos = vasprun.complete_dos
        self.vertical_lines = vertical_lines or []
        self.base_energy = base_energy
        self.energy_window = energy_window
        self.min_energy_idx = 1 if crop_first_value else 0
        self.max_energy_idx = len(self.complete_dos.energies)

    def make_dos_data(self):
        energies = self.complete_dos.energies.tolist()
        if self.energy_window:
            _min = self.get_min_energy_idx(energies)
            _max = self.get_max_energy_idx(energies)
            self.min_energy_idx = _min if _min else self.min_energy_idx
            self.max_energy_idx = _max if _max else self.max_energy_idx

        energies = energies[self.min_energy_idx: self.max_energy_idx + 1]

        return DosData(energies=energies,
                       total=np.array(self._total),
                       pdos=self._pdos,
                       vertical_lines=self.vertical_lines,
                       base_energy=self.base_energy)

    def get_max_energy_idx(self, energies):
        for i, e in enumerate(energies):
            if e > self.energy_window[1]:
                return i - 1

    def get_min_energy_idx(self, energies):
        for i, e in enumerate(energies):
            if e >= self.energy_window[0]:
                return i

    @property
    def _pdos(self):
        result = []
        for dos_by_site in self.complete_dos.pdos.values():
            pdos_kwargs = {}
            for orbital, dos_by_orbital in dos_by_site.items():
                pdos = np.array(
                    [dos_by_orbital[s] for s in [Spin.up, Spin.down]
                     if s in dos_by_orbital])

                pdos_kwargs[str(orbital)] = \
                    pdos[:, self.min_energy_idx:self.max_energy_idx + 1]
            try:
                result.append(PDos.from_dict(pdos_kwargs))
            except TypeError:
                logger.warning("Orbital DOSes are required. SET LORBIT = 11.")
                raise

        return result

    @property
    def _total(self):
        result = []
        for s in [Spin.up, Spin.down]:
            if s in self.complete_dos.densities:
                dos = self.complete_dos.densities[s]
                result.append(dos[self.min_energy_idx:self.max_energy_idx + 1])
        return result
