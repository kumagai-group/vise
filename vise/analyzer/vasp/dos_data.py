# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from pymatgen import Spin
from pymatgen.io.vasp import Vasprun

from vise.analyzer.dos_data import DosData
from vise.analyzer.dos_data import PDos


class VaspDosData(DosData):

    def __init__(self, vasprun: Vasprun, crop_first_value=False):
        self.complete_dos = vasprun.complete_dos
        self.crop_first_value = crop_first_value
        energies = self.complete_dos.energies.tolist()
        if self.crop_first_value:
            energies = energies[1:]
        super().__init__(energies=energies,
                         total=np.array(self._total),
                         pdos=self._pdos)

    @property
    def _pdos(self):
        result = []
        for dos_by_site in self.complete_dos.pdos.values():
            pdos_kwargs = {}
            for orbital, dos_by_orbital in dos_by_site.items():
                if self.crop_first_value:
                    pdos = [dos_by_orbital[s][1:] for s in [Spin.up, Spin.down]
                            if s in dos_by_orbital]
                else:
                    pdos = [dos_by_orbital[s] for s in [Spin.up, Spin.down]
                            if s in dos_by_orbital]

                pdos_kwargs[str(orbital)] = np.array(pdos)
            result.append(PDos(**pdos_kwargs))
        return result

    @property
    def _total(self):
        result = []
        for s in [Spin.up, Spin.down]:
            if s in self.complete_dos.densities:
                dos = np.copy(self.complete_dos.densities[s])
                if self.crop_first_value:
                    dos = dos[1:]
                result.append(dos)
        return result
