# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import List, Dict
import numpy as np
import matplotlib.pyplot as plt

from pymatgen import Spin


class DosPlotter:
    def __init__(self,
                 energies: List[float],
                 doses: List[Dict[str, Dict[Spin, List[float]]]]):
        self.energies = energies
        self.doses = doses
        self.plt = plt

    def construct_plot(self):

        _, axs = self.plt.subplots(2, 1, sharex=True)

        for i, dos in enumerate(self.doses):
            ax = axs[i]
            for name, dos_each_name in dos.items():
                for spin, dos_each_spin in dos_each_name.items():
                    dos_for_plot = [d * int(spin) for d in dos_each_spin]
                    ax.plot(self.energies, dos_for_plot)
#            self.plt.legend()

            ax.axhline(0)

