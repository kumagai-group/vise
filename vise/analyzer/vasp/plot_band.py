# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import re
from copy import deepcopy
from typing import Tuple, List
import numpy as np
from pymatgen import Spin

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp import Vasprun
from pymatgen.util.string import latexify
from vise.analyzer.plot_band import BandPlotInfo, BandInfo, XTicks, BandEdge


def greek_to_unicode(label: str) -> str:
    d = {"GAMMA": "Γ", "SIGMA": "Σ", "DELTA": "Δ"}
    for k, v in d.items():
        label = label.replace(k, v)
    return label


def italic_to_roman(label: str) -> str:
    return re.sub(r"([A-Z])_([0-9])", r"{\\rm \1}_\2", label)


class BandPlotInfoFromVasp:
    def __init__(self,
                 vasprun: Vasprun,
                 kpoints_filename: str,
                 vasprun2: Vasprun = None,
                 energy_window: List[float] = None):
        self.vasprun = vasprun
        self.kpoints_filename = kpoints_filename
        self.vasprun2 = vasprun2
        self.energy_window = energy_window

    def make_band_plot_info(self):
        bs = self.vasprun.get_band_structure(self.kpoints_filename,
                                             line_mode=True)
        plot_data = BSPlotter(bs).bs_plot_data(zero_to_efermi=False)
        self._composition = self.vasprun.final_structure.composition

        band_info = [BandInfo(band_energies=self._remove_spin_key(plot_data),
                              band_edge=self._band_edge(bs, plot_data),
                              fermi_level=bs.efermi)]

        if self.vasprun2:
            bs2 = self.vasprun2.get_band_structure(self.kpoints_filename,
                                                   line_mode=True)
            plot_data2 = BSPlotter(bs2).bs_plot_data(zero_to_efermi=False)
            band_info.append(
                BandInfo(band_energies=self._remove_spin_key(plot_data2),
                         band_edge=self._band_edge(bs2, plot_data2),
                         fermi_level=bs.efermi))

        return BandPlotInfo(band_info_set=band_info,
                            distances_by_branch=plot_data["distances"],
                            x_ticks=self._x_ticks(plot_data),
                            title=self._title)

    def _remove_spin_key(self, plot_data):
        """
        Pymatgen at 2020.11.11
         energy: A dict storing bands for spin up and spin down data
            {Spin:[np.array(nb_bands,kpoints),...]} as a list of discontinuous kpath
            of energies. The energy of multiple continuous branches are stored together.

        -> [branch][spin][band][k-point]
        """
        num_spin = len(plot_data["energy"])
        num_branch = len(plot_data["energy"]["1"])

        result = [[[] for _ in range(num_spin)] for __ in range(num_branch)]
        for spin_idx, (_, branch_energies) in enumerate(
                sorted(plot_data["energy"].items(),
                       key=lambda item: item[0],
                       reverse=True)):
            for branch_idx, branch_energy in enumerate(branch_energies):
                if self.energy_window:
                    removed_idxs = []
                    for i in range(len(branch_energy)):
                        _max = np.max(branch_energy[i, :])
                        _min = np.min(branch_energy[i, :])
                        if not self.in_energy(_max, _min):
                            removed_idxs.append(i)
                    x = np.delete(branch_energy, removed_idxs, axis=0).tolist()
                else:
                    x = branch_energy.tolist()

                result[branch_idx][spin_idx] = deepcopy(x)

        return result

    def in_energy(self, _max, _min):
        return _max >= self.energy_window[0] and _min <= self.energy_window[1]

    def _x_ticks(self, plot_data):
        labels = self._sanitize_labels(plot_data["ticks"]["label"])
        distances = plot_data["ticks"]["distance"]
        return XTicks(labels=labels, distances=distances)

    def _band_edge(self, bs, plot_data):
        if bs.is_metal():
            return None
        else:
            return BandEdge(
                vbm=plot_data["vbm"][0][1],
                cbm=plot_data["cbm"][0][1],
                vbm_distances=[i[0] for i in plot_data["vbm"]],
                cbm_distances=[i[0] for i in plot_data["cbm"]])

    @property
    def _title(self):
        return latexify(self._composition.reduced_formula)

    @staticmethod
    def _sanitize_labels(labels):
        return [italic_to_roman(greek_to_unicode(label)) for label in labels]
