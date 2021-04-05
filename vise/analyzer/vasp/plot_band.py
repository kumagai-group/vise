# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import re
from copy import deepcopy
from typing import List
import numpy as np

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp import Vasprun
from vise.analyzer.plot_band import BandPlotInfo, BandInfo, XTicks, BandEdge
from vise.analyzer.plot_brillouin_zone import BZPlotInfo
from vise.util.string import latexify


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
        self.bs = self.vasprun.get_band_structure(self.kpoints_filename,
                                                  line_mode=True)

    def make_band_plot_info(self):
        bs_plotter = BSPlotter(self.bs)
        plot_data = bs_plotter.bs_plot_data(zero_to_efermi=False)
        distances = [list(d) for d in plot_data["distances"]]
        self._composition = self.vasprun.final_structure.composition

        band_info = [BandInfo(band_energies=self._remove_spin_key(plot_data),
                              band_edge=self._band_edge(self.bs, plot_data),
                              fermi_level=self.bs.efermi)]

        if self.vasprun2:
            bs2 = self.vasprun2.get_band_structure(self.kpoints_filename,
                                                   line_mode=True)
            plot_data2 = BSPlotter(bs2).bs_plot_data(zero_to_efermi=False)
            band_info.append(
                BandInfo(band_energies=self._remove_spin_key(plot_data2),
                         band_edge=self._band_edge(bs2, plot_data2),
                         fermi_level=self.bs.efermi))

        x = bs_plotter.get_ticks_old()
        x_ticks = XTicks(_sanitize_labels(x["label"]), x["distance"])

        return BandPlotInfo(band_info_set=band_info,
                            distances_by_branch=distances,
                            x_ticks=x_ticks,
                            title=self._title)

    def make_bz_plot_info(self):
        rec_lat = self.vasprun.final_structure.lattice.reciprocal_lattice
        faces = [[[float(k) for k in j] for j in i] for i in rec_lat.get_wigner_seitz_cell()]
        labels = {}

        concat = False
        band_paths = []
        init_point = None
        for kpoint in self.bs.kpoints:
            if kpoint.label:
                c_coords = list(kpoint.cart_coords)
                f_coords = list(kpoint.frac_coords)
                label = greek_to_unicode(kpoint.label)
                labels[label] = {"cart": c_coords, "frac": f_coords}
                if concat is False and init_point:
                    band_paths.append([init_point, c_coords])
                init_point = c_coords
                concat = True
            else:
                concat = False

        return BZPlotInfo(faces, labels, band_paths, rec_lat.matrix.tolist())

    def _remove_spin_key(self, plot_data) -> List[List[List[List[float]]]]:
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


def _sanitize_label(label):
    return italic_to_roman(greek_to_unicode(label))


def _sanitize_labels(labels):
    return [_sanitize_label(label) for label in labels]
