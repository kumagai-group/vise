# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from itertools import tee
import re
from typing import List, Dict, Optional
import matplotlib.pyplot as plt
from dataclasses import dataclass

from pymatgen import Spin
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp import Vasprun

from vise.util.matplotlib import float_to_int_formatter


@dataclass()
class XTicks:
    # labels include "x$\\mid$y"
    labels: List[str]
    distances: List[float]


@dataclass
class BandEdge:
    vbm: float
    cbm: float
    vbm_distances: List[float]
    cbm_distances: List[float]


@dataclass
class BandInfo:
    band_energies_by_branch: List[Dict[Spin, List[List[float]]]]
    band_edge: Optional[BandEdge] = None
    fermi_level: Optional[float] = None


class BandPlotterDefaults:
    def __init__(self,
                 band_linestyle: Optional[List[str]] = None,
                 band_colors: Optional[List[str]] = None,
                 band_linewidth: float = 1.0,
                 band_edge_circle_size: int = 100
                 ):
        self.band_linestyle = band_linestyle or ["-"]
        self.band_colors = band_colors or ["red", "blue"]
        self.band_linewidth = band_linewidth
        self.band_edge_circle_size = band_edge_circle_size

    @property
    def band_args(self):
        return {"linestyle": self.band_linestyle[0],
                "color": self.band_colors[0],
                "linewidth": self.band_linewidth}

    @property
    def band_edge_circles_args(self):
        return {"marker": "o",
                "s": self.band_edge_circle_size}


class BandPlotter:

    def __init__(self,
                 band_info_set: List[BandInfo],
                 distances_by_branch: List[List[float]],
                 x_ticks: XTicks,
                 y_range: List[float],
                 defaults: Optional[BandPlotterDefaults] = BandPlotterDefaults()
                 ):

        assert distances_by_branch[0][0] == x_ticks.distances[0]
        assert distances_by_branch[-1][-1] == x_ticks.distances[-1]

        self.band_info_set = band_info_set
        self.distances_by_branch = distances_by_branch
        self.x_ticks = x_ticks
        self.y_range = y_range
        self.defaults = defaults
        self.plt = plt

    def construct_plot(self):
        self._set_band_set()
        self._set_x_range()
        self._set_y_range()
        self._set_axes()
        self._set_xticks()
        self._set_formatter()
        self.plt.tight_layout()

    def _set_band_set(self):
        for band_info in self.band_info_set:
            self._set_band_structures(band_info.band_energies_by_branch)
            if band_info.band_edge:
                self._set_band_edge(band_info.band_edge)
            if band_info.fermi_level:
                self._set_fermi_level(band_info.fermi_level)

    def _set_band_structures(self, energies_by_branch):
        for distances_of_a_branch, energies_of_a_branch \
                in zip(self.distances_by_branch, energies_by_branch):

            for spin, energies_of_a_spin in energies_of_a_branch.items():
                for energies_of_a_band in energies_of_a_spin:
                    self.plt.plot(distances_of_a_branch, energies_of_a_band,
                                  **self.defaults.band_args)

    def _set_band_edge(self, band_edge):
        self.plt.axhline(y=band_edge.vbm)
        for dist in band_edge.vbm_distances:
            self.plt.scatter(dist, band_edge.vbm,
                             **self.defaults.band_edge_circles_args)

        self.plt.axhline(y=band_edge.cbm)
        for dist in band_edge.cbm_distances:
            self.plt.scatter(dist, band_edge.cbm,
                             **self.defaults.band_edge_circles_args)

    def _set_fermi_level(self, fermi_level):
        self.plt.axhline(y=fermi_level)

    def _set_x_range(self):
        self.plt.xlim(self.distances_by_branch[0][0],
                      self.distances_by_branch[-1][-1])

    def _set_y_range(self):
        self.plt.ylim(self.y_range[0], self.y_range[1])

    def _set_axes(self):
        self.plt.xlabel("Wave vector")
        self.plt.ylabel("Energy (eV)")

    def _set_xticks(self):
        axis = self.plt.gca()
        axis.set_xticks(self.x_ticks.distances)
        axis.set_xticklabels(self.x_ticks.labels)
        for distance, label in zip(self.x_ticks.distances[1:-1],
                                   self.x_ticks.labels[1:-1]):
            linestyle = "-" if "\\mid" in label else "--"
            plt.axvline(x=distance, linestyle=linestyle)

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.yaxis.set_major_formatter(float_to_int_formatter)


def greek_to_unicode(label: str) -> str:
    d = {"GAMMA": "Γ", "SIGMA": "Σ", "DELTA": "Δ"}
    for k, v in d.items():
        label = label.replace(k, v)
    return label


def italic_to_roman(label: str) -> str:
    return re.sub(r"([A-Z])_([0-9])", r"{\\rm \1}_\2", label)


class VaspBandPlotter(BandPlotter):
    def __init__(self, vasprun: Vasprun, kpoints_filename: str):

        self.bs = vasprun.get_band_structure(kpoints_filename, line_mode=True)
        self.plot_data = BSPlotter(self.bs).bs_plot_data(zero_to_efermi=False)

        band_info = BandInfo(band_energies_by_branch=self.plot_data["energy"],
                             band_edge=self._band_edge,
                             fermi_level=self.bs.efermi)

        super().__init__(band_info_set=[band_info],
                         distances_by_branch=self.plot_data["distances"],
                         x_ticks=self._x_ticks,
                         y_range=[-10, 10])

    @property
    def _x_ticks(self):
        labels = self._sanitize_labels(self.plot_data["ticks"]["label"])
        distances = self.plot_data["ticks"]["distance"]
        x_ticks = XTicks(labels=labels, distances=distances)
        return x_ticks

    @property
    def _band_edge(self):
        if self.bs.is_metal():
            return None
        return BandEdge(vbm=self.plot_data["vbm"][0][1],
                        cbm=self.plot_data["cbm"][0][1],
                        vbm_distances=[i[0] for i in self.plot_data["vbm"]],
                        cbm_distances=[i[0] for i in self.plot_data["cbm"]])

    @staticmethod
    def _sanitize_labels(label_list):
        return [italic_to_roman(greek_to_unicode(l)) for l in label_list]
