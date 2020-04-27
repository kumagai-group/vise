# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from itertools import tee, groupby
import re
from typing import List, Dict, Optional
import matplotlib.pyplot as plt
from dataclasses import dataclass

from pymatgen import Spin


@dataclass()
class XTics:
    # labels include "x|y"
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
    band_energies: Dict[Spin, List[List[float]]]
    band_edge: Optional[BandEdge] = None
    fermi_level: Optional[float] = None


class BandPlotter:
    def __init__(self,
                 band_info_set: List[BandInfo],
                 distances: List[float],
                 x_ticks: XTics,
                 y_range: List[float],
                 ):
        self._band_info_set = band_info_set
        self._distances = distances
        self._x_ticks = x_ticks
        self._y_range = y_range

        assert self._distances[0] == self._x_ticks.distances[0]
        assert self._distances[-1] == self._x_ticks.distances[-1]

        self.plt = plt
        self._set_band_set()
        self._set_x_range()
        self._set_y_range()
        self._set_axes()
        self._set_xticks()
#        ax.yaxis.set_major_formatter(formatter)
#        self._set_major_fo()
        self.plt.tight_layout()

    def _set_band_set(self):
        for band_info in self._band_info_set:
            self._set_band_structures(band_info)
            if band_info.band_edge:
                self._set_band_edge(band_info.band_edge)
            if band_info.fermi_level:
                self._set_fermi_level(band_info.fermi_level)

    def _set_band_structures(self, band_info):
        for spin, band_energies in band_info.band_energies.items():
            for energies_per_band in band_energies:
                self.plt.plot(self._distances, energies_per_band)

    def _set_band_edge(self, band_edge):
        self.plt.axhline(y=band_edge.vbm)
        for dist in band_edge.vbm_distances:
            self.plt.scatter(dist, band_edge.vbm)

        self.plt.axhline(y=band_edge.cbm)
        for dist in band_edge.cbm_distances:
            self.plt.scatter(dist, band_edge.cbm)

    def _set_fermi_level(self, fermi_level):
        self.plt.axhline(y=fermi_level)

    def _set_x_range(self):
        self.plt.xlim(self._distances[0], self._distances[-1])

    def _set_y_range(self):
        self.plt.ylim(self._y_range[0], self._y_range[1])

    def _set_axes(self):
        self.plt.xlabel("Wave vector")
        self.plt.ylabel("Energy (eV)")

    def _set_xticks(self):
        axis = self.plt.gca()
        axis.set_xticks(self._x_ticks.distances)
        axis.set_xticklabels(self._x_ticks.labels)
        for distance, label in zip(self._x_ticks.distances[1:-1],
                                   self._x_ticks.labels[1:-1]):
            if "|" in label:
                plt.axvline(x=distance)
            else:
                plt.axvline(x=distance, linestyle="--")


def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def greek_to_unicode(label: str) -> str:
    d = {"GAMMA": "Γ", "SIGMA": "Σ", "DELTA": "Δ"}
    for k, v in d.items():
        label = label.replace(k, v)
    return label


def italic_to_roman(label: str) -> str:
    return re.sub(r"([A-Z])_([0-9])", r"${\\rm \1}_\2$", label)



# class Xtics:
#     def __init__(self, labels: List[str], distances: List[float]):
#         assert all(prior <= posterior
#                    for prior, posterior in pairwise(distances))
#
#         self.distances = []
#         self.labels = []
#         self.vertical_line_positions = []
#         distances_labels = list(zip(distances, labels))
#         for distance, group in groupby(distances_labels, key=lambda x: x[0]):
#             self.distances.append(distance)
#             labels = [g[1] for g in group]
#             label = "|".join(labels)
#             self.labels.append(label)
#             if len(labels) > 1:
#                 self.vertical_line_positions.append(distance)
#
#     @staticmethod
#     def sanitize_label(label):
#         return italic_to_roman(greek_to_unicode(label))


