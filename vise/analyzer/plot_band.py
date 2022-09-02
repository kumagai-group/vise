# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
from dataclasses import dataclass
from itertools import cycle
from typing import List, Optional, Tuple, Union, Dict

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MSONable
from num2words import num2words

from vise.error import ViseError
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.logger import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class XTicks(MSONable):
    """Info used for X-axis ticks in the band structure plot.

    labels: Special point names, e.g., "Γ", "C$\mid$${\rm C}_2$",
    distances: Distances to the corresponding labels. So, the length must be the
    same as that of labels.
    """
    labels: List[str]
    distances: List[float]


@dataclass
class BandEdge(MSONable):
    vbm: float
    cbm: float
    # The positions of the VBM and CBM denoted by total distances in the band
    # figure.
    vbm_distances: List[float]
    cbm_distances: List[float]


class BandInfo(MSONable):
    def __init__(self,
                 # [branch][spin][band][k-point]
                 # A branch is an area separated by a vertical bar.
                 band_energies: List[List[List[List[float]]]],
                 band_edge: Optional[BandEdge] = None,
                 fermi_level: Optional[float] = None,
                 # [branch][spin](k-point idx, [irrep names])
                 irreps: List[List[Tuple[int, List[str]]]] = None):
        self.band_energies = deepcopy(band_energies)
        self.band_edge = deepcopy(band_edge)
        self.fermi_level = fermi_level
        self.irreps = irreps

        if self.band_edge is None and self.fermi_level is None:
            raise ViseBandInfoError

    def slide_energies(self, base_energy):
        self._slide_band_energies(base_energy)
        self._slide_band_edge(base_energy)
        self._slide_fermi_level(base_energy)

    def _slide_band_energies(self, base_energy):
        new_array = []
        for band_energies_each_branch in self.band_energies:
            inner_array = []
            for band_energies_each_band in band_energies_each_branch:
                a = np.array(band_energies_each_band)
                inner_array.append((a - base_energy).tolist())
            new_array.append(inner_array)
        self.band_energies = new_array

    def _slide_fermi_level(self, base_energy):
        if self.fermi_level:
            self.fermi_level -= base_energy

    def _slide_band_edge(self, base_energy):
        if self.band_edge:
            self.band_edge.vbm -= base_energy
            self.band_edge.cbm -= base_energy

    @property
    def is_magnetic(self):
        return len(self.band_energies[0]) == 2

    def band_energy_region(self,
                           decision_width: float = 0.1,
                           bottom: float = None,
                           top: float = None,
                           offset: float = None) -> List[List[float]]:
        """Estimate continuously filled band energy region.

        Args:
            decision_width: Max range over which energy is considered continuous
            bottom: minimum energy to be considered
            top: maximum energy to be considered
            offset: subtract offset value from all the band energy

        Returns:
            energy regions: [[1st min, 1st max], [2nd min, 2nd max], ...]
        """
        result = []

        def add_boundary(lower, upper):
            if offset is not None:
                result.append([lower - offset, upper - offset])
            else:
                result.append([lower, upper])

        sorted_energies = sorted([energy
                                  for i in self.band_energies
                                  for j in i
                                  for k in j
                                  for energy in k])  # flatten nested list.
        if bottom is not None:
            sorted_energies = list(filter(lambda x: x >= bottom,
                                          sorted_energies))
        if top is not None:
            sorted_energies = list(filter(lambda x: x <= top,
                                          sorted_energies))

        prev_energy = sorted_energies.pop(0)
        lower_bound = prev_energy
        for energy in sorted_energies:
            if energy - prev_energy > decision_width:
                upper_bound = prev_energy
                add_boundary(lower_bound, upper_bound)
                lower_bound = energy  # update lower_bound
            prev_energy = energy
        else:
            upper_bound = energy  # last energy should be the upper bound.
            add_boundary(lower_bound, upper_bound)

        return result


class ViseBandInfoError(ViseError):
    pass


class BandMplSettings:
    def __init__(self,
                 colors: Optional[List[str]] = None,
                 linewidth: List[float] = None,
                 circle_size: int = 70,
                 circle_colors: Optional[List[str]] = None,
                 band_edge_line_width: float = 0.75,
                 band_edge_line_color: str = "black",
                 band_edge_line_style: str = "-.",
                 title_font_size: int = 15,
                 label_font_size: int = 15,
                 show_legend: bool = True,
                 legend_location: str = "lower right"
                 ):
        self.colors = colors or ['#E15759', '#4E79A7', '#F28E2B', '#76B7B2']
        self.linewidth = cycle(linewidth) if linewidth else cycle([1.0])

        self.circle_size = circle_size
        self.hline = {"linewidth": band_edge_line_width,
                      "color": band_edge_line_color,
                      "linestyle": band_edge_line_style}

        self.title_font_size = title_font_size
        self.label_font_size = label_font_size

        self.circle_colors = circle_colors or ["pink", "green"]
        self.circle_size = circle_size

        self.show_legend = show_legend
        self.legend = {"loc": legend_location}

    def band_structure(self, index, label):
        result = {"color": self.colors[index],
                  "linewidth": next(self.linewidth)}
        if index > 0:
            result["label"] = label
        return result

    def circle(self, index):
        return {"color": self.colors[index],
                "marker": "o",
                "s": self.circle_size}


@dataclass
class BandPlotInfo(MSONable, ToJsonFileMixIn):
    """Multiple BandInfo are accepted.

    Ex: include both the PBE band and GW band
    """
    band_infos: Dict[str, BandInfo]  # keys are subtitles.
    distances_by_branch: List[List[float]]
    x_ticks: XTicks
    title: str = None  # title of all plots

    @classmethod
    def from_dict(cls, d):
        """Needed for backward compatibility"""
        if "band_info_set" in d:
            d["band_infos"] = d.pop("band_info_set")

        if isinstance(d["band_infos"], list):
            d["band_infos"] = {str(i): x for i, x in enumerate(d["band_infos"])}
        return super().from_dict(d)

    def __post_init__(self):
        assert self.distances_by_branch[0][0] == self.x_ticks.distances[0]
        assert self.distances_by_branch[-1][-1] == self.x_ticks.distances[-1]

    def __add__(self, other: "BandPlotInfo"):
        assert self.distances_by_branch == other.distances_by_branch
        new_band_info_set = deepcopy(self.band_infos)
        new_band_info_set.update(other.band_infos)
        return BandPlotInfo(new_band_info_set, self.distances_by_branch,
                            self.x_ticks, self.title)


class BandPlotter:

    def __init__(self,
                 band_plot_info: BandPlotInfo,
                 energy_range: List[float],
                 base_energy: float = None,
                 base_energy_title: str = None,
                 mpl_defaults: Optional[BandMplSettings] = BandMplSettings()
                 ):
        # need deepcopy to avoid side effect caused by the energy shift
        band_plot_info = deepcopy(band_plot_info)
        self.band_infos = band_plot_info.band_infos
        self.distances_by_branch = band_plot_info.distances_by_branch
        self.x_ticks = band_plot_info.x_ticks

        self.energy_range = energy_range
        self.title = band_plot_info.title
        self.mpl_defaults = mpl_defaults
        self.plt = plt

        if base_energy is None:
            if base_energy_title is None:
                base_energy_title = next(iter(self.band_infos))
                logger.warning(f"Base energy is set to {base_energy_title}.")

            base_band_info = self.band_infos[base_energy_title]

            if base_band_info.band_edge is not None:
                base_energy = base_band_info.band_edge.vbm
            elif self.band_infos[base_energy_title].fermi_level:
                base_energy = base_band_info.fermi_level
            else:
                logger.warning(f"Band edge and Fermi level are absent in "
                               f"{base_energy_title}, so absolute energy"
                               f"is used for plot.")
                base_energy = 0.0

        self._slide_energies(base_energy)

    def _slide_energies(self, base_energy):
        for band_info in self.band_infos.values():
            band_info.slide_energies(base_energy)

    def construct_plot(self):
        self._add_band_set()
        self._set_figure_legend()
        self._set_x_range()
        self._set_y_range()
        self._set_labels()
        self._set_x_ticks()
        self._set_title()
        self._set_formatter()
        self.plt.tight_layout()

    def _add_band_set(self):
        for index, (band_name, band_info) in enumerate(self.band_infos.items()):
            self._add_band_structures(band_info, band_name, index)

            if band_info.band_edge:
                self._add_band_edge(band_info.band_edge, index)
            else:
                if band_info.fermi_level:
                    self._add_fermi_level(band_info.fermi_level)

    def _add_band_structures(self, band_info, band_name, index):
        mpl_args = self.mpl_defaults.band_structure(index, band_name)
        for distances, energies_each_branch \
                in zip(self.distances_by_branch,  band_info.band_energies):
            for spin_index, energies_by_spin in enumerate(energies_each_branch):
                if spin_index == 0:
                    mpl_args.pop("linestyle", None)
                else:
                    mpl_args["linestyle"] = ":"

                for energies_of_a_band in energies_by_spin:
                    self.plt.plot(distances, energies_of_a_band, **mpl_args)
                    mpl_args.pop("label", None)

    def _add_band_edge(self, band_edge, index):
        self.plt.axhline(y=band_edge.vbm, **self.mpl_defaults.hline)
        self.plt.axhline(y=band_edge.cbm, **self.mpl_defaults.hline)

        for dist in band_edge.vbm_distances:
            self.plt.scatter(dist, band_edge.vbm,
                             **self.mpl_defaults.circle(index))
        for dist in band_edge.cbm_distances:
            self.plt.scatter(dist, band_edge.cbm,
                             **self.mpl_defaults.circle(index))

    def _add_fermi_level(self, fermi_level):
        self.plt.axhline(y=fermi_level, **self.mpl_defaults.hline)

    def _set_figure_legend(self):
        if self.mpl_defaults.show_legend and len(self.band_infos) > 1:
            self.plt.legend(**self.mpl_defaults.legend)

    def _set_x_range(self):
        self.plt.xlim(self.x_ticks.distances[0], self.x_ticks.distances[-1])

    def _set_y_range(self):
        self.plt.ylim(self.energy_range[0], self.energy_range[1])

    def _set_labels(self):
        self.plt.xlabel("Wave vector", size=self.mpl_defaults.label_font_size)
        self.plt.ylabel("Energy (eV)", size=self.mpl_defaults.label_font_size)

    def _set_x_ticks(self):
        axis = self.plt.gca()
        axis.set_xticks(self.x_ticks.distances)
        axis.set_xticklabels(self.x_ticks.labels)
        for distance, label in zip(self.x_ticks.distances[1:-1],
                                   self.x_ticks.labels[1:-1]):
            linestyle = "-" if "\\mid" in label else "--"
            plt.axvline(x=distance, linestyle=linestyle)

    def _set_title(self):
        if self.title:
            self.plt.title(self.title, size=self.mpl_defaults.title_font_size)

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.yaxis.set_major_formatter(float_to_int_formatter)


