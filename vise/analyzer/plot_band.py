# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
from dataclasses import dataclass
from itertools import cycle
from typing import List, Optional, Dict

import matplotlib.pyplot as plt
from monty.json import MSONable

from vise.error import ViseError
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.logger import get_logger

logger = get_logger(__name__)


@dataclass
class Irrep(MSONable):
    frac_coords: List[float]
    symbols: List[str]
    energies: List[float]
    degeneracies: List[int]

    @property
    def irrep_info_set(self):
        return zip(self.symbols, self.energies, self.degeneracies)


@dataclass
class Irreps(MSONable):
    sg_num: int
    # key is special point name. Gamma = GM
    irreps: Dict[str, Irrep]

    def __call__(self):
        return self.irreps

    def get_distances(self, x_ticks: "XTicks") -> List[List[float]]:
        result = []
        for special_point_symbol in self.irreps:
            dist_list = []
            for label, distance in zip(x_ticks.labels, x_ticks.distances):
                if special_point_symbol in label:
                    dist_list.append(distance)
            result.append(dist_list)
        return result


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
class BandEdgeForPlot(MSONable):
    vbm: float
    cbm: float
    # The positions of the VBM and CBM denoted by total distances in the band
    # figure.
    vbm_distances: List[float]
    cbm_distances: List[float]


# backward compatibility
@dataclass
class BandEdge(BandEdgeForPlot):
    pass


@dataclass
class BandEnergyInfo(MSONable):
    # [branch][spin][band][k-point] A branch is an area in which
    # the k-points are continuous. Each branch is separated by a vertical bar.
    # We need to distinguish branch to draw continuous line in the area and
    # calculate the effective masses.
    band_energies: List[List[List[List[float]]]]
    band_edge: Optional[BandEdgeForPlot] = None
    fermi_level: Optional[float] = None
    irreps: Optional[Irreps] = None

    def __post_init__(self):
        if self.band_edge is None and self.fermi_level is None:
            raise ViseBandInfoError

    def slide_energies(self, base_energy):
        self._slide_band_energies(base_energy)
        self._slide_band_edge(base_energy)
        self._slide_fermi_level(base_energy)
        self._slide_irreps(base_energy)

    def _slide_band_energies(self, base_energy):
        self.band_energies = [[[[w - base_energy for w in x] for x in y]
                               for y in z] for z in self.band_energies]

    def _slide_band_edge(self, base_energy):
        if self.band_edge:
            self.band_edge.vbm -= base_energy
            self.band_edge.cbm -= base_energy

    def _slide_fermi_level(self, base_energy):
        if self.fermi_level:
            self.fermi_level -= base_energy

    def _slide_irreps(self, base_energy):
        if self.irreps:
            for i in self.irreps.irreps.values():
                i.energies = [e - base_energy for e in i.energies]

    @property
    def is_magnetic(self):
        return len(self.band_energies[0]) == 2

    def band_energy_region(self,
                           decision_width: float = 0.1,
                           bottom: float = None,
                           top: float = None,
                           offset: float = 0.0) -> List[List[float]]:
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
            result.append([lower - offset, upper - offset])

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


@dataclass
class BandPlotInfo(MSONable, ToJsonFileMixIn):
    """Multiple BandInfo with the same k-point path are accepted.

    Ex: include both the PBE band and GW band
    """
    band_energy_infos: Dict[str, BandEnergyInfo]  # keys are subtitles.
    distances_by_branch: List[List[float]]
    x_ticks: XTicks
    title: str = None  # title of all plots

    def __post_init__(self):
        assert self.distances_by_branch[0][0] == self.x_ticks.distances[0]
        assert self.distances_by_branch[-1][-1] == self.x_ticks.distances[-1]

    def __add__(self, other: "BandPlotInfo"):
        assert self.distances_by_branch == other.distances_by_branch
        new_band_info_set = deepcopy(self.band_energy_infos)
        new_band_info_set.update(other.band_energy_infos)
        return BandPlotInfo(new_band_info_set, self.distances_by_branch,
                            self.x_ticks, self.title)


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
                  "linewidth": next(self.linewidth),
                  "label": label}
        return result

    def circle(self, index):
        return {"color": self.colors[index],
                "marker": "o",
                "s": self.circle_size}


class BandMplPlotter:

    def __init__(self,
                 band_plot_info: BandPlotInfo,
                 energy_range: List[float],
                 base_energy: float = None,
                 base_energy_title: str = None,
                 mpl_defaults: Optional[BandMplSettings] = BandMplSettings()
                 ):
        # need deepcopy to avoid side effect caused by the energy shift
        band_plot_info = deepcopy(band_plot_info)
        self.band_energy_infos = band_plot_info.band_energy_infos
        self.distances_by_branch = band_plot_info.distances_by_branch
        self.x_ticks = band_plot_info.x_ticks

        self.energy_range = energy_range
        self.title = band_plot_info.title
        self.mpl_defaults = mpl_defaults
        self.plt = plt

        if base_energy is None:
            base_energy = get_base_energy(self.band_energy_infos, base_energy_title)
        slide_band_energies(self.band_energy_infos, base_energy)

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
        for index, (band_name, band_info) in enumerate(self.band_energy_infos.items()):
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
        if self.mpl_defaults.show_legend and len(self.band_energy_infos) > 1:
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


def get_base_energy(band_energy_infos, base_energy_title=None):
    if base_energy_title is None:
        base_energy_title = next(iter(band_energy_infos))
        if len(band_energy_infos) > 1:
            logger.warning(f"Base energy is set to {base_energy_title}.")

    base_band_info = band_energy_infos[base_energy_title]

    if base_band_info.band_edge:
        return base_band_info.band_edge.vbm

    return base_band_info.fermi_level


def slide_band_energies(band_energy_infos, base_energy):
    for band_info in band_energy_infos.values():
        band_info.slide_energies(base_energy)

