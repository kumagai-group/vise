# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import re
from itertools import product
from typing import List

from matplotlib import lines as mlines
from matplotlib import pyplot
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter, plot_brillouin_zone
from pymatgen.io.vasp import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import latexify_spacegroup, latexify

from vise.defaults import defaults
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter

logger = get_logger(__name__)


def greek_to_unicode(label: str) -> str:
    """Return the greek letter for the special points in reciprocal space."""
    d = {"GAMMA": "Γ",
         "SIGMA": "Σ",
         "DELTA": "Δ"}
    for k, v in d.items():
        label = label.replace(k, v)
    label = re.sub(r"([A-Z])_", r"{\\rm \1}_", label)  # italic to roman
    return label


def labels_to_unicode(d: dict) -> dict:
    """Recursively change the Greek label to unicode counterpart in dict."""
    new_d = {}
    for k, v in d.items():
        if isinstance(k, str):
            k = greek_to_unicode(k)
        if isinstance(v, str):
            v = greek_to_unicode(v)
        elif isinstance(v, dict):
            v = labels_to_unicode(v)
        new_d[k] = v

    return new_d


def make_bs_sym_lines(kpoints: str,
                      vaspruns: List[str]) -> List["BandStructureSymmLine"]:
    """Create VaspBandStructureSymmLine from vasp files.

    Args:
        kpoints (str):
            KPOINTS file names.
        vaspruns (List[str]):
            Vasprun.xml file names for the first band structure.

    Returns:
        BandStructureSymmLine instance.
    """
    def bs_sym_line(vr: str):
        bsl = Vasprun(vr).get_band_structure(kpoints, line_mode=True)
        st = bsl.structure
        bsl_dict = labels_to_unicode(bsl.as_dict())
        bsl = BandStructureSymmLine.from_dict(bsl_dict)
        bsl.structure = st  # structure is removed during dict round trip.
        return bsl

    return [bs_sym_line(v) for v in vaspruns]


class PrettyBSPlotter:

    def __init__(self,
                 bands: List[BandStructureSymmLine],
                 title: str = None,
                 absolute: bool = False,
                 y_range: List[float] = None,
                 legend: bool = False,
                 symprec: float = defaults.symmetry_length_tolerance,
                 angle_tolerance: float = defaults.symmetry_angle_tolerance
                 ) -> None:
        """Band structure plotter class

        Two bands are plotted in the same figure.

        Args:
            bands (List[BandStructureSymmLine]):
                List of BandStructureSymmLine object.
            absolute (bool):
                Whether to show the band structure in the absolute energies scale.
            y_range (List[float, float]):
                The min and max of energies range.
            legend (bool):
                Whether to show te figure legend.
            symprec (float):
                Symprec used for determining the space group.
            angle_tolerance (float):
                Angle tolerance used for determining the space group.
        """
        self.bs = bands
        self.absolute = absolute
        self.y_range = y_range
        self.legend = legend

        sga = SpacegroupAnalyzer(structure=bands[0].structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance)

        composition = bands[0].structure.composition
        comp_latex = latexify(composition.get_reduced_formula_and_factor()[0])
        sg_symbol = latexify_spacegroup(sga.get_space_group_symbol())
        sg_num = sga.get_space_group_number()

        self.title = title or f"{comp_latex}, SG: {sg_symbol} ({sg_num})"
        self.plotter = None

    def show(self, filename: str = None, format_type: str = "pdf") -> None:
        if self.plotter is None:
            raise ValueError

        if filename:
            self.plotter.savefig(filename, format=format_type)
        else:
            self.plotter.show()

    def get_plot(self, vbm_cbm_marker: bool = True) -> pyplot:
        """Get a matplotlib object for the band structure plot.

        Args:
            ylim (list):
                Specify the y-axis (energies) limits; by default None let the code
                choose. It is vbm-4 and cbm+4 if insulator efermi-10 and
                efermi+10 if set_metal_condition.
            vbm_cbm_marker (bool):
                Whether to show the VBM and CBM markers.
        Returns:
            Matplotlib pyplot
        """
        plt = pretty_plot(12, 8)
        bs_plotter = BSPlotter(self.bs[0])

        plt = bs_plotter._maketicks(plt)  # add vertical tics

        for b in self.bs:
            bs_data = bs_plotter.bs_plot_data(zero_to_efermi=not self.absolute)
            spin = [Spin.up, Spin.down] if b.is_spin_polarized else [Spin.up]
            num_distance = len(bs_data["distances"])
            num_band = len(bs_data["energies"])

            # bs_data["energies"][distance][spin][band]
            for d, e, s in product(range(num_distance), range(num_band), spin):
                data_x = bs_data["distances"][num_distance]
                data_y = [e[str(s)][num_band] for e in bs_data['energies'][d]]
                line_style = "b-" if s == Spin.up else "r--"
                plt.plot(data_x, data_y, line_style, linewidth=1)

            # Draw Fermi level, only if not zero and set_metal_condition
            if self._bs.is_metal():
                efermi = 0.0 if zero_to_efermi else self._bs.efermi
                plt.axhline(efermi, linewidth=0.75, color='k')
            else:
                band_gap = self._bs.get_band_gap()["energies"]
                x_position = (plt.xlim()[0] + plt.xlim()[1]) / 2
                # plot the band gap with an arrow.
                vbm = self.bs_plot_data(zero_to_efermi=zero_to_efermi)["vbm"][0][1]
                plt = self.add_band_gap(plt, vbm, band_gap, x_position)

            plt.xlim(0, bs_data['distances'][-1][-1])  # max range = last dist point

            if not self._bs.is_metal() and vbm_cbm_marker:
                for cbm in bs_data['cbm']:
                    plt.scatter(cbm[0], cbm[1], color='r', marker='o', s=100)
                for vbm in bs_data['vbm']:
                    plt.scatter(vbm[0], vbm[1], color='g', marker='o', s=100)

            if self.legend:
                if self._bs.is_spin_polarized:
                    handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                             label='up'),
                               mlines.Line2D([], [], linewidth=2, color='r',
                                             label='down', linestyle="--")]
                else:
                    handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                             label='bs')]
                plt.legend(handles=handles)

        # if self.y_range:
        #     plt.ylim(self.y_range)
        # else:
        #     if self.bs1.is_metal():
        #         fl = self.bs.efermi
        #         e_range = [fl - 10, fl + 10] if zero_to_efermi else [-10, 10]
        #         plt.ylim(e_range)
        #     else:
        #         plt.ylim(bs_data['vbm'][0][1] - 4, bs_data['cbm'][0][1] + 4)

        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.xlabel('Wave vector', fontsize=25)
        plt.ylabel('Energy (eV)', fontsize=25)
        if self.title:
            plt.title(self.title, size=24)

        # Modify 3.0 -> 3 in the axes.
        ax.yaxis.set_major_formatter(float_to_int_formatter)

        plt.tight_layout()
        return plt

    @staticmethod
    def add_band_gap(plt: pyplot,
                     vbm: float,
                     band_gap: float,
                     annotate_x: float,
                     line_color: str = "purple") -> pyplot:
        plt.hlines([vbm, vbm + band_gap], plt.xlim()[0], plt.xlim()[1],
                   line_color, linestyles='dashed')

        plt.annotate('', xy=(annotate_x, vbm),
                     xytext=(annotate_x, vbm + band_gap),
                     fontsize=11, color='black',
                     arrowprops=dict(edgecolor='black', arrowstyle='<|-|>',
                                     shrinkA=0, shrinkB=0))
        plt.annotate(f"{str(round(band_gap, 3))} eV",
                     (annotate_x * 1.03, vbm + band_gap * 0.4), fontsize=20)

        return plt

    def plot_compare(self,
                     another_plotter: "ModBSPlotter",
                     ylim: List[float] = None,
                     legend: bool = True,
                     zero_to_efermi: bool = True,
                     title: str = None,
                     ) -> pyplot:
        """Generate pyplot for comparing two band structures.

        Note: x- and y-ranges are determined based on the original data.

        Args:
            another_plotter (ModBSPlotter):
                Another BSPlotter class object.
            ylim (List[float]):
                Y-axis limit [y min, y max].
            legend (bool):
                If the figure legend is shown or not.
            zero_to_efermi (bool):
                If the energies is aligned to the Fermi level or not.
            title (str):
                Title of the plot.

        Returns:
            Matplotlib pyplot.
        """
        plt = self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi)
        data_orig = self.bs_plot_data(zero_to_efermi=zero_to_efermi)
        another_bs_data = \
            another_plotter.bs_plot_data(zero_to_efermi=zero_to_efermi)

        spin = [Spin.up, Spin.down] if self._bs.is_spin_polarized else [Spin.up]

        second_line_style = {Spin.up: 'c-', Spin.down: 'm--'}

        # plot second band
        second_band_line_width = 2
        for i in range(another_plotter._nb_bands):
            for d in range(len(data_orig['distances'])):
                for s in spin:
                    data_x = data_orig['distances'][d]
                    data_y = [e[str(s)][i] for e in another_bs_data['energies']][d]
                    plt.plot(data_x, data_y, second_line_style[s],
                             linewidth=second_band_line_width)

        # add second band gap and its arrow.
        if another_plotter._bs.is_metal() is False:
            band_gap = another_plotter._bs.get_band_gap()["energies"]
            x = (plt.xlim()[0] * 0.35 + plt.xlim()[1] * 0.65)
            vbm = another_plotter.bs_plot_data(
                zero_to_efermi=zero_to_efermi)["vbm"][0][1]
            plt = self.add_band_gap(plt, vbm, band_gap, x)

        if legend:
            if self._bs.is_spin_polarized:
                handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                         label='bs 1 up'),
                           mlines.Line2D([], [], linewidth=2, color='r',
                                         label='bs 1 down', linestyle="--")]
            else:
                handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                         label='bs 1')]
            if another_plotter._bs.is_spin_polarized:
                handles.extend([mlines.Line2D([], [], linewidth=2, color='c',
                                              label='bs 2 up'),
                                mlines.Line2D([], [], linewidth=2, color='m',
                                              linestyle="--",
                                              label='bs 2 down')])
            else:
                handles.extend([mlines.Line2D([], [], linewidth=2, color='c',
                                              label='bs 2')])
            plt.legend(handles=handles)

        if title:
            plt.title(title, size=20)

        plt.tight_layout()
        return plt

    def show_brillouin(self, filename: str = None):
        """Show the Brillouin zone """
        labels = {}
        for k in self._bs.kpoints:
            if k.label:
                labels[k.label + str(k.frac_coords)] = k.frac_coords

        plot_brillouin_zone(self._bs.lattice_rec,
                            labels=labels,
                            tight_layout=True,
                            title="Brillouin zone",
                            show=False if filename else True,
                            savefig=filename)

