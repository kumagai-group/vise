# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import List, Union

import numpy as np
import re
import scipy.interpolate as scint

from matplotlib import pyplot

from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.bandstructure import (
    get_reconstructed_band_structure)
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter, plot_brillouin_zone
from pymatgen.io.vasp import Kpoints, Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import latexify_spacegroup, latexify

from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from vise.util.logger import get_logger
from vise.util.matplotlib import formatter


logger = get_logger(__name__)


class PrettyBSPlotter:

    def __init__(self,
                 band: BandStructureSymmLine,
                 band2: BandStructureSymmLine = None,
                 absolute: bool = False,
                 y_range: List[float] = None,
                 legend: bool = False,
                 symprec: float = SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL,
                 smooth: bool = False) -> None:
        """Band structure plotter class

        Two bands are plotted in the same figure.

        Args:
            band (BandStructureSymmLine):
                A BandStructureSymmLine object.
            band2 (BandStructureSymmLine):
                Another BandStructureSymmLine object to be compared.
            absolute (bool):
                Whether to show the band structure in the absolute energy scale.
            y_range (List[float, float]):
                The min and max of energy range.
            legend (bool):
                Whether to show te figure legend.
            symprec (float):
                Symprec used for determining the space group.
            angle_tolerance (float):
                Angle tolerance used for determining the space group.
        """
        self.bs_plotter = ModBSPlotter(band)

        comp = band.structure.composition.get_reduced_formula_and_factor()[0]
        sga = SpacegroupAnalyzer(structure=band.structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance)

        comp = latexify(comp)
        sg_symbol = latexify_spacegroup(sga.get_space_group_symbol())
        sg_num = sga.get_space_group_number()

        kwargs = {"ylim": y_range,
                  "legend": legend,
                  "zero_to_efermi": not absolute,
                  "title": f"{comp}, SG: {sg_symbol} ({sg_num})",
                  "smooth": smooth}

        if not band2:
            self.plotter = self.bs_plotter.get_plot(**kwargs)
        else:
            bs_plotter2 = ModBSPlotter(band2)
            self.plotter = bs_plotter2.plot_compare(self.bs_plotter, **kwargs)

    @classmethod
    def from_vasp_files(cls,
                        kpoints_filenames: Union[List[str], str],
                        vasprun_filenames: Union[List[str], str],
                        vasprun2_filenames: Union[List[str], str] = None,
                        **kwargs) -> "PrettyBSPlotter":
        """Construct instance from vasp files.

        Compatible with the split multiple vasp calculations.

        Args:
            kpoints_filenames:
                List of KPOINTS file names.
            vasprun_filenames:
                List of vasprun.xml file names for the first band structure.
            vasprun2_filenames:
                List of vasprun.xml file names for the second band structure.
            **kwargs:
                Passed to the constructor. See __int__ docstrings.

        Returns:
            PrettyBSPlotter instance.
        """
        band = make_sym_line(kpoints_filenames, vasprun_filenames)

        if vasprun2_filenames:
            band2 = make_sym_line(kpoints_filenames, vasprun2_filenames)
        else:
            band2 = None

        return cls(band=band, band2=band2, **kwargs)

    def show(self, filename: str = None, format_type: str = "pdf") -> None:
        if filename:
            self.plotter.savefig(filename, format=format_type)
        else:
            self.plotter.show()


def greek_to_unicode(label):
    d = {"GAMMA": "Γ",
         "SIGMA": "Σ",
         "DELTA": "Δ"}
    for k, v in d.items():
        label = label.replace(k, v)
    label = re.sub(r"([A-Z])_", r"{\\rm \1}_", label)  # italic to roman
    return label


def labels_to_unicode(d: dict):
    """Recursively change the Greek label to unicode counterpart."""
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


def make_sym_line(kpoints_filenames: Union[List[str], str],
                  vasprun_filenames: Union[List[str], str]
                  ) -> "BandStructureSymmLine":
    """Create VaspBandStructureSymmLine from vasp files.

    Args:
        kpoints_filenames:
            List of KPOINTS file names.
        vasprun_filenames:
            List of vasprun.xml file names for the first band structure.

    Returns:
        PrettyBSPlotter instance.
    """

    if isinstance(kpoints_filenames, list):
        band = [Vasprun(v).get_band_structure(k, line_mode=True)
                for k, v in zip(kpoints_filenames, vasprun_filenames)]
        bsl = get_reconstructed_band_structure(band)
    else:
        v = Vasprun(vasprun_filenames)
        bsl = v.get_band_structure(kpoints_filenames, line_mode=True)

    st = bsl.structure
    bsl_dict = labels_to_unicode(bsl.as_dict())
    bsl = BandStructureSymmLine.from_dict(bsl_dict)
    bsl.structure = st  # structure is removed during dict round trip.

    return bsl


class ModBSPlotter(BSPlotter):

    def get_plot(self,
                 ylim: List[float] = None,
                 smooth: bool = False,
                 vbm_cbm_marker: bool = True,
                 smooth_tol: float = None,
                 legend: bool = False,
                 zero_to_efermi: bool = True,
                 title: str = None) -> pyplot:
        """Get a matplotlib object for the band structure plot.

        Original function is copied from pymatgen ver.2018.5.22.

        Args:
            ylim (list):
                Specify the y-axis (energy) limits; by default None let the code
                choose. It is vbm-4 and cbm+4 if insulator efermi-10 and
                efermi+10 if metal.
            zero_to_efermi (bool):
                Automatically subtract off the Fermi energy from the eigenvalues
                and plot (E-Ef).
            smooth (bool):
                Whether to interpolates the bands by a spline cubic
            vbm_cbm_marker (bool):
                Whether to show the VBM and CBM markers.
            smooth_tol (float):
                Tolerance for fitting spline to band data. Default is None such
                that no tolerance will be used.
            legend (bool):
                If the figure legend is shown or not.
            title (str):
                The title of the figure.

        Returns:
            Matplotlib pyplot
        """
        plt = pretty_plot(12, 8)

        band_line_width = 1
        spin = [Spin.up, Spin.down] if self._bs.is_spin_polarized else [Spin.up]
        line_style = {Spin.up: 'b-', Spin.down: 'r--'}

        # Ref:
        # http://pymatgen.org/modules/pymatgen/electronic_structure/plotter.html
        data = self.bs_plot_data(zero_to_efermi)
        if not smooth:
            for d in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    for s in spin:
                        data_x = data['distances'][d]
                        data_y = [data['energy'][d][str(s)][i][j]
                                  for j in range(len(data['distances'][d]))]
                        plt.plot(data_x, data_y, line_style[s],
                                 linewidth=band_line_width)
        else:
            # Interpolation failure can be caused by trying to fit an entire
            # band with one spline rather than fitting with piecewise splines
            # (splines are ill-suited to fit discontinuities).
            #
            # The number of splines used to fit a band is determined by the
            # number of branches (high symmetry lines) defined in the
            # BandStructureSymmLine object (see BandStructureSymmLine._branches)

            logger.info("Smoothing wth spline fitting has been performed.")

            warning = "WARNING! Distance / branch {d}, band {i} cannot be " \
                      "interpolated.\n See full warning in source.\n " \
                      "If this is not a mistake, try increasing smooth_tol.\n"\
                      "Current smooth_tol is {s}."

            for d in range(len(data['distances'])):
                # self._nb_bands means the number of bands.
                for i in range(self._nb_bands):
                    for s in spin:
                        # construct splines for each branch, e.g., Gamma-K
                        tck = scint.splrep(
                            data['distances'][d],
                            [data['energy'][d][str(s)][i][j]
                             for j in range(len(data['distances'][d]))],
                            s=smooth_tol)
                        step = (data['distances'][d][-1]
                                - data['distances'][d][0]) / 1000
                        xs = [x * step + data['distances'][d][0]
                              for x in range(1000)]
                        ys = [scint.splev(
                            x * step + data['distances'][d][0], tck, der=0)
                            for x in range(1000)]

                        for y in ys:
                            if np.isnan(y):
                                print(warning.format(d=str(d), i=str(i),
                                                     s=str(smooth_tol)))
                            break

                        plt.plot(xs, ys, line_style[s],
                                 linewidth=band_line_width)

        # change the label size.
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=20)
        # add vertical tics
        self._maketicks(plt)  # use superclass method.

        # Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ vector}$', fontsize=25)
        plt.ylabel(r'$\mathrm{Energy\ (eV)}$', fontsize=25)

        # Draw Fermi energy, only if not zero and metal
        if not zero_to_efermi and self._bs.is_metal():
            ef = self._bs.efermi
            plt.axhline(ef, linewidth=0.5, color='k')

        # X max range = last distance point
        x_max = data['distances'][-1][-1]
        plt.xlim(0, x_max)

        if ylim:
            plt.ylim(ylim)
        else:
            if self._bs.is_metal():
                fl = self._bs.efermi
                e_range = [fl - 10, fl + 10] if zero_to_efermi else [-10, 10]
                plt.ylim(e_range)
            else:
                plt.ylim(data['vbm'][0][1] - 4, data['cbm'][0][1] + 4)

        if not self._bs.is_metal() and vbm_cbm_marker:
            for cbm in data['cbm']:
                plt.scatter(cbm[0], cbm[1], color='r', marker='o', s=100)
            for vbm in data['vbm']:
                plt.scatter(vbm[0], vbm[1], color='g', marker='o', s=100)

        if legend:
            import matplotlib.lines as mlines
            if self._bs.is_spin_polarized:
                handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                         label='bs 1 up'),
                           mlines.Line2D([], [], linewidth=2, color='r',
                                         label='bs 1 down', linestyle="--")]
            else:
                handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                         label='bs 1')]
            plt.legend(handles=handles)

        # plot the band gap with an arrow.
        if self._bs.is_metal() is False:
            band_gap = self._bs.get_band_gap()["energy"]
            x_position = (plt.xlim()[0] + plt.xlim()[1]) / 2
            vbm = self.bs_plot_data(zero_to_efermi=zero_to_efermi)["vbm"][0][1]
            plt = self.add_band_gap(plt, vbm, band_gap, x_position)

        if title:
            plt.title(title, size=23)

        # Modify 3.0 -> 3 in the axes.
        ax.yaxis.set_major_formatter(formatter)

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
                     smooth: bool = False,
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
                If the energy is aligned to the Fermi level or not.
            title (str):
                Title of the plot.

        Returns:
            Matplotlib pyplot.
        """

        plt = self.get_plot(ylim=ylim,
                            zero_to_efermi=zero_to_efermi,
                            smooth=smooth)
        data_orig = self.bs_plot_data(zero_to_efermi=zero_to_efermi)
        data = another_plotter.bs_plot_data(zero_to_efermi=zero_to_efermi)

        spin = [Spin.up, Spin.down] if self._bs.is_spin_polarized else [Spin.up]

        second_line_style = {Spin.up: 'c-', Spin.down: 'm--'}

        # plot second band
        second_band_line_width = 2
        for i in range(another_plotter._nb_bands):
            for d in range(len(data_orig['distances'])):
                for s in spin:
                    data_x = data_orig['distances'][d]
                    data_y = [e[str(s)][i] for e in data['energy']][d]
                    plt.plot(data_x, data_y, second_line_style[s],
                             linewidth=second_band_line_width)

        # add second band gap and its arrow.
        if another_plotter._bs.is_metal() is False:
            band_gap = another_plotter._bs.get_band_gap()["energy"]
            x = (plt.xlim()[0] * 0.35 + plt.xlim()[1] * 0.65)
            vbm = another_plotter.bs_plot_data(
                zero_to_efermi=zero_to_efermi)["vbm"][0][1]
            plt = self.add_band_gap(plt, vbm, band_gap, x)

        import matplotlib.lines as mlines
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
        # get labels and lines
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

    def plot_brillouin(self):
        """Plot the Brillouin zone """
        self.show_brillouin()
