# -*- coding: utf-8 -*-

import numpy as np
from typing import List, Union
import scipy.interpolate as scint
from matplotlib import pyplot

from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.bandstructure \
    import get_reconstructed_band_structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp import Kpoints, Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import latexify_spacegroup, latexify

from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from vise.util.logger import get_logger


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class PrettyBSPlotter:

    def __init__(self,
                 band: BandStructureSymmLine,
                 band2: BandStructureSymmLine = None,
                 absolute: bool = False,
                 y_range: list = None,
                 legend: bool = False,
                 symprec: float = SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL):

        bs_plotter = ModBSPlotter(band)
        comp = latexify(
            band.structure.composition.get_reduced_formula_and_factor()[0])
        sga = SpacegroupAnalyzer(structure=band.structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance)

        sg_symbol = latexify_spacegroup(sga.get_space_group_symbol())
        sg_num = sga.get_space_group_number()

        kwargs = {"ylim": y_range,
                  "legend": legend,
                  "zero_to_efermi": not absolute,
                  "title": f"{comp}, SG: {sg_symbol} ({sg_num})"}

        if not band2:
            self.plotter = bs_plotter.get_plot(**kwargs)
        else:
            bs_plotter2 = ModBSPlotter(band2)
            self.plotter = bs_plotter2.plot_compare(bs_plotter, **kwargs)

    @classmethod
    def from_vasp_files(cls,
                        kpoints_filenames: Union[List[str], str],
                        vasprun_filenames: Union[List[str], str],
                        vasprun2_filenames: Union[List[str], str] = None,
                        **kwargs) -> "PrettyBSPlotter":

        band = cls.make_sym_line(kpoints_filenames, vasprun_filenames)
        band2 = cls.make_sym_line(kpoints_filenames, vasprun2_filenames) \
            if vasprun2_filenames else None

        return cls(band=band, band2=band2, **kwargs)

    @staticmethod
    def make_sym_line(kpoints_filenames: Union[List[str], str],
                      vasprun_filenames: Union[List[str], str]):
        if isinstance(kpoints_filenames, list):
            band = [VaspBandStructureSymmLine(k, v)
                    for k, v in zip(kpoints_filenames, vasprun_filenames)]
            return get_reconstructed_band_structure(band)
        else:
            return VaspBandStructureSymmLine(kpoints_filenames,
                                             vasprun_filenames)

    def show(self, filename: str = None, format_type: str = "pdf") -> None:
        if filename:
            self.plotter.savefig(filename, format=format_type)
        else:
            self.plotter.show()


class ModBSPlotter(BSPlotter):

    def get_plot(self, ylim=None, smooth=False, vbm_cbm_marker=True,
                 smooth_tol=None, legend=False, zero_to_efermi=True,
                 title=None):
        """Get a matplotlib object for the band structure plot.

        Original function is copied from PYMATGEN.2018.5.22.

        Args:
            zero_to_efermi (bool):
                Automatically subtract off the Fermi energy from the eigenvalues
                and plot (E-Ef).
            ylim (list):
                Specify the y-axis (energy) limits; by default None let the code
                choose. It is vbm-4 and cbm+4 if insulator efermi-10 and
                efermi+10 if metal.
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

            warning = "WARNING! Distance / branch {d}, band {i} cannot be " + \
                      "interpolated.\n" + \
                      "See full warning in source.\n" + \
                      "If this is not a mistake, try increasing " + \
                      "smooth_tol.\nCurrent smooth_tol is {s}."

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

        self._maketicks(plt)

        # Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ vector}$', fontsize=30)
        plt.ylabel(r'$\mathrm{Energy\ (eV)}$', fontsize=30)

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
            # default energy range
            if self._bs.is_metal():
                e_range = [-10, 10]
            else:
                e_range = [-4, 4]

            # metal case
            if self._bs.is_metal():
                if zero_to_efermi:
                    plt.ylim(e_range)
                else:
                    plt.ylim(self._bs.efermi + e_range[0],
                             self._bs.efermi + e_range[1])
            # insulator case
            else:
                plt.ylim(data['vbm'][0][1] + e_range[0],
                         data['cbm'][0][1] + e_range[1])

        if not self._bs.is_metal() and vbm_cbm_marker:
            for cbm in data['cbm']:
                plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                            s=100)
            for vbm in data['vbm']:
                plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                            s=100)

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

        plt.tight_layout()

        # plot the band gap with an arrow.
        if self._bs.is_metal() is False:
            band_gap = self._bs.get_band_gap()["energy"]
            x_position = (plt.xlim()[0] + plt.xlim()[1]) / 2
            vbm = self.bs_plot_data(zero_to_efermi=zero_to_efermi)["vbm"][0][1]
            plt = self.add_band_gap(plt, vbm, band_gap, x_position)

        if title:
            plt.title(title, size=20)

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
        # TODO: When adding the following annotation, the Gamma notation
        # TODO: disappears at x=0.
        plt.annotate(f"{str(round(band_gap, 3))} eV",
                     (annotate_x * 1.03, vbm + band_gap * 0.4), fontsize=20)

        return plt

    def plot_compare(self,
                     another_plotter: "ModBSPlotter",
                     ylim: list = None,
                     legend: bool = True,
                     zero_to_efermi: bool = True,
                     title: str = None) -> pyplot:
        """The x- and y-ranges are determined based on the original data.

        Args:
            another_plotter:
                Another BSPlotter class object.
            ylim:
                Y-axis limit [y min, y max].
            legend:
                If the figure legend is shown or not.
            zero_to_efermi:
                If the energy is aligned to the Fermi level or not.
            title:
                Title of the plot.
        """

        plt = self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi)
        data_orig = self.bs_plot_data(zero_to_efermi=zero_to_efermi)
        data = another_plotter.bs_plot_data(zero_to_efermi=zero_to_efermi)

        if self._bs.is_spin_polarized:
            spin = [Spin.up, Spin.down]
        else:
            spin = [Spin.up]

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

        return plt


class VaspBandStructureSymmLine(BandStructureSymmLine):

    def __init__(self,
                 kpoints_filename: str,
                 vasprun_filename: str,
                 is_projection=False):
        """

        Args:
            kpoints_filename:
                KPOINTS file name parsed.
            vasprun_filename:
                vasprun.xml file name parsed.
            is_projection:
                Whether to calculate the projection values.
        """
        kpoints = Kpoints.from_file(kpoints_filename)

        # BSVasprun cannot be used when the electronic convergence is checked.
        vasprun = Vasprun(filename=vasprun_filename,
                          parse_projected_eigen=is_projection)

        if vasprun.converged_electronic is False:
            logger.warning("SCF is not attained!!!!!")

        eigenvalues = vasprun.eigenvalues
        reciprocal_lattice = vasprun.final_structure.lattice.reciprocal_lattice
        efermi = vasprun.efermi

        # K-points with weight for Hamiltonian are omitted from the plot.
        first_index_wo_weight = sum(w > 0 for w in kpoints.kpts_weights)
        kpts_wo_weight = kpoints.kpts[first_index_wo_weight:]

        eigenvalues_wo_weight = {}
        for s in eigenvalues:
            # transpose is essential
            # When parsing vasprun.xml using BSVasprun,
            # eigenvalues[kpt-index][band-index] = [energy, occupation]
            # For BSPlotter
            # eigenvalues[band-index][kpt-index] = energy
            eigenvalues_wo_weight[Spin(s)] = \
                eigenvalues[Spin(s)][first_index_wo_weight:, :, 0].transpose()

        # Store label except for "None".
        labels_dict = {}
        for i, label in enumerate(kpoints.labels):
            if label != "None" and label is not None:
                if label == "GAMMA":
                    labels_dict["\u0393"] = kpoints.kpts[i]
                elif label == "SIGMA_0":
                    labels_dict["\u03A3" + "_0"] = kpoints.kpts[i]
                elif label == "DELTA_0":
                    labels_dict["\u0394" + "_0"] = kpoints.kpts[i]
                else:
                    labels_dict[label] = kpoints.kpts[i]

        structure = vasprun.final_structure

        # TODO: Write unittest for the projection.
        # TODO: Use projections for the plot.
        projections = None
        if is_projection:
            v = vasprun.projected_eigenvalues
            projections = {}
            for s in eigenvalues:
                # swapaxes is essential
                # When parsing vasprun.xml using BSVasprun,
                #  Vasprun.projected_eigenvalues[spin][kpoint index][band index]
                #                               [atom index][orbital_index]
                # For BSPlotter
                # projections: dict of orbital projections as {spin: ndarray}.
                # The indices of the ndarray are [band_index, kpoint_index,
                #                                 orbital_index, ion_index].
                projections[Spin(s)] = v[Spin(s)].swapaxes(0, 1).swapaxes(2, 3)

        super().__init__(kpoints=kpts_wo_weight,
                         eigenvals=eigenvalues_wo_weight,
                         lattice=reciprocal_lattice,
                         efermi=efermi,
                         labels_dict=labels_dict,
                         structure=structure,
                         projections=projections)
