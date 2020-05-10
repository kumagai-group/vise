# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from collections import OrderedDict, defaultdict
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import palettable
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import Dos
from pymatgen.electronic_structure.dos import add_densities
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.string import latexify_spacegroup, latexify

from vise.defaults import defaults
from vise.util.logger import get_logger

logger = get_logger(__name__)


class ViseDosPlotter(DosPlotter):

    def get_plot(self,
                 xlim: List[float] = None,
                 ylims: List[List[float]] = None,
                 cbm_vbm: List[float] = None,
                 legend: bool = True,
                 crop_first_value: bool = False,
                 title: str = None):
        """Get a matplotlib pyplot for the density of states.

        Override get_plot of DosPlotter.

        Args:
            xlim (list):
                Specifies the x-axis limits. Set to None for automatic
                determination.
            ylims (list):
                Specifies the y-axes limits. Two types of input.
                [[y1min, y1max], [y2min, y2max], ..]
            cbm_vbm (list):
                Specify cbm and vbm [cbm, vbm]
            legend (bool):
                Whether to show the figure legend.
            crop_first_value (bool):
                Whether to crop the fist DOS.
            title (str):
                Title of the figure

        Returns:
            Matplotlib pyplot.
        """
        # From 3 to 9
        ncolors = min(9, max(3, len(self._doses)))
        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        # The DOS calculated using VASP may hold a spuriously large value at the
        # first mesh to keep the consistency with the integrated DOS in DOSCAR
        # file. An example is shown below.
        # ------------- DOSCAR --------------------
        # 10 10 1 0
        # 0.1173120E+02 0.5496895E-09 0.5496895E-09 0.5496895E-09 0.5000000E-15
        # 1.000000000000000E-004
        # CAR
        # unknown system
        #     23.00000000 - 9.00000000 3201 6.62000004 1.00000000
        #     -9.000 0.6000E+03 0.6000E+03 0.6000E+01 0.6000E+01 <-- large DOS
        #     -8.990 0.0000E+00 0.0000E+00 0.6000E+01 0.6000E+01
        i = 1 if crop_first_value else 0

        all_densities = []
        all_energies = []
        for key, dos in self._doses.items():
            energies = dos['energies'][i:]
            densities = {Spin(k): v[i:] for k, v in dos['densities'].items()}
            all_energies.append(energies)
            all_densities.append(densities)

        # Make groups to be shown in the same figure.
        # Example, ZrTiSe4
        # keys = ['Total', 'Site:1 Zr-s', 'Site:1 Zr-p', 'Site:1 Zr-d',
        #         'Site:2 Ti-s', 'Site:2 Ti-p', 'Site:2 Ti-d', 'Site:3 Se-s',
        #         'Site:3 Se-p', 'Site:3 Se-d', 'Site:5 Se-s', 'Site:5 Se-p',
        #         'Site:5 Se-d']
        # grouped_keys =
        #           {"Total"; "Total",
        #            "Site:1": ["Site:1 Zr-s', 'Site:1 Zr-p', 'Site:1 Zr-d'],
        #            ...
        grouped_keys = OrderedDict()
        for k in self._doses.keys():
            first_word = k.split()[0]
            if first_word in grouped_keys:
                grouped_keys[first_word].append(k)
            else:
                grouped_keys[first_word] = [k]

        num_figs = len(grouped_keys)
        print("num_figs", num_figs)
        fig, axs = plt.subplots(num_figs, 1, sharex=True)


        if xlim:
            axs[0].set_xlim(xlim)

        n = 0
        for i, gk in enumerate(grouped_keys):
            for j, key in enumerate(grouped_keys[gk]):
                for spin in list(densities.keys()):
                    x = []
                    y = []
                    if spin in all_densities[n]:
                        density_list = int(spin) * all_densities[n][spin]
                        energies = all_energies[n]
                        x.extend(energies)
                        y.extend(density_list)

                    # Show legend only for spin up.
                    label = str(key) if spin == Spin.up else None
                    color = "black" if key == "Total" else colors[j % ncolors]
                    axs[i].plot(x, y, color=color, label=label, linewidth=2)
                n += 1

            # plot vertical lines for band edges or Fermi level
            if self.zero_at_efermi:
                # plot a line
                axs[i].axvline(0, color="black", linestyle="--", linewidth=0.5)
                if cbm_vbm:
                    axs[i].axvline(cbm_vbm[0] - cbm_vbm[1], color="black",
                                   linestyle="--", linewidth=0.5)
            else:
                axs[i].axvline(self._doses[key]['efermi'],
                               color="black", linestyle="--", linewidth=0.5)
                if cbm_vbm:
                    axs[i].axvline(cbm_vbm[0], color="black", linestyle="--",
                                   linewidth=0.5)

            if legend:
                axs[i].legend(loc="best", markerscale=0.1)
                leg = axs[i].get_legend()
                for legobj in leg.legendHandles:
                    legobj.set_linewidth(1.2)  # legend line width
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize=7)  # legend font size

            axs[i].axhline(0, color="black", linewidth=0.5)

        if ylims and len(ylims) == 2:
            if ylims[0][1] > 1.0:
                axs[0].set_ylim(ylims[0])
            for i in range(1, len(axs)):
                axs[i].set_ylim(ylims[1])
        elif ylims and len(ylims) == num_figs:
            for i in range(len(axs)):
                axs[i].set_ylim(ylims[i])
        elif ylims:
            raise ValueError("The number of y-ranges is not proper.")

        axs[-1].set_xlabel('Energy (eV)')
        axs[0].set_ylabel("Total DOS (1/eV)")
        for i in range(1, len(axs)):
            axs[i].set_ylabel("DOS (1/eV)")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                            wspace=0, hspace=0.1)

        if title:
            axs[0].title.set_text(title)

        return plt


def get_dos_plot(vasprun: str,
                 cbm_vbm: list = None,
                 pdos_type: str = "element",
                 specific: list = None,
                 orbital: bool = True,
                 xlim: list = None,
                 ymaxs: list = None,
                 zero_at_efermi: bool = True,
                 legend: bool = True,
                 crop_first_value: bool = True,
                 show_spg: bool = True,
                 symprec: float = defaults.symmetry_length_tolerance,
                 angle_tolerance: float = defaults.symmetry_angle_tolerancde):
    """Get

    Args:
        vasprun (str):
            vasprun.xml-type file name
        cbm_vbm (list):
            List of [cbm, vbm]. This is used when band edge is determined from
            the band structure calculation.
        pdos_type (str):
            Plot type of PDOS.
            "element": PDOS grouped by element type
            "site": PDOS grouped by equivalent sites
            "none": PDOS are not grouped, so all dos are shown individually.
        specific (list):
            Show specific PDOS. If list elements are integers,
            PDOS at particular sites are shown. If elements are shown,
            PDOS of particular elements are shown.
            ["1", "2"] --> At site 1 and 2 compatible with pdos_type = "none"
            ["Mg", "O"] --> Summed at Mg and O sites compatible with
                            pdos_type = "element"
        orbital (bool):
            Whether to show orbital decomposed PDOS.
        xlim (list):
            Specifies x-axis limits.
            Set to None for automatic determination.
        ymaxs (list):
            Specifies the maxima of absolute y-axis limits.
        zero_at_efermi (bool):
            Whether to show the plot in the absolute scale.
        legend (bool):
            Whether to show the figure legend.
        crop_first_value (bool):
            Whether to crop the fist DOS.
        show_spg (bool):
            Whether to show space group number in the title.
        symprec (float):
            Symprec for determining the equivalent sites.
        angle_tolerance (float):
            Angle tolerance for symmetry analysis in degree.
    """
    v = Vasprun(vasprun, ionic_step_skip=True, parse_eigen=False)
    if v.converged_electronic is False:
        logger.warning("Be aware SCF is not attained in the vasp calculation.")

    structure = v.final_structure
    complete_dos = v.complete_dos
    # check cbm
    if cbm_vbm is None:
        if v.complete_dos.get_gap() > 0.1:
            cbm_vbm = complete_dos.get_cbm_vbm()

    dos = OrderedDict()
    dos["Total"] = complete_dos  # CompleteDos has same format as Dos.

    if specific and specific[0].isdigit():
        if pdos_type is not "none":
            logger.warning(f"pdos_type is changed from {pdos_type} to none.")
        pdos_type = "none"

    elif specific and specific[0].isalpha():
        if pdos_type is not "none":
            logger.warning(f"pdos_type is changed from {pdos_type} to element.")
        pdos_type = "element"

    grouped_indices = defaultdict(list)
    sga = None
    if pdos_type == "element":
        for indices, s in enumerate(structure):
            grouped_indices[str(s.specie)].append(indices)

    elif pdos_type == "site":
        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        sga = SpacegroupAnalyzer(structure=structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance)
        symmetrized_structure = sga.get_symmetrized_structure()
        # equiv_indices = [[0], [1], [2, 3], [4, 5]]
        equiv_index_lists = symmetrized_structure.equivalent_indices

        for l in equiv_index_lists:
            specie = structure[l[0]].specie
            wyckoff = sga.get_symmetry_dataset()["wyckoffs"][l[0]]

            name = f"{specie} {wyckoff}"
            if name in grouped_indices:
                for i in range(2, 10):
                    name = f"{specie} {wyckoff}{i}"
                    if name in grouped_indices:
                        continue
                    else:
                        break
                else:
                    raise ValueError(
                        "The situation where more than 9 same Wyckoff sites "
                        "exist is not supported.")
            grouped_indices[name] = l
            logger.info({f"Atoms {grouped_indices[name]} grouped to {name}."})

    elif pdos_type == "none":
        for indices, s in enumerate(structure):
            name = f"{s.specie} site: {indices}"
            grouped_indices[name].append(indices)
    else:
        raise KeyError("The given pdos_type is not supported.")

    if specific:
        tmp = defaultdict(list)
        for key, value in grouped_indices.items():
            if pdos_type == "element" and key in specific:
                tmp[key] = value
            else:
                # type(index) is str
                index = ''.join(c for c in key if c.isdigit())
                if index in specific:
                    tmp[key] = value
        grouped_indices = tmp

    # efermi is set to VBM if exists.
    efermi = cbm_vbm[1] if cbm_vbm else complete_dos.efermi
    complete_dos.efermi = efermi
    energies = complete_dos.relative_energies

    for key, value in grouped_indices.items():
        for indices in value:
            site = structure[indices]
            if orbital:
                for orb, pdos in complete_dos.get_site_spd_dos(site).items():
                    # " " is used for grouping the plots.
                    if pdos_type == "none":
                        name = f"{key} {orb}"
                    else:
                        name = f"{key} #{len(value)} {orb}"
                    density = div_densities(pdos.densities, len(value))
                    if name in dos:
                        density = add_densities(dos[name].densities, density)
                    dos[name] = Dos(efermi, energies, density)
            else:
                name = f"{key}({len(key)})" if isinstance(key, list) else key
                pdos = complete_dos.get_site_dos(site)
                if name in dos:
                    density = add_densities(dos[name].densities, pdos.densities)
                    dos[name] = Dos(efermi, energies, density)
                else:
                    dos[name] = pdos

    # use complete_dos.efermi for total dos.
    plotter = ViseDosPlotter(zero_at_efermi=zero_at_efermi)
    plotter.add_dos_dict(dos)

    xlim = xlim or [-10, 10]

    if ymaxs:
        if v.incar.get("ISPIN", 1) == 2:
            ylims = [[-y, y] for y in ymaxs]
        else:
            ylims = [[0, y] for y in ymaxs]
    else:
        energies = complete_dos.relative_energies - efermi
        tdos_max = 1.1 * max_density(density=complete_dos.densities,
                                     energies=energies,
                                     xlim=xlim,
                                     crop_first_value=crop_first_value)
        if v.incar.get("ISPIN", 1) == 2:
            ylims = [[-tdos_max, tdos_max]]
        else:
            ylims = [[0, tdos_max]]

        pdos_max = 0.0
        for k, d in dos.items():
            if k == "Total":
                continue
            pdos_max = max(max_density(d.densities, energies, xlim), pdos_max)
        pdos_max *= 1.1

        if v.incar.get("ISPIN", 1) == 2:
            ylims.append([-pdos_max, pdos_max])
        else:
            ylims.append([0, pdos_max])

    x_str = f"{xlim[0]:7.2f} ~{xlim[1]:6.2f}"
    y_str = f"total: {ylims[0][0]:7.2f} ~{ylims[0][1]:6.2f}"
    if len(ylims) == 2:
        y_str += f", pdos: {ylims[1][0]:7.2f} ~{ylims[1][1]:6.2f}"

    logger.info(f"Dos plot, x-range (energies) {x_str}")
    logger.info(f"Dos plot, y-range (dos) {y_str}")

    comp = latexify(structure.composition.get_reduced_formula_and_factor()[0])
    if show_spg:
        if sga is None:
            sga = SpacegroupAnalyzer(structure, symprec=symprec)
        sg_num_str = str(sga.get_space_group_number())
        sg_symbol = latexify_spacegroup(sga.get_space_group_symbol())
        sg = f" {sg_symbol} ({sg_num_str})"
        logger.info(f"Space group: {sg}")
        title = f"{comp}, SG: {sg}"
    else:
        title = comp

    return plotter.get_plot(xlim=xlim,
                            ylims=ylims,
                            cbm_vbm=cbm_vbm,
                            legend=legend,
                            crop_first_value=crop_first_value,
                            title=title)


def div_densities(density: dict, denominator: float) -> Dict[Spin, np.ndarray]:
    """Method to divide density. """
    return {spin: np.array(val) / denominator for spin, val in density.items()}


def max_density(density: dict,
                energies: list,
                xlim: list,
                crop_first_value: bool = True) -> float:
    """Evaluate max value of the density of states in the given energies range.

    Args:
        density (dict):
            Note that the first value may contain a huge value when the
            lower limit of the calculation of density of states is larger than
            that of occupied states. Therefore, we need to crop the first value
            by default.
            {Spin.up: [...], Spin.down: [...] }
        energies (list):
            Energy mesh values.
        xlim (list):
             x-range with [x-min, x-max]
        crop_first_value (bool):
            Whether to crop the first value or not.

    Returns:
        Max value in the density within the given x-range.
    """
    values = []
    for density_in_each_spin in density.values():
        for i, (d, e) in enumerate(zip(density_in_each_spin, energies)):
            if crop_first_value and i == 0:
                continue
            if xlim[0] < e < xlim[1]:
                values.append(d)

    if not values:
        raise ValueError("DOS is empty at the given energies {0[0]} - {0[1]} "
                         "range.".format(xlim))

    return max(values)
