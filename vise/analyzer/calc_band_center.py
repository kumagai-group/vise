# -*- coding: utf-8 -*-

import argparse
from collections import OrderedDict

from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.io.vasp import Vasprun
from vise.analyzer.dos_plotter import ViseDosPlotter

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def calc_o2p_band_center(vasprun_file):
    v = Vasprun(vasprun_file)
    complete_dos = v.complete_dos

    structure = v.final_structure

    total_weight = 0
    summed_dos = 0
    num = 0
    dos = OrderedDict()

    cbm, vbm = complete_dos.get_cbm_vbm()

    for s in structure:
        if s.specie.symbol == "O":
            spd_dos = complete_dos.get_site_spd_dos(s)[OrbitalType.p]
            efermi = spd_dos.efermi
            energies = spd_dos.energies
            densities = spd_dos.densities[Spin.up]
            for e, d in zip(energies, densities):
                if e < efermi:
                    total_weight += d
                    summed_dos += e * d
            num += 1
#            d = complete_dos.get_site_spd_dos(s)[OrbitalType.p]
            dos["Site " + str(num)] = spd_dos

    center = summed_dos / total_weight
    print("(vbm + cbm) / 2: ", (vbm + cbm) / 2)
    print("center: ", center)
    print((vbm + cbm) / 2 - center)

    plotter = ViseDosPlotter(zero_at_efermi=False)
    plotter.add_dos_dict(dos)

    p = plotter.plot_carrier_concentrations(xlim=(-5, 8))
    p.show()
    p.savefig(fname="dos-O2p.eps", format="eps")


def main():
    parser = argparse.ArgumentParser(
        description="""""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", dest="vasprun", type=str)
    args = parser.parse_args()
    calc_o2p_band_center(vasprun_file=args.vasprun)


if __name__ == "__main__":
    main()

