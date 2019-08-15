# -*- coding: utf-8 -*-

import unittest

from pymatgen.electronic_structure.bandstructure \
    import get_reconstructed_band_structure

from vise.analyzer.band_plotter import ModBSPlotter, VaspBandStructureSymmLine

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class PlotBandStructureTest(unittest.TestCase):

    def test(self):
        band0 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/vasprun-finish.xml")
        band1 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/vasprun-finish.xml")
        band2 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/vasprun-finish.xml")
        band_a = get_reconstructed_band_structure([band0, band1, band2])
#        print(dir(band_a))
#        print(band_a.branches)
#        print(band_a.bands)
#        print(band_a.distance)
#        band_gap_a = band_a.get_band_gap()["distances"]
        bs_plotter_a = ModBSPlotter(band_a)
#        data = bs_plotter_a.bs_plot_data()
#        print(data["distances"])
#        p_a = bs_plotter_a.get_plot(zero_to_efermi=False, smooth=False)
#        p_a.show()
        band3 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/0/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/0/vasprun-finish.xml")
        band4 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/1/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/1/vasprun-finish.xml")
        band5 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/2/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/2/vasprun-finish.xml")
        band_b = get_reconstructed_band_structure([band3, band4, band5])
        # scissored_band_b = band_a.apply_scissor(band_gap_a)
        # bs_plotter_b = ModBSPlotter(scissored_band_b)
        bs_plotter_b = ModBSPlotter(band_b)
        p_b = bs_plotter_b.plot_compare(bs_plotter_a, zero_to_efermi=False)
#        p_b = bs_plotter_b.plot_compare(bs_plotter_a)
        p_b.show()




        # band0 = get_band_structure(kpoints_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/KPOINTS",
        #                            vasprun_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/vasprun-finish.xml",
        #                            is_projection=True,
        #                            poscar_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/POSCAR")
        # band1 = get_band_structure(kpoints_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/KPOINTS",
        #                            vasprun_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/vasprun-finish.xml",
        #                            is_projection=True,
        #                            poscar_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/POSCAR")
        # band2 = get_band_structure(kpoints_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/KPOINTS",
        #                            vasprun_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/vasprun-finish.xml",
        #                            is_projection=True,
        #                            poscar_name="/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/POSCAR")

        # band = get_reconstructed_band_structure([band0, band1, band2])
#        d = band0.get_projection_on_elements()
#        print(d)
#        plot_projected_band_structure(band)


if __name__ == '__main__':
    unittest.main()
