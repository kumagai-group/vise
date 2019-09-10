# -*- coding: utf-8 -*-

from pymatgen.util.testing import PymatgenTest

from pymatgen.electronic_structure.core import Spin

from vise.analyzer.dos_plotter import ModDosPlotter, get_dos_plot, \
    divide_densities, max_density

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class PlotDosTest(PymatgenTest):

    def test(self):
#        dos = get_dos_plot(vasprun_file="vasprun.xml", orbital=True)
        dos = get_dos_plot(vasprun_file="vasprun.xml",
#                           sites=[1, 2, 3, 55, 56],
                           orbital=True,
                           xlim=[-20, 5],
#                           ymaxs=[500, 5, 5, 5, 1, 1],
                           ymaxs=[500, 5],
                           zero_at_efermi=False,
                           legend=True)
#                           zero_at_efermi=True)
        dos.show()

#         band0 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/0/vasprun-finish.xml")
#         band1 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/1/vasprun-finish.xml")
#         band2 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2/2/vasprun-finish.xml")
#         band_a = get_reconstructed_band_structure([band0, band1, band2])
# #        print(dir(band_a))
#         print(band_a.branches)
# #        print(band_a.bands)
# #        print(band_a.distance)
# #        band_gap_a = band_a.get_band_gap()["distances"]
#         bs_plotter_a = ModBSPlotter(band_a)
#         data = bs_plotter_a.bs_plot_data()
#         print(data["distances"])
#         p_a = bs_plotter_a.get_plot(zero_to_efermi=False, smooth=True)
#         p_a.show()
#         band3 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/0/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/0/vasprun-finish.xml")
#         band4 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/1/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/1/vasprun-finish.xml")
#         band5 = VaspBandStructureSymmLine("/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/2/KPOINTS", "/Users/kuma/my_programs/pydefect_test_MgSe/unitcell/band2_scan/2/vasprun-finish.xml")
#         band_b = get_reconstructed_band_structure([band3, band4, band5])
#         # scissored_band_b = band_a.apply_scissor(band_gap_a)
#         # bs_plotter_b = ModBSPlotter(scissored_band_b)
#         bs_plotter_b = ModBSPlotter(band_b)
#         p_b = bs_plotter_b.plot_compare(bs_plotter_a, zero_to_efermi=False)
# #        p_b = bs_plotter_b.plot_compare(bs_plotter_a)
#         p_b.show()




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


class MaxDensityTest(PymatgenTest):
    def setUp(self) -> None:
        self.density = \
            {Spin.up: [0.0, 1.0, 2.0, 3.0], Spin.down: [4.0, 8.0, 6.0, 7.0]}
        self.energies = [0.0, 10.0, 20.0, 230.0]
        self.xlim = [5.0, 25.0]
        self.crop_first_value = True

    def test(self):
        actual = max_density(self.density, self.energies, self.xlim, self.crop_first_value)
        print(actual)
