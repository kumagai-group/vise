# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.analyzer.band_plotter import (
    PrettyBSPlotter, ModBSPlotter, VaspBandStructureSymmLine)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class PlotBandStructureTest(ViseTest):
    def setUp(self) -> None:
        self.mgo_band = VaspBandStructureSymmLine("MgO_band_KPOINTS",
                                                  "MgO_band_vasprun.xml",
                                                  is_projection=True)
        self.mgo_nhse_band = \
            VaspBandStructureSymmLine("MgO_band_KPOINTS",
                                      "MgO_band_nhse_vasprun.xml")

    def test(self):
        plotter = PrettyBSPlotter(band=self.mgo_band, band2=self.mgo_nhse_band,
                                  absolute=True)
        plotter.show()


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

