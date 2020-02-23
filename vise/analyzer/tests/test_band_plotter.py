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
                                  absolute=False)
        plotter.show()

