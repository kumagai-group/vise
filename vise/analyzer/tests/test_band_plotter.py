# -*- coding: utf-8 -*-

from pathlib import Path

from vise.util.testing import ViseTest
from vise.analyzer.band_plotter import (
    PrettyBSPlotter, ModBSPlotter, VaspBandStructureSymmLine)

parent_dir = Path(__file__).parent


class PlotTwoBandStructureTest(ViseTest):
    def setUp(self) -> None:
        self.mgo_band = VaspBandStructureSymmLine(
            parent_dir / "MgO_band_KPOINTS",
            parent_dir / "MgO_band_vasprun.xml", is_projection=True)
        self.mgo_nhse_band = \
            VaspBandStructureSymmLine(parent_dir / "MgO_band_KPOINTS",
                                      parent_dir / "MgO_band_nhse_vasprun.xml")

    def test(self):
        plotter = PrettyBSPlotter(band=self.mgo_band, band2=self.mgo_nhse_band,
                                  absolute=False)
        plotter.show()


class PlotFerromagneticBandStructureTest(ViseTest):
    def setUp(self) -> None:
        self.ko2_band = VaspBandStructureSymmLine(
            kpoints_filename=(parent_dir / "KO2_band_KPOINTS"),
            vasprun_filename=(parent_dir / "KO2_band_vasprun.xml"),
            is_projection=True)

    def test(self):
        plotter = PrettyBSPlotter(band=self.ko2_band, absolute=False)
        plotter.show()

