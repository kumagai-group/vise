# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import unittest
import warnings

from vise.util.testing import ViseTest
from vise.analyzer.band_plotter import (
    PrettyBSPlotter, ModBSPlotter, VaspBandStructureSymmLine)

parent_dir = Path(__file__).parent


class PlotTwoBandStructureTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        mgo_band = VaspBandStructureSymmLine(
            parent_dir / "MgO_band_KPOINTS",
            parent_dir / "MgO_band_vasprun.xml", is_projection=True)
        mgo_nhse_band = \
            VaspBandStructureSymmLine(parent_dir / "MgO_band_KPOINTS",
                                      parent_dir / "MgO_band_nhse_vasprun.xml")
        plotter = PrettyBSPlotter(band=mgo_band, band2=mgo_nhse_band,
                                  absolute=False)
        plotter.show()


class PlotFerromagneticBandStructureTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        ko2_band = VaspBandStructureSymmLine(
            kpoints_filename=(parent_dir / "KO2_band_KPOINTS"),
            vasprun_filename=(parent_dir / "KO2_band_vasprun.xml"),
            is_projection=True)
        plotter = PrettyBSPlotter(band=ko2_band, absolute=False)
        plotter.show()

