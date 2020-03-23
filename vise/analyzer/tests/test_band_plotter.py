# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import unittest
import warnings

from vise.util.testing import ViseTest
from vise.analyzer.band_plotter import (
    labels_to_unicode, PrettyBSPlotter, ModBSPlotter, make_sym_line)

parent_dir = Path(__file__).parent


class LabelsToUnicodeTest(ViseTest):
    def test(self):
        d = {"GAMMA": {"SIGMA": {"DELTA": "S_0"}}}
        self.assertEqual({'Γ': {'Σ': {'Δ': '{\\rm S}_0'}}},
                         labels_to_unicode(d))


class PlotTwoBandStructureTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        kpoints = parent_dir / "MgO_band_KPOINTS"
        vasprun1 = parent_dir / "MgO_band_vasprun.xml"
        vasprun2 = parent_dir / "MgO_band_nhse_vasprun.xml"
        plotter = PrettyBSPlotter.from_vasp_files(kpoints,
                                                  vasprun1,
                                                  vasprun2,
                                                  absolute=False)
        plotter.show()

    def test_from_vasp_files(self):
        kpoint = parent_dir / "MgO_band_KPOINTS"
        vasprun = parent_dir / "MgO_band_vasprun.xml"
        plotter = PrettyBSPlotter.from_vasp_files(kpoint, vasprun)
        plotter.show()

    def test2(self):
        mgo_band = make_sym_line(
            parent_dir / "MgO_band_KPOINTS",
            parent_dir / "MgO_band_vasprun.xml")
        plotter = PrettyBSPlotter(band=mgo_band, absolute=False)
        plotter.bs_plotter.plot_brillouin()


class PlotFerromagneticBandStructureTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        ko2_band = make_sym_line(
            kpoints_filenames=(parent_dir / "KO2_band_KPOINTS"),
            vasprun_filenames=(parent_dir / "KO2_band_vasprun.xml"))
        plotter = PrettyBSPlotter(band=ko2_band, absolute=False)
        plotter.show()

