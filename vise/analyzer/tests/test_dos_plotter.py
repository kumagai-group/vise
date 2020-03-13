# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.electronic_structure.core import Spin
from pathlib import Path
import tempfile
import unittest
import warnings

from vise.analyzer.dos_plotter import get_dos_plot, max_density
from vise.util.testing import ViseTest

parent_dir = Path(__file__).parent


class ViseDosPlotterNormalTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        with tempfile.TemporaryDirectory() as tmp_dirname:
            dos_vasprun = parent_dir / "MgO_dos_vasprun.xml"
            self.mgo_dos = get_dos_plot(vasprun=dos_vasprun)
            self.mgo_dos.show()
            self.mgo_dos.savefig(Path(tmp_dirname) / "MgO_dos.pdf")


class ViseDosPlotterWithArgsTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        nan3_dos_mod = get_dos_plot(vasprun=parent_dir / "NaN3_dos_vasprun.xml",
                                    pdos_type="none",
                                    orbital=False,
                                    xlim=[-3, 5],
                                    ymaxs=[10, 5],
                                    zero_at_efermi=False,
                                    specific=["1"],
                                    crop_first_value=False,
                                    show_spg=False)
        nan3_dos_mod.show()
        with tempfile.TemporaryDirectory() as tmp_dirname:
            nan3_dos_mod.savefig(Path(tmp_dirname) / "NaN3_dos.pdf")


class ViseDosPlotterFerromagneticTest(ViseTest):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    @classmethod
    def tearDownClass(cls):
        warnings.simplefilter("default")

    @unittest.skipIf(not ViseTest.DISPLAY_DIAGRAM, ViseTest.no_display_reason)
    def test(self):
        mno_dos = get_dos_plot(vasprun=parent_dir / "MnO_uniform_vasprun.xml")
        mno_dos.show()
        with tempfile.TemporaryDirectory() as tmp_dirname:
            mno_dos.savefig(Path(tmp_dirname) / "MnO_dos.pdf")


class MaxDensityTest(ViseTest):
    def test(self):
        density = {Spin.up: [10.0, 1.0, 2.0, 3.0],
                   Spin.down: [4.0, 8.0, 6.0, 7.0]}
        energies = [0.0, 10.0, 20.0, 230.0]
        xlim = [5.0, 25.0]
        crop_first_value = True

        expected = 8.0
        actual = max_density(density,
                             energies,
                             xlim,
                             crop_first_value)
        self.assertEqual(expected, actual)
