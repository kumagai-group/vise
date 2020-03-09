# -*- coding: utf-8 -*-

import os
from pymatgen.electronic_structure.core import Spin
from pathlib import Path

from vise.analyzer.dos_plotter import get_dos_plot, max_density
from vise.util.testing import ViseTest

parent_dir = Path(__file__).parent


class ViseDosPlotterNormalTest(ViseTest):

    def setUp(self) -> None:
        self.mgo_dos = get_dos_plot(vasprun=parent_dir / "MgO_dos_vasprun.xml")

    def test(self):
        self.mgo_dos.show()
        self.mgo_dos.savefig(parent_dir / "MgO_dos.pdf")

    def tearDown(self) -> None:
        try:
            os.remove(parent_dir / "MgO_dos.pdf")
        except OSError:
            pass


class ViseDosPlotterWithArgsTest(ViseTest):

    def setUp(self) -> None:
        self.nan3_dos_mod = \
            get_dos_plot(vasprun=parent_dir / "NaN3_dos_vasprun.xml",
                         pdos_type="none",
                         orbital=False,
                         xlim=[-3, 5],
                         ymaxs=[10, 5],
                         zero_at_efermi=False,
                         specific=["1"],
                         crop_first_value=False,
                         show_spg=False)

    def test(self):
        self.nan3_dos_mod.show()
        self.nan3_dos_mod.savefig(parent_dir / "NaN3_dos.pdf")

    def tearDown(self) -> None:
        try:
            os.remove(parent_dir / "NaN3_dos.pdf")
        except OSError:
            pass


class ViseDosPlotterFerromagneticTest(ViseTest):

    def setUp(self) -> None:
        self.mno_dos = \
            get_dos_plot(vasprun=parent_dir / "MnO_uniform_vasprun.xml")

    def test(self):
        self.mno_dos.show()
        self.mno_dos.savefig(parent_dir / "MnO_dos.pdf")

    def tearDown(self) -> None:
        try:
            os.remove(parent_dir / "MnO_dos.pdf")
        except OSError:
            pass


class MaxDensityTest(ViseTest):
    def setUp(self) -> None:
        self.density = {Spin.up: [10.0, 1.0, 2.0, 3.0],
                        Spin.down: [4.0, 8.0, 6.0, 7.0]}
        self.energies = [0.0, 10.0, 20.0, 230.0]
        self.xlim = [5.0, 25.0]
        self.crop_first_value = True

    def test(self):
        expected = 8.0
        actual = max_density(self.density,
                             self.energies,
                             self.xlim,
                             self.crop_first_value)
        self.assertEqual(expected, actual)
