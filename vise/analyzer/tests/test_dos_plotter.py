# -*- coding: utf-8 -*-

import os
from pymatgen.electronic_structure.core import Spin
from pathlib import Path

from vise.analyzer.dos_plotter import get_dos_plot, max_density
from vise.util.testing import ViseTest

parent_dir = Path(__file__).parent


class ViseDosPlotterTest(ViseTest):

    # Must split pyplot constructors. Otherwise, mgo_dos.pdf shows NaN3 results.
    def test_normal(self):
        mgo_dos = get_dos_plot(vasprun=parent_dir / "MgO_dos_vasprun.xml")
        mgo_dos.savefig("mgo_dos.pdf")
        os.remove("mgo_dos.pdf")

    def test_w_args(self):
        nan3_dos_mod = get_dos_plot(vasprun=parent_dir / "NaN3_dos_vasprun.xml",
                                    pdos_type="none",
                                    orbital=False,
                                    xlim=[-3, 5],
                                    ymaxs=[10, 5],
                                    zero_at_efermi=False,
                                    specific=["1"],
                                    crop_first_value=False,
                                    show_spg=False,
                                    )
        nan3_dos_mod.savefig("nan3_dos.pdf")
        os.remove("nan3_dos.pdf")


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
