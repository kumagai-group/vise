# -*- coding: utf-8 -*-


from pymatgen.electronic_structure.core import Spin

from vise.analyzer.dos_plotter import ViseDosPlotter, get_dos_plot, \
    divide_densities, max_density
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ViseDosPlotterTest(ViseTest):

    def setUp(self) -> None:
        self.mgo_dos = get_dos_plot(vasprun_file="MgO_dos_vasprun.xml")
        self.nan3_dos_mod = get_dos_plot(vasprun_file="NaN3_dos_vasprun.xml",
                                        pdos_type="none",
                                        orbital=False,
                                        xlim=[-3, 5],
                                        ymaxs=[10, 5],
                                        zero_at_efermi=False,
                                        specific=["1"],
                                        crop_first_value=False,
                                        show_spg=False,
                                        )

    def test_normal(self):
        self.mgo_dos.show()

    def test_w_args(self):
        self.nan3_dos_mod.show()


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
