# -*- coding: utf-8 -*-

import unittest
import numpy as np
from unittest.mock import patch

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element
from pymatgen.analysis.phase_diagram import (
    PDEntry, PhaseDiagram, CompoundPhaseDiagram, PDPlotter)

from vise.chempotdiag.chem_pot_diag import ChemPotDiag, sort_coords
from vise.util.testing import ViseTest


DISABLE_DISPLAY_DIAGRAM = False


class TestChemPotDiag(ViseTest):

    def setUp(self) -> None:

        mg = PDEntry(Composition("Mg"), -1.0)
        mg2 = PDEntry(Composition("Mg"), -0.5)
        ca = PDEntry(Composition("Ca"), -2.0)
        sr = PDEntry(Composition("Sr"), -3.0)
        o = PDEntry(Composition("O"), -4.0)

        camg = PDEntry(Composition("Ca2Mg2"), -16.0)  # rel -10.0
        camgo = PDEntry(Composition("CaMgO"), -17.0)  # rel -10.0
        camgo2 = PDEntry(Composition("CaMgO2"), -11.0)  # rel 0.0
        camgsro4 = PDEntry(Composition("CaMgSrO4"), -32.0)  # rel -10

        self.pd_1d = PhaseDiagram([mg, mg2])
        self.pd_2d = PhaseDiagram([mg, ca, camg])
        self.pd_3d = PhaseDiagram([mg, ca, camg, o, camgo, camgo2])
        self.pd_4d = PhaseDiagram([mg, ca, camg, o, camgo, sr, camgsro4])

        self.cpd_1d = ChemPotDiag.from_phase_diagram(pd=self.pd_1d,
                                                     target_comp="Mg")
        self.cpd_2d = ChemPotDiag.from_phase_diagram(pd=self.pd_2d,
                                                     target_comp="CaMg")
        self.cpd_3d = ChemPotDiag.from_phase_diagram(pd=self.pd_3d,
                                                     target_comp="CaMgO")
        self.cpd_3d_unstable = \
            ChemPotDiag.from_phase_diagram(pd=self.pd_3d,
                                           target_comp="CaMgO2",
                                           allow_unstable_target_chempot=True)
        self.cpd_4d = ChemPotDiag.from_phase_diagram(pd=self.pd_4d,
                                                     target_comp="SrCaMgO4")
#        target_comp="CaMgSrO4")

        self.comp_pd = CompoundPhaseDiagram(
            entries=[mg, ca, o, camg, camgo, camgo2],
            terminal_compositions=[Composition("CaMg"),
                                   Composition("Ca"),
                                   Composition("O")])

        self.comp_cpd = ChemPotDiag.from_phase_diagram(
            pd=self.comp_pd,
            target_comp="MgCaO2",
            allow_unstable_target_chempot=True)

    def test_cpd_1d(self):
        self.assertEqual([Element.Mg], self.cpd_1d.elements)
        self.assertEqual(1, self.cpd_1d.dim)
        self.assertEqual([[0.0]], self.cpd_1d.vertices)
        self.assertEqual("Mg", self.cpd_1d.target_comp)
        self.assertEqual({'A': [0.0]}, self.cpd_1d.target_comp_chempot)
        self.assertEqual({'A': [-1.0]}, self.cpd_1d.target_comp_abs_chempot)

    # def test_pd_1d_draw(self):
    #     pdp = PDPlotter(self.pd_1d, show_unstable=10)
    #     pdp.show()

    def test_cpd_1d_draw(self):
        with self.assertRaises(NotImplementedError):
            self.cpd_1d.draw_diagram()

    def test_cpd_2d(self):
        self.assertEqual([Element.Ca, Element.Mg], self.cpd_2d.elements)
        self.assertEqual(2, self.cpd_2d.dim)
        self.assertEqual([[0.0, -5.0], [-5.0, 0.0]], self.cpd_2d.vertices)
#        self.assertEqual([[-5.0, 0.0], [0.0, -5.0]], self.cpd_2d.vertices)
        self.assertEqual("CaMg", self.cpd_2d.target_comp)
        self.assertEqual({'A': [0.0, -5.0], 'B': [-5.0, 0.0]},
                         self.cpd_2d.target_comp_chempot)
        self.assertEqual({'A': [-2.0, -6.0], 'B': [-7.0, -1.0]},
                         self.cpd_2d.target_comp_abs_chempot)

    @unittest.skipIf(DISABLE_DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_cpd_2d_draw(self):
        self.cpd_2d.draw_diagram()

    def test_cpd_2d_draw_show(self):
        path = "vise.chempotdiag.chem_pot_diag.plt.show"
        with patch(path) as show_patch:
            self.cpd_2d.draw_diagram()
            assert show_patch.called

    def test_cpd_2d_draw_savefig(self):
        path = "vise.chempotdiag.chem_pot_diag.plt.savefig"
        with patch(path) as show_patch:
            self.cpd_2d.draw_diagram(filename="a.pdf")
            show_patch.assert_called_once_with("a.pdf")

    def test_cpd_3d(self):
        self.assertEqual([Element.Ca, Element.Mg, Element.O],
                         self.cpd_3d.elements)
        self.assertEqual(3, self.cpd_3d.dim)
        self.assertEqual([[-10.0, 0.0, 0.0], [-5.0, 0.0, -5.0],
                          [0.0, -10.0, 0.0], [0.0, -5.0, -5.0]],
                         self.cpd_3d.vertices)
        self.assertEqual("CaMgO", self.cpd_3d.target_comp)
        self.assertEqual({'A': [-10.0, 0.0, 0.0], 'B': [-5.0, 0.0, -5.0],
                          'C': [0.0, -10.0, 0.0], 'D': [0.0, -5.0, -5.0]},
                         self.cpd_3d.target_comp_chempot)

    @unittest.skipIf(DISABLE_DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_cpd_3d_draw(self):
        self.cpd_3d.draw_diagram()

    def test_cpd_3d_unstable(self):
        self.assertEqual({'A': [-10.0, 0.0, 0.0], 'B': [0.0, -10.0, 0.0]},
                         self.cpd_3d_unstable.target_comp_chempot)

    @unittest.skipIf(DISABLE_DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_cpd_3d_unstable_draw(self):
        self.cpd_3d_unstable.draw_diagram()

    def test_cpd_4d(self):
        self.assertEqual([Element.Sr, Element.Ca, Element.Mg, Element.O],
                         self.cpd_4d.elements)
        self.assertEqual([[0.0, -10.0, 0.0, 0.0], [0.0, 0.0, -10.0, 0.0],
                          [0.0, -5.0, 0.0, -5.0], [0.0, 0.0, -5.0, -5.0]],
                         self.cpd_4d.vertices)

    def test_pd_4d_plot(self):
        pdp = PDPlotter(self.pd_4d)
        pdp.show()

    def test_cpd_4d(self):
        # print(self.comp_cpd.elements)
        # print(self.comp_cpd.el_ref_list)
        # print(self.comp_cpd.dim)
        # print(self.comp_cpd.vertices)
        # print(self.comp_cpd.qhull_entries)
        # print(self.comp_cpd.comp_facets)
        # print(self.comp_cpd.target_comp)
        pdp = PDPlotter(self.comp_pd)
        pdp.show()

    def test_comp_cpd_4d(self):
#        self.comp_pd.draw_diagram()
        print(self.comp_cpd.target_comp_abs_chempot)
        self.comp_cpd.draw_diagram()


class TestSortCoords(ViseTest):
    def setUp(self) -> None:
        # x + 2y + 3z = 4
        self.coords = \
            np.array([[3, 2, -1], [-1, 2, 0], [-6, -1, 4], [1, -3, 3]])

    def test_sort_coords(self):
#        print(self.coords)
        print(sort_coords(self.coords))

