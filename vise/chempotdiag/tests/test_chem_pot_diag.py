#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.

import unittest
from pathlib import Path
import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter

from vise.chempotdiag.chem_pot_diag import ChemPotDiag, sort_coords
from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE, MOLECULE_SUFFIX
from vise.util.testing import ViseTest


DISABLE_DISPLAY_DIAGRAM = False


class TestChemPotDiag(ViseTest):

    def setUp(self) -> None:

        mg = PDEntry(Composition("Mg"), -1.0)
        ca = PDEntry(Composition("Ca"), -2.0)
        sr = PDEntry(Composition("Sr"), -1.0)
        o = PDEntry(Composition("O"), -3.0)
        camg2 = PDEntry(Composition("Ca2Mg4"), -25.0)
        ca2mgo4 = PDEntry(Composition("Ca2MgO4"), -60.0)
        ca2mgsro4 = PDEntry(Composition("Ca2MgSrO4"), -160.0)

        phase_diagram = PhaseDiagram([mg, ca, o, camg2, ca2mgo4])
#        phase_diagram = PhaseDiagram([mg, ca, sr, o, camg2, ca2mgo4, ca2mgsro4])
#        phase_diagram = PhaseDiagram([mg, ca, camg2])
        self.cpd = ChemPotDiag.from_phase_diagram(
            pd=phase_diagram, target_composition="CaMg2")

    def test_msonable(self):
        self.assertMSONable(self.cpd)

    @unittest.skipIf(DISABLE_DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_cpd(self):
        print(self.cpd.target_comp_chempot)

        #        pdp = PDPlotter(phase_diagram)
#        PDPlotter(self.phase_diagram).plot_chempot_range_map([Element("Mg"), Element("Ca")])
#        pdp.show()
#        print(self.phase_diagram)
#        print(self.phase_diagram.facets)
#        print(self.cpd.target_comp_chempot)
        print(self.cpd.target_comp_abs_chempot)
#        print(cpd.vertices)
#        print(cpd.absolute_chem_pot(0))
        self.cpd.draw_diagram()


class TestSortCoords(ViseTest):
    def setUp(self) -> None:
        # x + 2y + 3z = 4
        self.coords = \
            np.array([[3, 2, -1], [-1, 2, 0], [-6, -1, 4], [1, -3, 3]])

    def test_sort_coords(self):
#        print(self.coords)
        print(sort_coords(self.coords))

