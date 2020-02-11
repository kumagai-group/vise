#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.

from pathlib import Path

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter

from vise.chempotdiag.compound import (
    Compound, CompoundsList, DummyCompoundForDiagram)
from vise.chempotdiag.chem_pot_diag import ChemPotDiag
from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE, MOLECULE_SUFFIX
from vise.util.testing import ViseTest


#MP_TEST = True

class TestChemPotDiag(ViseTest):

    def setUp(self) -> None:

        mg = PDEntry(Composition("Mg"), -1.0)
        ca = PDEntry(Composition("Ca"), -2.0)
#        o = PDEntry(Composition("O"), -3.0)
        camg2 = PDEntry(Composition("Ca2Mg4"), -25.0)
#        camgo2 = PDEntry(Composition("Ca2MgO4"), -60.0)

#        self.phase_diagram = PhaseDiagram([mg, ca, o, camg2, camgo2])
        self.phase_diagram = PhaseDiagram([mg, ca, camg2])

    def test_cpd(self):
        pdp = PDPlotter(self.phase_diagram)
#        PDPlotter(self.phase_diagram).plot_chempot_range_map([Element("Mg"), Element("Ca")])
#        pdp.show()
        print(self.phase_diagram)
        print(self.phase_diagram.facets)
        cpd = ChemPotDiag(pd=self.phase_diagram)
        print(cpd.vertices)
        print(cpd.absolute_chem_pot(0))
        cpd.draw_diagram(remarked_comp=Composition("Ca2Mg4"))
