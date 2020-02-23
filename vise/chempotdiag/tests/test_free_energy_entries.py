# -*- coding: utf-8 -*-

import os
from pathlib import Path
import shutil

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element

from vise.chempotdiag.free_energy_entries import (
    FreeEnergyEntry, FreeEnergyEntrySet, ConstrainedFreeEnergyEntrySet)
from vise.util.testing import ViseTest


DISABLE_DISPLAY_DIAGRAM = False
parent_dir = Path(__file__).parent


class TestFreeEnergyEntry(ViseTest):
    def setUp(self) -> None:
        self.mgo = FreeEnergyEntry("MgO", total_energy=-10)
        self.o2 = FreeEnergyEntry("O2",
                                  total_energy=-2,
                                  zero_point_vib=-10,
                                  free_e_shift=-100,
                                  data={"id": 10},
                                  attribute="new")

    def test_dict(self):
        expected = self.o2.as_dict()
        actual = FreeEnergyEntry.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)

    def test_msonable(self):
        self.assertMSONable(self.mgo)
        self.assertMSONable(self.o2)

    def test_print(self):
        expected = "PDEntry : O2 with composition O2. Energy: -112.000. Zero point vib energy: -10.000. Free energy contribution: -100.000. Data: {'id': 10}"
        self.assertEqual(expected, str(self.o2))


class TestFreeEnergyEntrySet(ViseTest):
    def setUp(self) -> None:
        mgo = FreeEnergyEntry("MgO", total_energy=-10)
        self.mg = FreeEnergyEntry("Mg", total_energy=-10)
        self.o2 = FreeEnergyEntry("O2",
                                  total_energy=-2,
                                  zero_point_vib=-10,
                                  free_e_shift=-100,
                                  data={"id": 10},
                                  attribute="new")

        self.entry_set = FreeEnergyEntrySet([mgo, self.o2])

    def test_msonable(self):
        self.assertMSONable(self.entry_set)

    def test_add_discard(self):
        es = self.entry_set
        self.assertEqual(2, len(es))
        es.add(self.mg)
        self.assertEqual(3, len(es))
        es.discard(self.o2)
        self.assertEqual(2, len(es))

    def test_csv(self):
        self.entry_set.to_csv("entries.csv")
        actual = open("entries.csv").read()
        expected = """Name,Mg,O,Energy
MgO,1.0,1.0,-10
O2,0,2.0,-112
"""
        self.assertEqual(expected, actual)
        os.remove("entries.csv")

    def test_from_vasp_files(self):
        paths = [parent_dir / "vasp_Mg",
                 parent_dir / "mol_O2",
                 parent_dir / "vasp_MgO"]
        entry_set = FreeEnergyEntrySet.from_vasp_files(paths,
                                                       parse_gas=True,
                                                       temperature=200)
        expected = {"Mg", "MgO", "O2"}
        actual = {str(e.name) for e in entry_set}

        self.assertEqual(expected, actual)
        expected = {-5.18827492, -12.51242284, -9.20808429}
        actual = {e.energy for e in entry_set}
        self.assertEqual(expected, actual)
        print(entry_set)

    def test_from_mp(self):
        entry_set = FreeEnergyEntrySet.from_mp(["Mg", "O"])
        expected = {"Mg", "MgO", "O2"}
        actual = {str(e.name) for e in entry_set}
        self.assertEqual(expected, actual)

        expected = {-39.42964153, -14.43901464, -11.96804153}
        actual = {e.energy for e in entry_set}
        self.assertEqual(expected, actual)
        print(entry_set)

    def test_dict(self):
        expected = self.o2.as_dict()
        actual = FreeEnergyEntry.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)

    def test_msonable(self):
        self.assertMSONable(self.entry_set)

    def test_print(self):
        print(self.entry_set)


class TestConstrainedFreeEnergyEntrySet(ViseTest):
    def setUp(self) -> None:

        mg = FreeEnergyEntry("Mg", total_energy=-1)
        ca = FreeEnergyEntry("Ca", total_energy=-1)
        mgo = FreeEnergyEntry("MgO", total_energy=-12)
        cao = FreeEnergyEntry("CaO", total_energy=-12)
        o2 = FreeEnergyEntry("O2", total_energy=-2)
        mgcao = FreeEnergyEntry("MgCaO", total_energy=-33)
        self.entry_set = FreeEnergyEntrySet([mg, ca, mgo, cao, o2, mgcao])

    def test_constraint(self):
        entry_set = ConstrainedFreeEnergyEntrySet.\
            from_entry_set(self.entry_set, {"Ca": -2})
        expected = {Element.Mg, Element.O}
        actual = set()
        for e in entry_set.entries:
            actual.update(e.composition.elements)
        self.assertEqual(expected, actual)

        expected = FreeEnergyEntry(composition=Composition({"Mg": 1, "O": 1}),
                                   total_energy=-33-(-2),
                                   name="CaMgO")
        actual = [e for e in entry_set if e.name == "CaMgO"][0]
        self.assertEqual(str(expected), str(actual))


