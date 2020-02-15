# -*- coding: utf-8 -*-

import unittest
import numpy as np
from copy import deepcopy
from pathlib import Path

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element

from vise.chempotdiag.free_energy_entries import (
    FreeEnergyEntry, FreeEnergyEntrySet, ConstrainedFreeEnergyEntrySet)
from vise.util.testing import ViseTest


DISABLE_DISPLAY_DIAGRAM = False


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
        expected = "PDEntry : O2 with composition O2 and energy = -112.0000"
        self.assertEqual(expected, str(self.o2))


class TestFreeEnergyEntrySet(ViseTest):
    def setUp(self) -> None:
        mgo = FreeEnergyEntry("MgO", total_energy=-10)
        self.o2 = FreeEnergyEntry("O2",
                                  total_energy=-2,
                                  zero_point_vib=-10,
                                  free_e_shift=-100,
                                  data={"id": 10},
                                  attribute="new")

        self.mg = FreeEnergyEntry("Mg", total_energy=-10)
        self.entry_set = FreeEnergyEntrySet([mgo, self.o2])

    def test_msonable(self):
        self.assertMSONable(self.entry_set)

    def test_add_discard(self):
        es = self.entry_set
        print(len(es))
        es.add(self.mg)
        print(len(es))
        es.discard(self.o2)
        print([str(i) for i in es])
        print(len(es))

    def test_csv(self):
        self.entry_set.to_csv("entries.csv")
        a = FreeEnergyEntrySet.from_csv("entries.csv")
        print([str(i) for i in a])

    def test_from_vasp_files(self):
        paths = [Path(".") / "vasp_Mg",
                 Path(".") / "vasp_O2",
                 Path(".") / "vasp_MgO"]
        entry_set = FreeEnergyEntrySet.from_vasp_files(paths, parse_gas=True)
        print([e for e in entry_set])

    def test_from_mp(self):
        entry_set = FreeEnergyEntrySet.from_mp(["Mg", "O"])
        print([e for e in entry_set])


class TestConstrainedFreeEnergyEntrySet(ViseTest):
    def setUp(self) -> None:

        mg = FreeEnergyEntry("Mg", total_energy=-1)
        ca = FreeEnergyEntry("Ca", total_energy=-1)
        mgo = FreeEnergyEntry("MgO", total_energy=-12)
        cao = FreeEnergyEntry("CaO", total_energy=-12)
        o2 = FreeEnergyEntry("O2", total_energy=-2)
        mgcao = FreeEnergyEntry("MgCaO", total_energy=-33)
        self.entry_set = FreeEnergyEntrySet([mg, ca, mgo, cao, o2, mgcao])

    def test_(self):
        entry_set = ConstrainedFreeEnergyEntrySet.from_entry_set(self.entry_set,
                                                                 {"Ca": -2})
        print(entry_set.entries)
        print(entry_set.original_entries)

