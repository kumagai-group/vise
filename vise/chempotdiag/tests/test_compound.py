#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
from pymatgen.core.periodic_table import Element
from pathlib import Path
import unittest

from vise.chempotdiag.compound import (
    Compound, CompoundsList, DummyCompoundForDiagram)
from vise.chempotdiag.gas import Gas
from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE, MOLECULE_SUFFIX
from vise.util.testing import ViseTest


MP_TEST = True


class TestCompound(ViseTest):

    def setUp(self) -> None:

        self.chempot_dir = self.TEST_FILES_DIR / "chempotdiag"
        self.o2_dir = self.chempot_dir / "O2molecule"
        self.mg_dir = self.chempot_dir / "Mg"
        self.mgo_dir = self.chempot_dir / "MgO"

    def test_from_vasp_files(self):
        o2 = Compound.from_vasp_files(self.o2_dir, "vasprun.xml.finish")
        expected_o2_energy = -10.26404723 / 2
        self.assertAlmostEqual(expected_o2_energy, o2.energy, 5)

    @unittest.skipIf(not MP_TEST, "not test mp_database")
    def test_from_mp_database(self):
        mgo = Compound.from_mp_database(mp_id="mp-1265")
        expected_mgo_energy = -6.335165765
        self.assertAlmostEqual(expected_mgo_energy, mgo.energy, 5)

    def test_composition_vector(self):
        mgo = Compound.from_vasp_files(self.mgo_dir, "vasprun.xml.finish")
        actual = mgo.composition_vector(["Mg", "N"])
        self.assertAlmostEqual(0.5, actual[0])  # Mg
        self.assertAlmostEqual(0.0, actual[1])  # N

    def test_composition_vector(self):
        mg = Compound.from_vasp_files(self.mg_dir, "vasprun.xml.finish")
        mg.standard_energy = 3.0
        actual = mg.standardized_energy
        expected = -3.42786289 / 2 - 3.0
        self.assertAlmostEqual(expected, actual, 5)


class TestDummyCompound(ViseTest):
    def setUp(self) -> None:
        self.chempot_dir = self.TEST_FILES_DIR / "chempotdiag"
        self.mg_dir = self.chempot_dir / "Mg"

    def test_dummy(self):
        d = DummyCompoundForDiagram.construct_boundary(Element("Mg"), -10)

        self.assertEqual(str(d.elements[0]), "Mg")
        self.assertEqual(d.energy, -10)
        self.assertEqual(d.composition[Element("Mg")], -1)
        with self.assertRaises(NotImplementedError):
            DummyCompoundForDiagram.from_vasp_files(path=self.mg_dir)

    def test_from_vasp_files(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        vasprun_paths = [d+VASPRUN_NAME for d in DFT_DIRECTORIES]
        #  from outcar
        cl: CompoundsList = \
            CompoundsList.from_vasp_files(poscar_paths,
                                          outcar_paths)
        o2: Compound
        mg: Compound
        mgo: Compound
        o2, mg, mgo = tuple(cl)
        self.assertEqual(o2.gas, Gas.O2)
        self.assertAlmostEqual(o2.energy, O2_OUTCAR_VAL, 5)
        self.assertAlmostEqual(mg.energy, MG_OUTCAR_VAL, 5)
        self.assertAlmostEqual(mgo.energy, MGO_OUTCAR_VAL, 5)
        self.assertSetEqual(cl.elements, {Element("Mg"), Element("O")})

        #  from vasprun.xml
        cl_vasprun = CompoundsList.from_vasp_files(poscar_paths,
                                                   vasprun_paths,
                                                   fmt="vasprun")
        o2_vasprun: Compound
        mg_vasprun: Compound
        mgo_vasprun: Compound
        o2_vasprun, mg_vasprun, mgo_vasprun = tuple(cl_vasprun)
        self.assertAlmostEqual(o2.energy, o2_vasprun.energy, 5)
        self.assertAlmostEqual(mg.energy, mg_vasprun.energy, 5)

    def test_standardize_energy(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        cl = CompoundsList.from_vasp_files(poscar_paths,
                                           outcar_paths)
        o2, mg, mgo = tuple(cl)
        for temp, pres in [
                           # only consider zero point energy
                           (0, None),
                           # default: REFERENCE_PRESSURE
                           (ROOM_TEMPERATURE, None),
                           # Specify both T and P
                           (ROOM_TEMPERATURE, 1e+6)]:
            if pres:
                arr, elem_en = cl.standard_energy_array(temp, {"O2": pres})
            else:
                if temp:
                    pres = REFERENCE_PRESSURE
                    arr, elem_en = \
                        cl.standard_energy_array(temp,
                                                 {"O2": REFERENCE_PRESSURE})
                else:
                    arr, elem_en = cl.standard_energy_array(temp, pres)

            if temp:
                self.assertAlmostEqual(elem_en[Element("O")],
                                       O2_OUTCAR_VAL +
                                       Gas.O2.energy_shift(temp, pres))
            else:
                self.assertAlmostEqual(elem_en[Element("O")],
                                       O2_OUTCAR_VAL +
                                       Gas.O2.energy_shift(temp, pres))
            self.assertAlmostEqual(elem_en[Element("Mg")], MG_OUTCAR_VAL)
            for c in arr:
                c: Compound
                if c.composition.is_element:
                    self.assertAlmostEqual(0, c.free_energy(temp, pres), 5)
                else:  # c.composition = MgO
                    # eV / atom, not eV, then divided by 2
                    expected = \
                        MGO_OUTCAR_VAL - \
                        (O2_OUTCAR_VAL +
                         o2.gas.energy_shift(temp, pres) +
                         MG_OUTCAR_VAL
                         ) / 2
                    self.assertAlmostEqual(expected,
                                           c.free_energy(temp, pres), 5)

        # impossible to standardize due to lacking of simple O2.
        cl = CompoundsList.from_vasp_files(
            poscar_paths[1:], outcar_paths[1:])
        with self.assertRaises(ValueError):
            cl.standard_energy_array(ROOM_TEMPERATURE, None)

    def test_gas_shift(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        cl = CompoundsList.from_vasp_files(poscar_paths,
                                           outcar_paths)
        for temp, pres in [
            # only consider zero point energy
            (0, None),
            # default: REFERENCE_PRESSURE
            (ROOM_TEMPERATURE, None),
            # Specify both T and P
            (ROOM_TEMPERATURE, 1e+6)]:
            p = pres if pres else REFERENCE_PRESSURE
            actual_shift = cl.gas_energy_shifts(
                temperature=temp,
                pressure={"O2": p})
            self.assertAlmostEqual(Gas.O2.energy_shift(temp, p),
                                   actual_shift[0], 5)
            self.assertAlmostEqual(0, actual_shift[1], 5)  # Mg
            self.assertAlmostEqual(0, actual_shift[2], 5)  # MgO
            actual_free_energies = cl.free_energies(temperature=temp,
                                                    pressure={"O2": p})
            # O2
            self.assertAlmostEqual(Gas.O2.energy_shift(temp, p) + O2_OUTCAR_VAL,
                                   actual_free_energies[0], 5)
            # Mg
            self.assertAlmostEqual(actual_free_energies[1], MG_OUTCAR_VAL, 5)
            # MgO
            self.assertAlmostEqual(actual_free_energies[2], MGO_OUTCAR_VAL, 5)

