#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
from pymatgen.core.periodic_table import Element
from pathlib import Path

from vise.chempotdiag.compound import (
    Compound, CompoundsList, DummyCompoundForDiagram)
from vise.chempotdiag.gas import Gas
from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE, MOLECULE_SUFFIX
from vise.util.testing import ViseTest


class TestCompound(ViseTest):

    def setUp(self) -> None:
        self.FILENAME_2D = self.TEST_FILES_DIR + "energy_2d.txt"

        self.chempot_dir = self.TEST_FILES_DIR / "chempotdiag"
        # For read DFT test. We don't check these files are physically proper.
        self.o2_dir = self.chempot_dir / "O2"
        self.mg_dir = self.chempot_dir / "Mg"
        self.mgo_dir = self.chempot_dir / "MgO"

        # value of energy(sigma->0)
        # eV / atom, not eV, then divided by 2
        self.O2_OUTCAR_VAL = -10.26404723 / 2
        self.MGO_OUTCAR_VAL = -12.52227962 / 2
        self.MG_OUTCAR_VAL = -3.42786289 / 2

    def test_dummy(self):
        d = DummyCompoundForDiagram.construct_boundary(Element("Mg"), -10)

        self.assertEqual(str(d.elements[0]), "Mg")
        self.assertAlmostEqual(d.energy, -10, 5)
        self.assertAlmostEqual(d.composition[Element("Mg")], -1, 5)
        with self.assertRaises(TypeError):
            o = f"{DFT_DIRECTORIES[1]}/{OUTCAR_NAME}"
            p = f"{DFT_DIRECTORIES[1]}/{POSCAR_NAME}"
            DummyCompoundForDiagram.from_vasp_calculation_files(o, p)

    def test_from_vasp_files(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        vasprun_paths = [d+VASPRUN_NAME for d in DFT_DIRECTORIES]
        #  from outcar
        cl: CompoundsList = \
            CompoundsList.from_vasp_calculations_files(poscar_paths,
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
        cl_vasprun = CompoundsList.from_vasp_calculations_files(poscar_paths,
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
        cl = CompoundsList.from_vasp_calculations_files(poscar_paths,
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
                arr, elem_en = cl.energy_standardized_array(temp, {"O2": pres})
            else:
                if temp:
                    pres = REFERENCE_PRESSURE
                    arr, elem_en = \
                        cl.energy_standardized_array(temp,
                                                     {"O2": REFERENCE_PRESSURE})
                else:
                    arr, elem_en = cl.energy_standardized_array(temp, pres)

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
        cl = CompoundsList.from_vasp_calculations_files(
            poscar_paths[1:], outcar_paths[1:])
        with self.assertRaises(ValueError):
            cl.energy_standardized_array(ROOM_TEMPERATURE, None)

    def test_gas_shift(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        cl = CompoundsList.from_vasp_calculations_files(poscar_paths,
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

