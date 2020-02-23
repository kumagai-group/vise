#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.

from pathlib import Path

from vise.chempotdiag.gas import (
    FundamentalFrequencies, ShomateThermodynamicsFunction, Gas)

from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE
from vise.util.testing import ViseTest

molecules_dir = Path(__file__).parent / ".." / "molecules"


class TestFundamentalFrequencies(ViseTest):

    def test_from_yaml(self):
        nh3_yaml = molecules_dir / "NH3" / "fundamental_frequencies.yaml"

        nh3 = FundamentalFrequencies.from_yaml(nh3_yaml)
        expected = [3506, 1022, 3577, 1691]
        actual = [f.wave_number for f in nh3.frequencies]
        self.assertListEqual(expected, actual)

        expected = [1, 1, 2, 2]
        actual = [f.degeneration for f in nh3.frequencies]
        self.assertListEqual(expected, actual)


class TestShomateThermodynamicsFunction(ViseTest):

    def setUp(self):
        o2_file = molecules_dir / "O2" / "shomate_nist.dat"
        no2_file = molecules_dir / "NO2" / "shomate_nist.dat"
        self.o2 = ShomateThermodynamicsFunction.from_nist_table(o2_file)
        self.no2 = ShomateThermodynamicsFunction.from_nist_table(no2_file)

    def test_heat_capacity(self):
        """
        Reference:
        Malcom W. Chase Jr. NIST-JANAF Thermodynamical Tables Fourth Edition,
        J. Phys. Chem. Ref. Data, Monograph 9
        https://srd.nist.gov/JPCRD/jpcrdM9.pdf
        """
        self.assertAlmostEqual(self.o2.heat_capacity(ROOM_TEMPERATURE),
                               29.376, 1)
        self.assertAlmostEqual(self.o2.heat_capacity(400), 30.106, 1)
        self.assertAlmostEqual(self.o2.heat_capacity(1500), 36.544, 1)
        self.assertAlmostEqual(self.o2.heat_capacity(5500), 43.426, 1)

    def test_standard_enthalpy(self):
        """
        Expected value:
        Malcom W. Chase Jr. NIST-JANAF Thermodynamical Tables Fourth Edition,
        J. Phys. Chem. Ref. Data, Monograph 9
        https://srd.nist.gov/JPCRD/jpcrdM9.pdf
        """
        self.assertAlmostEqual(self.o2.standard_enthalpy(ROOM_TEMPERATURE),
                               0, 1)
        self.assertAlmostEqual(self.o2.standard_enthalpy(400), 3.025 / 1000, 1)
        self.assertAlmostEqual(self.o2.standard_enthalpy(1500), 40.599 / 1000, 1)
        self.assertAlmostEqual(self.o2.standard_enthalpy(5500),
                               202.267 / 1000, 1)

    def test_standard_entropy(self):
        """
        Expected value:
        Malcom W. Chase Jr. NIST-JANAF Thermodynamical Tables Fourth Edition,
        J. Phys. Chem. Ref. Data, Monograph 9
        https://srd.nist.gov/JPCRD/jpcrdM9.pdf
        """
        self.assertAlmostEqual(self.o2.standard_entropy(ROOM_TEMPERATURE),
                               205.147, 1)
        self.assertAlmostEqual(self.o2.standard_entropy(400), 213.871, 1)
        self.assertAlmostEqual(self.o2.standard_entropy(1500), 258.068, 1)
        self.assertAlmostEqual(self.o2.standard_entropy(5500), 309.639, 1)

    def test_pressure_term(self):
        self.assertAlmostEqual(self.o2.r_ln_p_p0(REFERENCE_PRESSURE), 0)

    def test_zero_point(self):
        self.assertAlmostEqual(self.o2.enthalpy_zero, -8.683 / 1000)


class TestGas(ViseTest):

    def setUp(self):
        self.o2 = Gas.O2

    def test_read_properties(self):
        self.assertEqual(self.o2.properties["incar"]["NUPDOWN"], 2)

    def test_str(self):
        self.assertEqual(str(self.o2), "O2")

    def test_thermodynamics(self):
        """
        Just check Gas class's thermodynamical quantity is equal to
        ShomateThermodynamicsFunction.
        In TestShomateThermodynamicalFunction, these values are compared to
        those from NIST=JANAF Thermodynamical Tables.

        """
        o2_file = molecules_dir / "O2" / "shomate_nist.dat"
        o2_shomate = ShomateThermodynamicsFunction.from_nist_table(o2_file)

        self.assertEqual(self.o2.heat_capacity(400),
                         o2_shomate.heat_capacity(400))
        self.assertEqual(self.o2.standard_enthalpy(400),
                         o2_shomate.standard_enthalpy(400))
        self.assertEqual(self.o2.standard_entropy(400),
                         o2_shomate.standard_entropy(400))

    def test_temperature_range(self):
        self.assertEqual(self.o2.min_temperature, 100)
        self.assertEqual(self.o2.max_temperature, 6000)
        self.assertEqual(self.o2.temperature_range, (100, 6000))

    def test_zero_point(self):
        # Expected value is script of Prof. Kumagai
        self.assertAlmostEqual(Gas.O2.zero_point_vibrational_energy, 0.098/2, 2)
        self.assertAlmostEqual(Gas.N2.zero_point_vibrational_energy, 0.146/2, 2)

    def test_energy_shift(self):
        
        for t in range(0, 1001, 100):
            for p in [REFERENCE_PRESSURE/100, REFERENCE_PRESSURE]:
                if t == 0:
                    expected = self.o2.zero_point_vibrational_energy
                else:
                    expected = self.o2.zero_point_vibrational_energy + \
                               self.o2.vib_rot_term(t) + \
                               t * self.o2.r_ln_p_p0_ev_per_atom(p)

                actual = self.o2.energy_shift(temperature=t,
                                              pressure=p)
                self.assertAlmostEqual(expected, actual, 5)

