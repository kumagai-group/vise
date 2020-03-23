# -*- coding: utf-8 -*-
from pathlib import Path

from vise.chempotdiag.gas import Gas

from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE
from vise.util.testing import ViseTest

molecules_dir = Path(__file__).parent / ".." / "molecules"


class TestGas(ViseTest):

    def setUp(self):
        self.o2 = Gas("O2", temperature=300)
        self.nh3 = Gas("NH3", temperature=1000)

    def test_o2 (self):
        self.assertAlmostEqual(self.o2.trans_free_energy, -0.40862436734991214)
        self.assertAlmostEqual(self.o2.rot_free_energy, -0.1107263561129291)
        self.assertAlmostEqual(self.o2.zero_point_vibrational_energy,
                               0.09720351919699763)
        self.assertAlmostEqual(self.o2.vib_free_energy, -1.401901140122062e-05)
        self.assertAlmostEqual(self.o2.spin_free_energy, -0.02840132465202343)
        self.assertAlmostEqual(self.o2.energy_shift, -0.2252812739646341)

    def test_nh3(self):
        self.assertAlmostEqual(self.nh3.rot_free_energy, -0.5283431627639921)

    def test_assert_error(self):
        with self.assertRaises(ValueError):
            low_temp_nh3 = Gas("NH3", temperature=100)

