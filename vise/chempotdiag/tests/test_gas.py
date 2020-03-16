# -*- coding: utf-8 -*-
from pathlib import Path

from vise.chempotdiag.gas import Gas

from vise.config import ROOM_TEMPERATURE, REFERENCE_PRESSURE
from vise.util.testing import ViseTest

molecules_dir = Path(__file__).parent / ".." / "molecules"


class TestGas(ViseTest):

    def setUp(self):
        self.o2 = Gas("O2", temperature=300)
        self.nh3 = Gas("NH3", temperature=100)

    def test_thermodynamics(self):
        """ """
        print("")
        print("o2")
        print(self.o2.trans_free_energy)
        print(self.o2.rot_free_energy())
        print(self.o2.zero_point_vibrational_energy)
        print(self.o2.vib_free_energy)
        print(self.o2.spin_free_energy)
        print(self.o2.energy_shift)
        print("nh3")
        print(self.nh3.trans_free_energy)
        print(self.nh3.rot_free_energy())
        print(self.nh3.zero_point_vibrational_energy)
        print(self.nh3.vib_free_energy)
        print(self.nh3.spin_free_energy)
        print(self.nh3.energy_shift)
        print(self.nh3.chem_pot_shift)

    def test_nh3_2(self):
        for i in range(300, 1001, 10):
            nh3 = Gas("NH3", temperature=i)
            print(i, nh3.chem_pot_shift)
