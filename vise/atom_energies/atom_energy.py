# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from monty.serialization import loadfn
from vise.util.enum import ExtendedEnum
from vise.util.logger import get_logger

path = Path(__file__).parent / "datasets"
logger = get_logger(__name__)

mp_energies = loadfn(path / "mp_vasp544_atom_energy.yaml")


class AtomEnergyType(ExtendedEnum):
    pbe = "pbe"
    pbesol = "pbesol"
    hse = "hse"

    @property
    def energies(self):
        if self.pbe:
            return loadfn(path / "vise_pbe_vasp544_atom_energy.yaml")
        elif self.pbesol:
            return loadfn(path / "vise_pbesol_vasp544_atom_energy.yaml")
        elif self.hse:
            return loadfn(path / "vise_hse_vasp544_atom_energy.yaml")
