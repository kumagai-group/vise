# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from monty.serialization import loadfn
from vise.util.enum import ExtendedEnum
from vise.util.logger import get_logger

path = Path(__file__).parent

mp_energies = loadfn(path / "mp_atom_energy.yaml")


