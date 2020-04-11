# -*- coding: utf-8 -*-

from pathlib import Path

from monty.serialization import loadfn

unoccupied_bands = loadfn(Path(__file__).parent / "unoccupied_bands.yaml")
magmom = loadfn(Path(__file__).parent / "magmom.yaml")
