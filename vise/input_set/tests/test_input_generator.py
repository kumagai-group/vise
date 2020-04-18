# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

import pytest

from pymatgen.core.structure import Structure

from vise.input_set.input_generator import ViseInputSet
from vise.input_set.task import Task
from vise.input_set.xc import Xc


@pytest.fixture()
def sc_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    results = Structure(lattice=lattice, species=["H"], coords=coords)

    return results


def test_integration(sc_structure):
    input_set = ViseInputSet(sc_structure, task=Task.structure_opt, xc=Xc.pbe)
    input_set.create_input_files(Path.cwd())



