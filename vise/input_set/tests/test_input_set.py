# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import pytest

from pymatgen.core.structure import Structure

from vise.input_set.input_options import ClassifiedInputOptions
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.input_set.task import Task
from vise.input_set.xc import Xc


@pytest.fixture()
def input_options():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    structure = Structure(lattice=lattice, species=["H", "Li"], coords=coords)
    return ClassifiedInputOptions(input_structure=structure, task=Task.structure_opt, xc=Xc.pbe)


def test_integration(input_options, mocker):
#    mock = mocker.patch("vise.input_set.input_set.InputOptions")
    input_set = VaspInputFiles(input_options)
    input_set.create_input_files(Path.cwd())




