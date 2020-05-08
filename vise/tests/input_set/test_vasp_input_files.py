# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

from pymatgen.core.structure import Structure

from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.task import Task
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.input_set.xc import Xc


def test_integration(tmpdir):
    tmpdir.chdir()
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    structure = Structure(lattice=lattice, species=["Mn", "O"], coords=coords)
    input_options = CategorizedInputOptions(
        structure=structure, task=Task.structure_opt, xc=Xc.pbe)

    input_set = VaspInputFiles(input_options)
    input_set.create_input_files(Path.cwd())




