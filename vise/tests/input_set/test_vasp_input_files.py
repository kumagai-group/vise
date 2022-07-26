# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

import pytest
from pymatgen.core.structure import Structure

from vise import __version__
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.task import Task
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.input_set.xc import Xc


@pytest.fixture
def vasp_input_files():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    structure = Structure(lattice=lattice, species=["Mn", "O"], coords=coords)
    input_options = CategorizedInputOptions(
        structure=structure, task=Task.structure_opt, xc=Xc.pbe, charge=1)

    return VaspInputFiles(input_options, overridden_incar_settings={"NSW": 2})


def test_integration(tmpdir, vasp_input_files):
    tmpdir.chdir()
    vasp_input_files.create_input_files(Path.cwd())


def test_vise_log(vasp_input_files):
    expected = {"version": __version__,
                "task": str(Task.structure_opt),
                "xc": str(Xc.pbe),
                "input_options": {"charge": 1},
                "user_incar_settings": {"NSW": 2}}
    assert vasp_input_files.vise_log == expected



