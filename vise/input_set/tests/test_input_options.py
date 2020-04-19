# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

import pytest

from pymatgen.core.structure import Structure

from vise.input_set.input_options import InputOptions, ViseInputOptionsError
from vise.input_set.task import Task
from vise.input_set.xc import Xc


@pytest.fixture()
def sc_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    results = Structure(lattice=lattice, species=["H", "Li"], coords=coords)

    return results


def test_no_options(sc_structure):
    opts = InputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe)

    structure_kpoints_keys = set(opts.structure_kpoints_options.keys())
    assert structure_kpoints_keys == {'initial_structure', 'task'}

    potcar_keys = set(opts.potcar_options.keys())
    assert potcar_keys == {"xc"}

    incar_settings_keys = set(opts.incar_settings_options.keys())
    assert incar_settings_keys == {'xc', 'task'}


def test_initial_structures(sc_structure):
    opts = InputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe)
    assert str(opts.initial_structure.composition) == "H1 Li1"


def test_raise_error(sc_structure):
    with pytest.raises(ViseInputOptionsError):
        InputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                     fake_option="")


def test_structure_kpoints_options(sc_structure):
    opts = InputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                        kpt_density=100)

    assert opts.structure_kpoints_options["kpt_density"] == 100


def test_potcar_options(sc_structure):
    opts = InputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                        override_potcar_set={"A": "B"})

    assert opts.potcar_options["override_potcar_set"] == {"A": "B"}


def test_incar_settings_options(sc_structure):
    opts = InputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                        charge=100.0)

    assert opts.incar_settings_options["charge"] == 100.0


@pytest.fixture()
def input_options():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    results = Structure(lattice=lattice, species=["H", "Li"], coords=coords)

    return results


# def test_integration(sc_structure):
#      input_set = ViseInputSet(sc_structure, task=Task.structure_opt, xc=Xc.pbe)
#      input_set.create_input_files(Path.cwd())



