# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise import __version__
from vise.defaults import defaults
from vise.input_set.input_options import (
    CategorizedInputOptions, ViseInputOptionsError)
from vise.input_set.task import Task
from vise.input_set.vise_log import ViseLog
from vise.input_set.xc import Xc


def test_no_options(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe)

    structure_kpoints_keys = set(opts.structure_kpoints_options.keys())
    assert structure_kpoints_keys == {"kpt_density", "initial_structure",
                                      "task"}

    potcar_keys = set(opts.potcar_options.keys())
    assert potcar_keys == {"xc"}

    incar_settings_keys = set(opts.incar_settings_options.keys())
    assert incar_settings_keys == {'xc', 'task'}


def test_initial_structures(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe)
    assert str(opts.initial_structure.composition) == "H1"


def test_raise_error(sc_structure):
    with pytest.raises(ViseInputOptionsError):
        CategorizedInputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                                fake_option="")


def test_structure_kpoints_options(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt,
                                   xc=Xc.pbe, kpt_density=100)

    assert opts.structure_kpoints_options["kpt_density"] == 100


def test_potcar_options(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                                   overridden_potcar={"A": "B"})

    assert opts.potcar_options["overridden_potcar"] == {"A": "B"}


def test_incar_settings_options(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                                   charge=100.0)

    assert opts.incar_settings_options["charge"] == 100.0


def test_initial_structure(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt, xc=Xc.pbe,
                                   charge=100.0)

    assert opts.initial_structure == sc_structure


def test_insulator_kpt_density(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt,
                                   xc=Xc.pbe, vbm_cbm=[0, 0.5])
    actual = opts.structure_kpoints_options["kpt_density"]
    expected = defaults.insulator_kpoint_density
    assert actual == expected

    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt,
                                   xc=Xc.pbe, band_gap=1.0)
    actual = opts.structure_kpoints_options["kpt_density"]
    expected = defaults.insulator_kpoint_density
    assert actual == expected

    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt,
                                   xc=Xc.pbe, vbm_cbm=[0, 0.5],
                                   kpt_density=100)
    assert opts.structure_kpoints_options["kpt_density"] == 100


def test_parameter_dict(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.structure_opt,
                                   xc=Xc.pbe, vbm_cbm=[0, 0.5])
    assert opts.vise_log == ViseLog(
        version=__version__, task=Task.structure_opt, xc=Xc.pbe,
        input_options={"vbm_cbm": [0, 0.5],
                       "kpt_density": defaults.insulator_kpoint_density})


def test_defect_kpt_density(sc_structure):
    opts = CategorizedInputOptions(sc_structure, task=Task.defect, xc=Xc.pbe)
    actual = opts.structure_kpoints_options["kpt_density"]
    expected = defaults.defect_kpoint_density
    assert actual == expected

