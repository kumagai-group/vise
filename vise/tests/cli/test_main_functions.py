# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
import pytest
from pathlib import Path

from argparse import Namespace
from pymatgen.core.structure import Structure

from vise.cli.main_functions import get_poscar_from_mp, vasp_set

from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.input_set.kpoints_mode import KpointsMode


def test_get_poscar_from_mp(tmpdir):
    args = Namespace(mpid="mp-110", poscar="POSCAR")

    tmpdir.chdir()
    get_poscar_from_mp(args)
    file = tmpdir.join("POSCAR")
    expected = """Mg1
1.0
2.922478 0.000000 -1.033252
-1.461239 2.530940 -1.033252
0.000000 0.000000 3.099756
Mg
1
direct
0.000000 0.000000 0.000000 Mg
"""
    assert file.read() == expected


@pytest.fixture
def structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    results = Structure(lattice=lattice, species=["Mn", "Mn"], coords=coords)
    return results


option_args = {"structure": structure,
               "task": Task.structure_opt,
               "xc": Xc.pbe,
               "kpt_density": 1.0,
               "override_potcar_set": {"Mn": "Mn_pv"},
               "charge": 2.0}

default_args = deepcopy(option_args)
default_args.update({"user_incar_settings": None,
                     "prev_dir": None,
                     "options": None,
                     "file_transfer_type": None,
                     "uniform_kpt_mode": False,
                     })


test_data = [
    ({}, {}, {}, {}),
    ({"user_incar_settings": {"key": "value"}},
     {"key": "value"},
     {},
     {}),
    ({"options": {"key": "value"}},
     {},
     {"key": "value"},
     {}),
    ({"uniform_kpt_mode": True},
     {},
     {"kpt_mode": KpointsMode.uniform},
     {}),
    ({"prev_dir": "a", "file_transfer_type": {"file": "c"}},
     {},
     {},
     {"x": "y"})
]


@pytest.mark.parametrize("modified_settings,"
                         "overridden_incar_settings,"
                         "overridden_options_args,"
                         "prior_info", test_data)
def test_user_incar_settings(mocker,
                             modified_settings,
                             overridden_incar_settings,
                             overridden_options_args,
                             prior_info):
    args = deepcopy(default_args)
    args.update(modified_settings)

    prior_info = mocker.patch("vise.cli.main_functions.PriorInfoFromCalcDir")
    options = mocker.patch("vise.cli.main_functions.CategorizedInputOptions")
    vif = mocker.patch("vise.cli.main_functions.VaspInputFiles")

    prior_info.return_value.input_options_kwargs = prior_info

    name_space = Namespace(**args)
    vasp_set(name_space)

    incar_settings = defaults.user_incar_settings
    incar_settings.update(overridden_incar_settings)

    options_args = deepcopy(option_args)
    options_args.update(overridden_options_args)
    options_args.update(prior_info)

    vif.assert_called_once_with(options.return_value, incar_settings)
    options.assert_called_once_with(**options_args)


