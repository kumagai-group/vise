# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
import tempfile
from pathlib import Path

import pytest
from pymatgen.core.structure import Structure

from vise.input_set.prior_info import PriorInfo, PriorInfoFromCalcDir
from vise.input_set.task import Task
from vise.input_set.xc import Xc


@pytest.fixture
def nonmagnetic_insulator():
    return PriorInfo(energy_per_atom=-0.5,
                     band_gap=1.2,
                     total_magnetization=0.000001,
                     data_source="Materials Project",
                     is_cluster=False,
                     magnetization_criterion=0.001,
                     band_gap_criterion=0.1,
                     incar={"NUPDOWN": 2})


def test_round_trip_dict(nonmagnetic_insulator):
    d = nonmagnetic_insulator.as_dict()
    prior_info_from_dict = PriorInfo.from_dict(d)
    assert prior_info_from_dict.as_dict() == d


def test_round_trip_yaml(nonmagnetic_insulator):
    tmp_file = tempfile.NamedTemporaryFile()
    nonmagnetic_insulator.dump_yaml(tmp_file.name)
    prior_info_from_json = PriorInfo.load_yaml(tmp_file.name)
    assert prior_info_from_json.as_dict() == nonmagnetic_insulator.as_dict()


def test_round_trip_json(nonmagnetic_insulator):
    tmp_file = tempfile.NamedTemporaryFile()
    nonmagnetic_insulator.dump_json(tmp_file.name)
    prior_info_from_json = PriorInfo.load_json(tmp_file.name)
    assert prior_info_from_json.as_dict() == nonmagnetic_insulator.as_dict()


def test_properties(nonmagnetic_insulator):
    assert nonmagnetic_insulator.is_magnetic is False
    assert nonmagnetic_insulator.has_band_gap is True
    assert nonmagnetic_insulator.is_metal is False


"""
TODO
* Allow to set directory Path.
* Construct FileTransfers. 
* Parse all the previous options.
* Consider if parse_prev_calc is changed to a function or not.

DONE
* Extract fixture from the test.
* Parse the cwd and get a structure from an existing POSCAR
* Parse previous directory and get structure
* Allow the different POSCAR name.
* Parse cwd, get band gap and set vbm_cbm
* Parse cwd, get magnetic moment and set is_magnetization.
"""


@pytest.fixture()
def sc_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    results = Structure(lattice=lattice, species=["H", "Li"], coords=coords)

    return results

@pytest.fixture
def input_options(sc_structure):
    with tempfile.TemporaryDirectory() as prev_calc_dir:
        os.chdir(prev_calc_dir)
        sc_structure.to(fmt="POSCAR", filename="CONTCAR-relaxed")
        vasp_files = PriorInfoFromCalcDir(contcar="CONTCAR-relaxed",
                                          prev_dir_path=Path(prev_calc_dir))
        return vasp_files.generate_input_options(task=Task.structure_opt,
                                                 xc=Xc.pbe,
                                                 charge=1.0)


def test_get_structure_from_prev_dir(input_options, sc_structure, mocker):
    # mock = mocker.patch("vise.input_set.input_options.vbm_cbm_from_vasprun")
    # mock.return_value = [0.0, 1.0]
    # mock_mag = mocker.patch("vise.input_set.input_options.is_magnetic_from_outcar")
    # mock_mag.return_value = True

    assert input_options.initial_structure == sc_structure
    assert input_options.incar_settings_options["charge"] == 1.0
    assert input_options.incar_settings_options["vbm_cbm"] == [0.0, 1.0]
    assert input_options.incar_settings_options["is_magnetization"] is True


# def test_get_structure_from_prev_dir(sc_structure):
#         with tempfile.TemporaryDirectory() as tmp_from:
#         os.chdir(tmp_from)
#         sc_structure.to(fmt="POSCAR", filename="CONTCAR-relaxed")
#         opts = CategorizedInputOptions.parse_prev_calc(
#             task=Task.structure_opt,
#             xc=Xc.pbe,
#             contcar="CONTCAR-relaxed",
#             charge=1.0)

# assert opts.initial_structure == sc_structure
# assert opts.incar_settings_options["charge"] == 1.0

