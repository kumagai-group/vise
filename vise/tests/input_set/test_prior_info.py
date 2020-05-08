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
from vise.defaults import defaults


@pytest.fixture
def nonmagnetic_insulator():
    return PriorInfo(energy_per_atom=-0.5,
                     band_gap=1.2,
                     vbm_cbm=[1.0, 2.0],
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

#
# @pytest.fixture()
# def sc_structure():
#     lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#     coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
#     results = Structure(lattice=lattice, species=["H", "Li"], coords=coords)
#
#     return results
#
#
# @pytest.fixture
# def input_options(sc_structure, tmpdir):
#     prev_calc_dir = tmpdir
#     os.chdir(prev_calc_dir)
#     sc_structure.to(fmt="POSCAR", filename="CONTCAR-relaxed")
#     vasp_files = PriorInfoFromCalcDir(contcar="CONTCAR-relaxed",
#                                       prev_dir_path=Path(prev_calc_dir))
#     return vasp_files.generate_input_options(task=Task.structure_opt,
#                                              xc=Xc.pbe,
#                                              charge=1.0)


def test_get_structure_from_prev_dir_actual_files(test_data_files, tmpdir):
    tmpdir.chdir()
    prior_info = PriorInfoFromCalcDir(prev_dir_path=test_data_files,
                                      vasprun="MnO_uniform_vasprun.xml",
                                      outcar="MnO_uniform_OUTCAR",
                                      file_transfer_type={"test_one_line": "c"})
    prior_info.file_transfers.transfer(Path(tmpdir))

    assert prior_info.energy_per_atom == -8.024678125
    assert pytest.approx(prior_info.band_gap) == 0.4702
    assert prior_info.vbm_cbm == [4.6666, 5.1368]
    assert prior_info.total_magnetization == 5.0000019
    assert tmpdir.join("test_one_line").read() == "test"


"""
TODO
* Construct FileTransfers. 
* Parse all the previous options.
* Consider if parse_prev_calc is changed to a function or not.

DONE

"""


