# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import tempfile

import pytest

from vise.input_set.prior_info import PriorInfo, prior_info_from_calc_dir


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


def test_get_structure_from_prev_dir_actual_files(test_data_files):
    prior_info = prior_info_from_calc_dir(prev_dir_path=test_data_files,
                                          vasprun="MnO_uniform_vasprun.xml",
                                          outcar="MnO_uniform_OUTCAR")

    assert prior_info.energy_per_atom == -8.024678125
    assert pytest.approx(prior_info.band_gap) == 0.4702
    assert prior_info.vbm_cbm == [4.6666, 5.1368]
    assert prior_info.total_magnetization == 5.0000019




