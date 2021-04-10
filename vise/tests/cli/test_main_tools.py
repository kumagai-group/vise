# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from vise.cli.main_tools import potcar_str2dict, list2dict


def test_none():
    assert potcar_str2dict(None) == {}


def test_list():
    expected = {"Mg": "Mg_pv", "O": "O_h"}
    actual = potcar_str2dict(["Mg_pv", "O_h"])
    assert actual == expected


def test_str():
    expected = {"Mg": "Mg_pv"}
    actual = potcar_str2dict("Mg_pv")
    assert actual == expected


def test_mutliple_potcars_for_same_element_error():
    with pytest.raises(ValueError):
        potcar_str2dict(["Mg_pv", "Mg"])


def test_incorrect_element_potcar_error():
    with pytest.raises(ValueError):
        potcar_str2dict(["MgHH", "Mg"])


@pytest.fixture
def key_candidates():
    return ["ENCUT", "MAGMOM", "LWAVE", "EFIELD_PEAD"]


def test_dict2list(key_candidates):
    flattened_list = ["ENCUT", "500", "MAGMOM", "4", "4.0", "LWAVE", "F"]
    actual = list2dict(flattened_list, key_candidates)
    expected = {"ENCUT": 500, "MAGMOM": [4, 4.0], "LWAVE": False}
    assert actual == expected


def test_dict2list_2(key_candidates):
    flattened_list = ["EFIELD_PEAD", "0.1", "0.1", "0.1"]
    actual = list2dict(flattened_list, key_candidates)
    expected = {"EFIELD_PEAD": [0.1, 0.1, 0.1]}
    assert actual == expected


def test_fail(key_candidates):
    flattened_list = ["ENCAT", "500"]
    with pytest.raises(ValueError):
        list2dict(flattened_list, key_candidates)


def test_fail2(key_candidates):
    flattened_list = ["ENCUT", "500", "MAGMOM"]
    with pytest.raises(ValueError):
        list2dict(flattened_list, key_candidates)


