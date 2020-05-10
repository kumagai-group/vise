# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.input_set.datasets.potcar_set import PotcarSet


def test_potcar_list_normal():
    potcar_list = PotcarSet.normal.potcar_dict()
    assert potcar_list["H"] == "H"


def test_potcar_list_gw():
    potcar_list = PotcarSet.gw.potcar_dict()
    assert potcar_list["Li"] == "Li_GW"


def test_potcar_list_non():
    potcar_list = PotcarSet.mp_relax_set.potcar_dict()
    assert potcar_list["Fr"] is None
    assert potcar_list["Sr"] == "Sr_sv"


def test_potcar_set():
    assert PotcarSet.normal.overridden_potcar_dict()["Li"] == "Li"


def test_potcar_set_override():
    actual = PotcarSet.normal.overridden_potcar_dict({"Li": "Li_test"})["Li"]
    assert actual == "Li_test"