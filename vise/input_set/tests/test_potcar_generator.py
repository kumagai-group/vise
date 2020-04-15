# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pymatgen import Element
from vise.input_set.potcar_generator import generate_potcar, NoPotcarError


@pytest.fixture()
def mg_symbol():
    return [str(Element("Mg"))]


@pytest.fixture()
def max_z_Cm_symbol():
    return [str(Element.from_Z(96))]


@pytest.fixture()
def Bk_symbol():
    return [str(Element.from_Z(97))]


def test_normal_mg(mg_symbol):
    potcar = generate_potcar(mg_symbol)
    expected = "  PAW_PBE Mg 13Apr2007                   "
    assert str(potcar).split("\n")[0] == expected


def test_mp_relax_set_mg(mg_symbol):
    potcar = generate_potcar(mg_symbol, potcar_set_name="mp_relax_set")
    expected = "  PAW_PBE Mg_pv 13Apr2007                "
    assert str(potcar).split("\n")[0] == expected


def test_gw(mg_symbol):
    potcar = generate_potcar(mg_symbol, potcar_set_name="gw")
    expected = "  PAW_PBE Mg_GW 13Apr2007                    "
    assert str(potcar).split("\n")[0] == expected


def test_normal_largest_z(max_z_Cm_symbol):
    potcar = generate_potcar(max_z_Cm_symbol)
    expected = "  PAW_PBE Cm 17Jan2011                   "
    assert str(potcar).split("\n")[0] == expected


def test_not_exist(Bk_symbol):
    with pytest.raises(NoPotcarError):
        generate_potcar(Bk_symbol)


