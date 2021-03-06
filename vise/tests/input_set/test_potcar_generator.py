# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen.core import Element

from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.potcar_generator import generate_potcar, ViseNoPotcarError
from vise.input_set.xc import Xc


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
    potcar = generate_potcar(mg_symbol, Xc.pbe)
    expected = "  PAW_PBE Mg 13Apr2007".split()
    assert str(potcar).split("\n")[0].split() == expected


def test_mp_relax_set_mg(mg_symbol):
    potcar = generate_potcar(mg_symbol,
                             Xc.pbe,
                             potcar_set=PotcarSet.mp_relax_set)
    expected = "  PAW_PBE Mg_pv 13Apr2007".split()
    assert str(potcar).split("\n")[0].split() == expected


def test_gw(mg_symbol):
    potcar = generate_potcar(mg_symbol, Xc.pbe, potcar_set=PotcarSet.gw)
    expected = "  PAW_PBE Mg_GW 13Apr2007                    "
    assert str(potcar).split("\n")[0] == expected


def test_normal_largest_z(max_z_Cm_symbol):
    potcar = generate_potcar(max_z_Cm_symbol, Xc.pbe)
    expected = "  PAW_PBE Cm 17Jan2011                   "
    assert str(potcar).split("\n")[0] == expected


def test_not_exist(Bk_symbol):
    with pytest.raises(ViseNoPotcarError):
        generate_potcar(Bk_symbol, Xc.pbe)


