# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.xc import Xc


def test_xc_str():
    assert "pbe" == str(Xc.pbe)


def test_xc_from_string():
    assert Xc.from_string("pbe") == Xc.pbe
    assert Xc.from_string("perdew-zunger81") == Xc.lda


def test_xc_from_string_raise_error():
    with pytest.raises(AttributeError):
        Xc.from_string("not_exist")


def test_potcar_functional():
    assert Xc.scan.potcar_functional == "PBE_54"
    assert Xc.lda.potcar_functional == "LDA"


def test_xc_is_lda_or_gga():
    assert Xc.pbe.is_lda_or_gga
    assert Xc.hse.is_lda_or_gga is False


def test_xc_is_hybrid_functional():
    assert Xc.hse.is_hybrid_functional is True


def test_xc_is_local_or_semilocal():
    assert Xc.pbe.is_local_or_semilocal is True
    assert Xc.scan.is_local_or_semilocal is True
    assert Xc.hse.is_local_or_semilocal is False


def test_xc_is_nonlocal():
    assert Xc.scan.is_nonlocal is False
    assert Xc.hse.is_nonlocal is True


def test_incar_algo():
    assert Xc.hse.incar_algo == "Damped"
    assert Xc.pbe.incar_algo == "Normal"


def test_incar_lwave():
    assert Xc.hse.incar_lwave is True
    assert Xc.pbe.incar_lwave is False


def test_gga_optional():
    assert Xc.pbesol.incar_gga_optional == "PS"
    assert Xc.pbe.incar_gga_optional is None


def test_metagga_optional():
    assert Xc.scan.incar_metagga_optional == "SCAN"
    assert Xc.pbe.incar_metagga_optional is None
