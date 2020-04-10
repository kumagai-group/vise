# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.xc import Xc


def test_xc_str():
    assert "pbe" == str(Xc.pbe)


def test_xc_from_string():
    assert Xc.pbe == Xc.from_string("pbe")
    assert Xc.lda == Xc.from_string("perdew-zunger81")


def test_xc_from_string_raise_error():
    with pytest.raises(AttributeError):
        Xc.from_string("not_exist")


def test_xc_is_lda_or_gga():
    assert Xc.pbe.is_lda_or_gga
    assert Xc.hse.is_lda_or_gga is False


def test_xc_is_hybrid_functional():
    assert Xc.hse.is_hybrid_functional


def test_xc_is_local_or_semilocal():
    assert Xc.pbe.is_local_or_semilocal
    assert Xc.scan.is_local_or_semilocal
    assert Xc.hse.is_local_or_semilocal is False


def test_xc_is_nonlocal():
    assert Xc.scan.is_nonlocal is False
    assert Xc.hse.is_nonlocal
