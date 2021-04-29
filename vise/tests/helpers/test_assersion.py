# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Tuple

import pytest
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


@dataclass
class Data:
    x: float


@dataclass
class TestData:
    a: float = None
    b: List[float] = None
    c: List[Tuple[float]] = None
    d: Data = None


def test_assert_dataclass_almost_equal_float():
    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), 1.0)

    assert_dataclass_almost_equal(TestData(1.0), TestData(1.0))

    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), TestData(2.0))

    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), TestData(1.01), digit=2)

    assert_dataclass_almost_equal(TestData(1.0), TestData(1.01), digit=1)


def test_assert_dataclass_almost_equal_iterable():
    assert_dataclass_almost_equal(TestData(b=[1.0]), TestData(b=[1.0]))

    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(b=[1.0]), TestData(b=[2.0]))


def test_assert_dataclass_almost_equal_nested_iterable():
    assert_dataclass_almost_equal(TestData(c=[(1.0,)]), TestData(c=[(1.0,)]))
    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(c=[(1.0,)]),
                                      TestData(c=[(2.0,)]), digit=2)


def test_assert_dataclass_almost_equal_dataclass_attribute():
    assert_dataclass_almost_equal(TestData(d=Data(1.0)), TestData(d=Data(1.0)))
    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(d=Data(1.0)), TestData(d=Data(1.01)), digit=2)

    assert_dataclass_almost_equal(TestData(d=Data(1.0)), TestData(d=Data(1.001)), digit=2)


