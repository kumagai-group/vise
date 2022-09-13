# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from copy import copy
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pytest
from monty.json import MSONable
from pymatgen.core import IStructure, Lattice
from vise.tests.helpers.assertion import assert_dataclass_almost_equal, \
    assert_structure_almost_same, assert_msonable


@dataclass
class Data:
    x: float


def test_assert_msonable_round_trip_cover_all_attributes():
    @dataclass
    class TestMSONableData(MSONable):
        x: float = None
        y: str = None

        def as_dict(self):
            return {"x": self.x}

        @classmethod
        def from_dict(cls, d):
            kwargs = copy(d)
            kwargs.pop("@module", None)
            kwargs.pop("@class", None)
            return cls(**d)

    with pytest.raises(AssertionError):
        assert_msonable(TestMSONableData(1.0, "a"))


@dataclass
class TestData:
    a: float = None
    b: List[float] = None
    c: List[Tuple[float]] = None
    d: Data = None


def test_assert_structure_almost_same():
    coords = [[0.0, 0.0, 0.0]]
    s1 = IStructure(lattice=Lattice.cubic(1.0), species=["H"], coords=coords)
    coords = [[0.0, 0.0, 0.0]]

    s2 = IStructure(lattice=Lattice.cubic(1.000001), species=["H"],
                    coords=coords)
    assert_structure_almost_same(s1, s2)

    s3 = IStructure(lattice=Lattice.cubic(1.001), species=["H"], coords=coords)
    with pytest.raises(AssertionError):
        assert_structure_almost_same(s1, s3)


def test_assert_dataclass_almost_equal_float():
    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), 1.0)

    assert_dataclass_almost_equal(TestData(1.0), TestData(1.0))

    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), TestData(2.0))

    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), TestData(1.01), digit=2)

    assert_dataclass_almost_equal(TestData(1.0), TestData(1.01), digit=1)


def test_assert_dataclass_almost_equal_w_different_type_attr():
    with pytest.raises(AssertionError):
        assert_dataclass_almost_equal(TestData(1.0), TestData(np.float32(1.0)))


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


