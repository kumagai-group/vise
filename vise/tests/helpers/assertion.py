# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import dataclasses
import json

from monty.json import MSONable, MontyDecoder
from monty.serialization import loadfn
from collections.abc import Iterable


def assert_msonable(obj):
    assert isinstance(obj, MSONable)
    assert obj.as_dict() == obj.__class__.from_dict(obj.as_dict()).as_dict()
    assert json.loads(obj.to_json(), cls=MontyDecoder)


def assert_json_roundtrip(obj, tmpdir):
    tmpdir.chdir()
    obj.to_json_file("a.json")
    actual = loadfn("a.json").as_dict()
    expected = obj.as_dict()
    for k, v in actual.items():
        try:
            assert v == expected[k]
        except AssertionError:
            print(k)
            print(v)
            print(expected[k])


def assert_dataclass_almost_equal(actual, expected, digit=8):
    assert actual.__class__ == expected.__class__

    for k, v1 in actual.__dict__.items():
        print(f"key {k}")
        v2 = expected.__getattribute__(k)
        assert_attribute_almost_same(v1, v2, digit)


def assert_attribute_almost_same(v1, v2, digit):
    if isinstance(v1, float):
        try:
            assert round(v1, digit) == round(v2, digit)
        except AssertionError:
            print(v1, v2)
            raise
    elif isinstance(v1, tuple) or isinstance(v1, list):
        for vv1, vv2 in zip(v1, v2):
            assert_attribute_almost_same(vv1, vv2, digit)
    elif dataclasses.is_dataclass(v1):
        assert_dataclass_almost_equal(v1, v2, digit)
    else:
        try:
            assert v1 == v2
        except AssertionError:
            print(v1, v2)
            raise
