# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import dataclasses
import json
import sys
from difflib import Differ
from pathlib import Path

from monty.json import MSONable, MontyDecoder
from monty.serialization import loadfn
import numpy as np
from py._path.local import LocalPath
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    ElementComparator
from pymatgen.core import IStructure
from vise.util.mix_in import ToYamlFileMixIn


def assert_structure_almost_same(s1: IStructure, s2: IStructure):
    matcher = StructureMatcher(ltol=0.0001,
                               stol=0.0001,
                               angle_tol=0.0001,
                               primitive_cell=False,
                               scale=False,
                               attempt_supercell=False,
                               allow_subset=False,
                               comparator=ElementComparator())
    assert matcher.fit(s1, s2, skip_structure_reduction=True)


def assert_msonable(obj):
    assert isinstance(obj, MSONable)
    rounded_obj = obj.__class__.from_dict(obj.as_dict())
    try:
        assert obj == rounded_obj
        assert obj.as_dict() == rounded_obj.as_dict()
        assert json.loads(obj.to_json(), cls=MontyDecoder)
    except AssertionError:
        print(obj)
        print(rounded_obj)
        for k, v in obj.__dict__.items():
            print(k, v, type(v))
        for k, v in rounded_obj.__dict__.items():
            print(k, v, type(v))
        raise


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


def assert_yaml_roundtrip(obj: ToYamlFileMixIn,
                          tmpdir: LocalPath,
                          expected_text: str,
                          compare_dict: bool = True,
                          compare_items: bool = True):
    tmpdir.chdir()
    obj.to_yaml_file("a.yaml")
    actual_text = Path("a.yaml").read_text()
    try:
        assert actual_text == expected_text
    except AssertionError:
        print(tmpdir)
        print(actual_text)
        a = actual_text.split("\n")
        b = expected_text.split("\n")
        sys.stdout.writelines(list(Differ().compare(a, b)))
        raise

    if compare_dict:
        try:
            actual = obj.from_yaml("a.yaml").as_dict()
            expected = obj.as_dict()
        except AttributeError:
            from dataclasses import asdict
            actual = asdict(obj.from_yaml("a.yaml"))
            expected = asdict(obj)
        try:
            assert len(actual) == len(expected)
        except AssertionError:
            print(tmpdir)
            raise
    else:
        actual = obj.from_yaml("a.yaml")
        print("actual", actual)

    if compare_items is True:
        for k, v in actual.items():
            try:
                assert v == expected[k]
            except AssertionError:
                print(tmpdir)
                print(f"key: {k}, actual: {v}, expected: {expected[k]}")
                raise


def assert_dataclass_almost_equal(actual, expected, digit=8):
    assert actual.__class__ == expected.__class__

    for k, v1 in actual.__dict__.items():
        try:
            v2 = expected.__getattribute__(k)
            assert_attribute_almost_same(v1, v2, digit)
        except AssertionError:
            print(f"Error raised key {k}")


def assert_attribute_almost_same(v1, v2, digit):
    try:
        assert v1.__class__ == v2.__class__
    except AssertionError:
        print("type(v1), type(v2)")
        print(type(v1), type(v2))
        raise

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
            if isinstance(v1, np.ndarray):
                np.testing.assert_array_almost_equal(v1, v2)
            else:
                assert v1 == v2
        except AssertionError:
            print(v1, v2)
            raise
