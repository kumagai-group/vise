# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import  dataclasses
import sys
from difflib import Differ
from pathlib import Path

from monty.json import MSONable
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

        if dataclasses.is_dataclass(obj):
            try:
                assert obj == rounded_obj
            except ValueError: # numpy array cannot be compared.
                pass
        assert obj.as_dict() == rounded_obj.as_dict()

    except AssertionError:
        print(obj)
        print(rounded_obj)
        print("compare items")
        for k, v in obj.__dict__.items():
            rounded = rounded_obj.__dict__[k]
            print(rounded)
            print(k, v, rounded, type(v), type(rounded))
            try:
                print(v == rounded)
            except ValueError:
                print((v == rounded).all())
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
        print("actual text")
        print(actual_text)
        a = actual_text.split("\n")
        b = expected_text.split("\n")
        sys.stdout.writelines(list(Differ().compare(a, b)))
        raise

    if compare_dict:
        try:
            actual = obj.__class__.from_yaml("a.yaml").as_dict()
            expected = obj.as_dict()
            # When as_dict is not implemented, actual becomes None.
            if actual is None:
                raise AttributeError
        except AttributeError:
            from dataclasses import asdict
            actual = asdict(obj.from_yaml("a.yaml"))
            expected = asdict(obj)
        try:
            assert len(actual) == len(expected)
        except (AssertionError, TypeError):
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


def assert_dataclass_almost_equal(actual, expected, digit=8, check_is_subclass=False):
    assert actual.__class__ == expected.__class__

    for k, v1 in actual.__dict__.items():
        try:
            v2 = expected.__getattribute__(k)
            assert_attribute_almost_same(v1, v2, digit, check_is_subclass)
        except AssertionError:
            print(f"Error raised key {k}")
            raise


def assert_attribute_almost_same(v1, v2, digit, check_is_subclass=False):
    try:
        if check_is_subclass:
            assert issubclass(v1.__class__, v2.__class__)
        else:
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
            try:
                assert_attribute_almost_same(vv1, vv2, digit, check_is_subclass)
            except AssertionError:
                print("original keys")
                print("v1", v1)
                print("v2", v2)
                raise
    elif dataclasses.is_dataclass(v1):
        try:
            assert_dataclass_almost_equal(v1, v2, digit, check_is_subclass)
        except AssertionError:
            print(v1, v2)
            raise
    else:
        try:
            if isinstance(v1, np.ndarray):
                try:
                    np.testing.assert_array_almost_equal(v1, v2)
                except AssertionError:
                    print(v1, v2)
                    raise
            else:
                try:
                    assert v1 == v2
                except AssertionError:
                    print(v1, v2)
                    raise
        except AssertionError:
            print(v1, v2)
            raise
