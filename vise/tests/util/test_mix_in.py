# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from pathlib import Path

from monty.json import MSONable
from monty.serialization import loadfn
from vise.util.mix_in import ToJsonFileMixIn, ToYamlFileMixIn, ToCsvFileMixIn


@dataclass
class TestTest(MSONable, ToJsonFileMixIn, ToCsvFileMixIn, ToYamlFileMixIn):
    a: str = "aaa"

    def as_dict(self):
        return {"a": self.a}

    @classmethod
    def from_yaml(cls, filename):
        d = loadfn(filename)
        return cls(**d)

    @property
    def csv_column_names(self):
        return "a"

    @property
    def csv_data(self):
        return [["aaa"]]


def test_to_json_file_mix_in(tmpdir):
    actual = TestTest()
    tmpdir.chdir()
    actual.to_json_file()
    actual = loadfn("test_test.json")
    assert actual == TestTest()


def test_to_csv_file_mix_in(tmpdir):
    test = TestTest()
    tmpdir.chdir()
    test.to_csv_file()
    actual = Path("test_test.csv").read_text()
    expected = """a
aaa
"""
    assert actual == expected


def test_to_yaml_file_mix_in(tmpdir):
    actual = TestTest()
    tmpdir.chdir()
    print(tmpdir)
    actual.to_yaml_file()
    actual = TestTest.from_yaml("test_test.yaml")
    assert actual == TestTest()

