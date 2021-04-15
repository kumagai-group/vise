# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

from monty.json import MSONable
from monty.serialization import loadfn
from vise.util.mix_in import ToJsonFileMixIn, ToYamlFileMixIn


@dataclass
class TestTest(MSONable, ToJsonFileMixIn, ToYamlFileMixIn):
    a: str = "aaa"

    def to_yaml(self):
        return f"a: {self.a}\n"

    @classmethod
    def from_yaml(cls, filename):
        d = loadfn(filename)
        return cls(**d)


def test_to_json_file_mix_in(tmpdir):
    actual = TestTest()
    tmpdir.chdir()
    actual.to_json_file()
    actual = loadfn("test_test.json")
    assert actual == TestTest()


def test_to_yaml_file_mix_in(tmpdir):
    actual = TestTest()
    tmpdir.chdir()
    print(tmpdir)
    actual.to_yaml_file()
    actual = TestTest.from_yaml("test_test.yaml")
    assert actual == TestTest()

