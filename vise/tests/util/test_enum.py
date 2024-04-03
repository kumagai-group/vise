# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

from vise.util.enum import ExtendedEnum


def test_extended_enum():

    @dataclass
    class Test(ExtendedEnum):
        a: str = "a"

    assert Test("a") == Test.a
    assert Test.values() == ["a"]

