# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from enum import unique, Enum


@unique
class ExtendedEnum(Enum):
    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value

    @classmethod
    def from_string(cls, s: str):
        for m in cls:
            if m.value == s:
                return m
        raise AttributeError(f"{str(s)} is not proper member for {cls}.\n",
                             f"Supported names:\n {cls.names_string()}")

    @classmethod
    def names_string(cls):
        return ", ".join([e.name for e in cls])

    @classmethod
    def name_list(cls):
        return [e for e in cls]

    @classmethod
    def name_str_list(cls):
        return [str(e) for e in cls]