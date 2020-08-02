# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.util.logger import get_logger


def test_get_logger():
    assert str(get_logger("a")) == "<Logger a (DEBUG)>"
