# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.util.matplotlib import my_formatter, float_to_int_formatter
import matplotlib.pyplot as plt


def test_formatter():
    assert my_formatter(0.0000049, None) == 0
    assert my_formatter(0.0000051, None) == 5.1e-6
    assert my_formatter("a", None) == "a"


def test_plot():
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    axis = plt.gca()
    axis.xaxis.set_major_formatter(float_to_int_formatter)
    plt.show()
