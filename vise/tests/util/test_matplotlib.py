# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.util.testing import ViseTest
from vise.util.matplotlib import my_formatter


class MyFormatterTest(ViseTest):
    def test(self):
        self.assertEqual(0, my_formatter(0.0, None))
        self.assertEqual(0.1, my_formatter(0.1, None))
        self.assertEqual("a", my_formatter("a", None))

