# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.util.testing import ViseTest
from vise.util.logger import get_logger


class GetLoggerTest(ViseTest):
    def test_get_logger(self):
        print(get_logger("a"))

#        self.assertEqual(True, False)

