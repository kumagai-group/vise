# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from vise.database.symmetry import transmat_primitive2standard
from vise.util.testing import ViseTest


class TransmatPrimitive2StandardTest(ViseTest):
    def test(self):
        expected = np.array([[ 1,  0,  0],
                             [ 0,  1,  1],
                             [ 0, -1,  1]], dtype=int)
        actual = transmat_primitive2standard("A")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 1, -1,  0],
                             [ 1,  1,  0],
                             [ 0,  0,  1]], dtype=int)
        actual = transmat_primitive2standard("C")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 1, -1,  0],
                             [ 0,  1, -1],
                             [ 1,  1,  1]], dtype=int)
        actual = transmat_primitive2standard("R")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 0,  1,  1],
                             [ 1,  0,  1],
                             [ 1,  1,  0]], dtype=int)
        actual = transmat_primitive2standard("I")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[-1,  1,  1],
                             [ 1, -1,  1],
                             [ 1,  1, -1]], dtype=int)
        actual = transmat_primitive2standard("F")
        self.assertArrayEqual(expected, actual)
