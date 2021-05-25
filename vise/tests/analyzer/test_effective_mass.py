# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest

from monty.serialization import loadfn
from vise.analyzer.effective_mass import EffectiveMass, eigvals_and_vecs, \
    lowest_eigval_and_vecs


@pytest.fixture
def em():
    return EffectiveMass(p=[[[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]],
                         n=[[[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]],
                         temperature=300,
                         concentrations=[10**18])


def test_str(em):
    expected = """temperature: 300
------------------------------
concentration: 1e+18
p:
-  -  -
1  0  0
0  1  0
0  0  1
-  -  -
n:
-  -  -
2  0  0
0  2  0
0  0  2
-  -  -"""
    assert str(em) == expected


def test_effective_mass_json_file_mixin(em, tmpdir):
    tmpdir.chdir()
    em.to_json_file()
    actual = loadfn("effective_mass.json")
    assert actual.p == em.p


def test_effective_mass(em):
    actual = em.effective_mass("p", 10**18)
    assert actual == em.p[0]
    actual = em.average_mass("p", 10**18)
    assert actual == 1.0
    actual = em.minimum_mass("p", 10**18)
    assert actual == 1.0


def test_lowest_eigval_and_vec(em):
    matrix = np.array([[1.0, 0.2, 0.1], [0.2, 2.0, 0.3], [0.1, 0.3, 3.0]])
    actual_vals, actual_vecs = eigvals_and_vecs(matrix)
    expected_vals = [0.961, 1.945, 3.094]
    vec = np.array([-46.194, 8.6, 1.0])
    norm = np.linalg.norm(vec, ord=2)
    expected_vec0 = vec / norm
    np.testing.assert_array_almost_equal(actual_vals, expected_vals, decimal=3)
    np.testing.assert_array_almost_equal(actual_vecs[0], expected_vec0)

    actual_lowest_eigval, actual_eigvecs = lowest_eigval_and_vecs(matrix)
    expected_lowest_eigval = 0.9606010675519356
    expected_eigvecs = np.array([[-0.98,  0.18,  0.02]])
    assert actual_lowest_eigval == expected_lowest_eigval
    np.testing.assert_array_almost_equal(actual_eigvecs, expected_eigvecs)


