# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os

import pytest
from vise.input_set.datasets.make_atom_vasp_set import make_atom_vasp_set, \
    make_atom_mp_relax_set
from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.xc import Xc

#
# def test(tmpdir, mocker):
#     tmpdir.chdir()
#     a = PotcarSet.normal
#     make_atom_vasp_set(a, Xc.pbesol)
#
#
#
# def test_atom_mp_relax_set(tmpdir):
#     print(tmpdir)
#     tmpdir.chdir()
#     make_atom_mp_relax_set()


"""
TODO
- Make H dir with potcar inside

DONE
"""
