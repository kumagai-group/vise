# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pathlib import Path

from argparse import Namespace
from pymatgen.core.structure import Structure

from vise.cli.main_functions import get_poscar_from_mp
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc


def test_get_poscar_from_mp(tmpdir):
    args = Namespace(mpid="mp-110", poscar="POSCAR")

    tmpdir.chdir()
    get_poscar_from_mp(args)
    file = tmpdir.join("POSCAR")
    expected = """Mg1
1.0
2.922478 0.000000 -1.033252
-1.461239 2.530940 -1.033252
0.000000 0.000000 3.099756
Mg
1
direct
0.000000 0.000000 0.000000 Mg
"""
    assert file.read() == expected


def test_create_actual_files():
    poscar = """Mg1
1.0
2.922478 0.000000 -1.033252
-1.461239 2.530940 -1.033252
0.000000 0.000000 3.099756
Mg
1
direct
0.000000 0.000000 0.000000 Mg
"""
    structure = Structure.from_str(poscar, fmt="POSCAR")
    options = CategorizedInputOptions(structure, Task(defaults.task), Xc(defaults.xc))
    vif = VaspInputFiles(options)
    vif.create_input_files(Path.cwd())


"""
TODO:
- Create Vasp Set of default set.

"""