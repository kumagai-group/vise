# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pymatgen.core import Element
from vise.cli.main_util_functions import make_atom_poscars


def test_make_atom_poscars(mocker):
    args = Namespace(dirname=Path("a"), elements=[Element.H, Element.He])
    mock = mocker.patch("vise.cli.main_utils_functions.make_atom_poscar_dirs")
    make_atom_poscars(args)
    mock.assert_called_once_with(Path("a"), [Element.H, Element.He])



