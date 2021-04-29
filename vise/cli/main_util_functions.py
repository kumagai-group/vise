# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace

from vise.atom_energies.make_atom_vasp_set import make_atom_poscar_dirs
from vise.util.logger import get_logger

logger = get_logger(__name__)


def make_atom_poscars(args: Namespace) -> None:
    make_atom_poscar_dirs(args.dirname, args.elements)


