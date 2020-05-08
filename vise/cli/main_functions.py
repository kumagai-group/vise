# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace

from pymatgen.ext.matproj import MPRester


def get_poscar_from_mp(args: Namespace) -> None:
    s = MPRester().get_structure_by_material_id(args.mpid)
    s.to(fmt="poscar", filename=args.poscar)
