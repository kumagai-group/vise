#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

from vise.input_set.incar import ViseIncar
from vise.input_set.make_kpoints import MakeKpoints

from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.inputs import VaspInput


def create_atom_dirs(task, xc):

    lattice = Lattice([[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    for z in range(1, 92):
        elem = str(Element.from_Z(z))
        d = Path(elem)
        s = Structure(lattice, [elem], [[0.5, 0.5, 0.5]])

        mp_set = MPStaticSet(structure=s)
        incar = ViseIncar.from_dict(mp_set.incar.as_dict())
        potcar = mp_set.potcar
        make_kpoints = MakeKpoints(mode="manual_set", structure=s)
        make_kpoints.make_kpoints()
        kpoints = make_kpoints.kpoints
        poscar = Poscar(s)
        vasp_input = VaspInput(incar, kpoints, poscar, potcar)

        d.mkdir()
        vasp_input.write_input(d)



