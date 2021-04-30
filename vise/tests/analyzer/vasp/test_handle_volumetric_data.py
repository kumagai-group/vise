# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from pathlib import Path

from vise.analyzer.vasp.handle_volumetric_data import \
    light_weight_vol_text
from pymatgen import Structure, Lattice
from pymatgen.io.vasp import Chgcar
import numpy as np


def test_light_weight_vol_text(tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    structure = Structure(Lattice.cubic(1), species=["O"], coords=[[0]*3])
    chgcar = Chgcar(structure, data={"total": np.array([[[0.0, 0.1, 0.2, -1e-5]]])})
    actual = light_weight_vol_text(volumetric_data=chgcar,
                                   border_fractions=[0.45, 0.55])
    expected = """O1
1.0
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
O
1
direct
0.000000 0.000000 0.000000 O

1 1 4
0 1 2 0"""
    assert actual == expected

