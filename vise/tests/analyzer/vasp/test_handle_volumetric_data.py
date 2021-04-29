# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from pathlib import Path

from vise.analyzer.vasp.handle_volumetric_data import write_little_weight_vol_data
from pymatgen import Structure, Lattice
from pymatgen.io.vasp import Chgcar
import numpy as np


def test_write_little_weight_vol_data(tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    structure = Structure(Lattice.cubic(1), species=["O"], coords=[[0]*3])
    chgcar = Chgcar(structure, data={"total": np.array([[[0.0, 0.1, 0.2, -1e-5]]])})
    write_little_weight_vol_data(volumetric_data=chgcar,
                                 filename=Path("CHGCAR"),
                                 border_fracs=[0.45, 0.55])

    actual = Chgcar.from_file("CHGCAR")
    expected = Chgcar(structure, data={"total": np.array([[[0, 1, 2, 0]]])})
    np.testing.assert_almost_equal(actual.data["total"], expected.data["total"])

