# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from math import prod
from pathlib import Path
from typing import List

import fire
import numpy as np
from pydefect.cli.vasp.make_defect_charge_info import make_spin_charges
from pymatgen.io.vasp import Chgcar, Poscar, VolumetricData
from scipy.constants import physical_constants as pc
from vise.analyzer.vesta.vesta_file import add_density
from vise.util.logger import get_logger

logger = get_logger(__name__)



_minor = 1e-3
default_border_fracs = [0.1, 0.5, 0.8]
atomic_unit_to_angstrom = pc["atomic unit of length"][0] * 10 ** 10


def write_little_weight_vol_data(volumetric_data: VolumetricData,
                                 filename: Path,
                                 border_fracs: List[float] = None):
    data = np.zeros(prod(volumetric_data.dim), dtype=int)
    normalized_values = (volumetric_data.data["total"]
                         / np.max(volumetric_data.data["total"])) + _minor

    border_fracs = border_fracs or default_border_fracs
    for border in border_fracs:
        # transpose needed as vasp is based on column measure (Fortran)
        data += (normalized_values > border).T.flatten()

    lines = [Poscar(volumetric_data.structure).get_string()]
    lines.append(" ".join([str(d) for d in volumetric_data.dim]))
    lines.append(" ".join(data.astype(str)))
    filename.write_text("\n".join(lines))


def calc_isurfs(volume, border_fracs):
    isurfs = np.array(border_fracs) / volume * len(border_fracs)
    return np.round(isurfs, 5).tolist()


def add_little_weight_vol_to_vesta(volumetric_file: str,
                                   vesta_file: Path = None,
                                   to_vesta_file: Path = None,
                                   isurfs: List[float] = None):
    add_density(vesta_file, to_vesta_file, isurfs, volumetric_file)
