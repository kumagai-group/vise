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
default_border_fractions = [0.1, 0.5, 0.8]
au_to_angstrom = pc["atomic unit of length"][0] * 10 ** 10  # = 0.529177210903


def write_light_weight_vol_data(volumetric_data: VolumetricData,
                                filename: Path,
                                border_fractions: List[float] = None):
    data = np.zeros(prod(volumetric_data.dim), dtype=int)
    normalized_values = (volumetric_data.data["total"]
                         / np.max(volumetric_data.data["total"])) + _minor

    border_fractions = border_fractions or default_border_fractions
    for border in border_fractions:
        # transpose needed as vasp is based on column measure (Fortran)
        data += (normalized_values > border).T.flatten()

    lines = [Poscar(volumetric_data.structure).get_string(),
             " ".join([str(d) for d in volumetric_data.dim]),
             " ".join(data.astype(str))]
    filename.write_text("\n".join(lines))


def calc_isurfs(border_fractions: List[float], is_chg: bool, volume: float
                ) -> List[float]:
    """Calc ISURFS values used in the VESTA files

    This is valid only for light-weighted volumetric data.
    """
    # Since max value is set to len(border_fracs), isurfs is multiplied by it.
    isurfs = np.array(border_fractions) * len(border_fractions)
    if is_chg:
        # VESTA uses atomic unit in length.
        isurfs /= (volume / (au_to_angstrom ** 3))
    return np.round(isurfs, 5).tolist()

