# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from numpy import prod
from typing import List

import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Poscar, VolumetricData, Chgcar
from vise.util.logger import get_logger

logger = get_logger(__name__)

_minor = 1e-3
default_border_fractions = [0.1, 0.5, 0.8]


def make_spin_charges(chgcar: Chgcar) -> List[Chgcar]:
    result = [Chgcar(chgcar.structure, {"total": chgcar.spin_data[Spin.up]})]
    if "diff" in chgcar.data:
        result.append(
            Chgcar(chgcar.structure, {"total": chgcar.spin_data[Spin.down]}))
    return result


def light_weight_vol_text(volumetric_data: VolumetricData,
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
    return "\n".join(lines)
    # filename.write_text("\n".join(lines))


