# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from typing import List, Tuple

import numpy as np
from pymatgen.io.vasp import Kpoints

from vise.analyzer.irrep import Irreps, Irrep
from vise.error import ViseError
from vise.util.logger import get_logger

logger = get_logger(__name__)

try:
    from irrep.bandstructure import BandStructure
except ImportError:
    logger.warning(f"To find irreps, one needs to install the irrep package.")
    raise


def special_points_from_kpoints(kpoints_filename: str) \
        -> Tuple[List[str], List[int]]:

    kpoints = Kpoints.from_file(kpoints_filename)
    special_points, kpt_indices = [], []

    for idx, label in enumerate(kpoints.labels, 1):
        if label in [None, "None"] or label in special_points:
            continue
        special_points.append(label)
        kpt_indices.append(idx)

    special_points = [x.replace("GAMMA", "GM") for x in special_points]

    return special_points, kpt_indices


def make_irreps_from_wavecar(special_point_characters: List[str],
                             kpt_indices: List[int],
                             wavecar_filename: str = "WAVECAR",
                             poscar_filename: str = "POSCAR",
                             plane_wave_cutoff: float = 50.0,
                             degeneracy_threshold: float = 0.01) -> Irreps:

    bs = BandStructure(fWAV=wavecar_filename,
                       fPOS=poscar_filename,
                       Ecut=plane_wave_cutoff,
                       kplist=np.array(kpt_indices),
                       spinor=False,
                       searchUC=True)
    characters = bs.write_characters(degen_thresh=degeneracy_threshold,
                                     kpnames=special_point_characters)
    irrep_dict = {}
    for character, c_kpt, kpt in \
            zip(special_point_characters, characters["k-points"], bs.kpoints):
        symbols = []
        for irrep in c_kpt["irreps"]:
            try:
                irrep_str = find_irrep(irrep)
            except ViseNoIrrepError:
                irrep_str = "Unknown"
            symbols.append(irrep_str)

        irrep_dict[character] = Irrep(tuple(kpt.K),
                                      symbols,
                                      c_kpt["energies"].tolist(),
                                      c_kpt["dimensions"])

    return Irreps(bs.spacegroup.number, irrep_dict)


class ViseNoIrrepError(ViseError):
    pass


def find_irrep(d: dict, threshold: float = 0.99):
    for k, v in d.items():
        if v[0] > threshold:
            return k
    logger.warning(f"""Any irrep could not be found. 
Threshold: {threshold}
Characters: {d}""")
    raise ViseNoIrrepError
