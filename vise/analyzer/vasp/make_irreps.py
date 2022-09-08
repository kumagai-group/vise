# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from typing import List

import numpy as np
from pymatgen.io.vasp import Kpoints

from vise.analyzer.irrep import Irreps, Irrep, Character
from vise.error import ViseError
from vise.util.logger import get_logger

logger = get_logger(__name__)

try:
    from irrep.bandstructure import BandStructure
except ImportError:
    logger.warning(f"To find irreps, one needs to install the irrep package.")
    raise


def special_points_from_kpoints(kpoints_filename):
    kpoints = Kpoints.from_file(kpoints_filename)
    special_points, kpt_idx = [], []

    for idx, label in enumerate(kpoints.labels, 1):
        if label in [None, "None"] or label in special_points:
            continue
        special_points.append(label)
        kpt_idx.append(idx)

    special_points = [x.replace("GAMMA", "GM") for x in special_points]

    return special_points, kpt_idx


def make_irreps_from_wavecar(special_point_chars: List[str],
                             kplist: List[int],
                             wavecar_filename: str = "WAVECAR",
                             poscar_filename: str = "POSCAR",
                             expansion_cutoff: float = 50.0,
                             degenthresh: float = 0.01):  # in eV

    bandstr = BandStructure(fWAV=wavecar_filename,
                            fPOS=poscar_filename,
                            Ecut=expansion_cutoff,
                            kplist=np.array(kplist),
                            spinor=False,
                            searchUC=True)
    characters = bandstr.write_characters(degen_thresh=degenthresh,
                                          kpnames=special_point_chars)
    irrep = {}
    for spc, ckpt, kpt in zip(special_point_chars, characters["k-points"],
                              bandstr.kpoints):
        chars = []
        for e, i, d in zip(ckpt["energies"], ckpt["irreps"], ckpt["dimensions"]):
            try:
                irrep_str = find_irrep(i)
            except ViseNoIrrepError:
                irrep_str = "Unknown"
            chars.append(Character(irrep_str, e, d))

        irrep[spc] = Irrep(tuple(kpt.K), chars)

    return Irreps(bandstr.spacegroup.number, irrep)


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
