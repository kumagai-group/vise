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
    result = set(kpoints.labels)
    result.discard(None)
    result.discard("None")
    return {s.replace("GAMMA", "GM") for s in result}


def make_irreps_from_wavecar(wavecar_filename: str = "WAVECAR",
                             poscar_filename: str = "POSCAR",
                             expansion_cutoff: float = 50.0,
                             special_point_chars: List[str] = None,
                             degenthresh=0.01):  # in eV

    bandstr = BandStructure(fWAV=wavecar_filename,
                            fPOS=poscar_filename,
                            Ecut=expansion_cutoff,
                            kplist=np.array([14]),
                            spinor=False,
                            searchUC=True)
    characters = bandstr.write_characters(degen_thresh=degenthresh,
                                          kpnames=special_point_chars)
    irrep = {}
    for spc, kpt in zip(special_point_chars, characters["k-points"]):
        chars = []
        for e, i, d in zip(kpt["energies"], kpt["irreps"], kpt["dimensions"]):
            chars.append(Character(find_irrep(i), e, d))
        irrep[spc] = Irrep(tuple(bandstr.kpoints[0].K), chars)

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
