# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from typing import List, Tuple

import numpy as np
from irreptables import IrrepTable
from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints

from vise.analyzer.plot_band import Irrep, Irreps
from vise.analyzer.vasp.plot_band import greek_to_unicode
from vise.error import ViseError
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

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


def make_irreps_from_wavecar(special_point_symbols: List[str],
                             kpt_indices: List[int],
                             sg_num: int = None,
                             wavecar_filename: str = "WAVECAR",
                             poscar_filename: str = "POSCAR",
                             plane_wave_cutoff: float = 50.0,
                             degeneracy_threshold: float = 0.01,
                             **bs_kwargs) -> Irreps:

    if sg_num is None:
        structure = Structure.from_file(poscar_filename)
        sg_analyzer = StructureSymmetrizer(structure=structure)
        sg_num = sg_analyzer.sg_number

    logger.info("We set spinor is False.")
    irrep_table = IrrepTable(sg_num, spinor=False)
    symbols = {i.kpname for i in irrep_table.irreps}

    listed_k_indices, listed_symbols = [], []
    for s, k_index in zip(special_point_symbols, kpt_indices):
        if s in symbols:
            listed_symbols.append(s)
            listed_k_indices.append(k_index)
        else:
            logger.info(f"{s} is not listed in the IrrepTable.")

    bs = BandStructure(fWAV=wavecar_filename,
                       fPOS=poscar_filename,
                       Ecut=plane_wave_cutoff,
                       kplist=np.array(listed_k_indices),
                       spinor=False,
                       **bs_kwargs)

    characters = bs.write_characters(degen_thresh=degeneracy_threshold,
                                     kpnames=listed_symbols)

    irrep_dict = {}
    for symbol, c_kpt, kpt in \
            zip(listed_symbols, characters["k-points"], bs.kpoints):
        symbols = []
        for irrep in c_kpt["irreps"]:
            try:
                irrep_str = find_irrep(irrep)
            except ViseNoIrrepError:
                irrep_str = "Unknown"
            symbols.append(greek_to_unicode(irrep_str))

        irrep_dict[greek_to_unicode(symbol)] = \
            Irrep(kpt.K.tolist(),
                  symbols,
                  c_kpt["energies"].tolist(),
                  c_kpt["dimensions"].tolist())

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
