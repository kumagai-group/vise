#!/usr/bin/env python

from typing import Tuple, Optional
import math
from pymatgen.io.vasp.outputs import Vasprun


def band_gap_properties(vasprun: Vasprun,
                        digit: int = 1) -> Optional[Tuple]:
    """
    Args:
        vasprun (Vasprun):
        digit (int):
    """
    weight = vasprun.actual_kpoints_weights
    for v in vasprun.eigenvalues.values():
        for i in range(len(v[0])):
            data_along_k = v[:, i]
            # array([-12.8753,   1.    ]) Average of (energy, occupation)
            average = sum([data_along_k[1] * w for v, w in zip(v, weight)])
            frac_occu = round(average[1], digit) - average[1]
            print(frac_occu)
            if not math.isclose(frac_occu, 0, rel_tol=0.01):
                band_gap = {'energy': 0.0, 'direct': False, 'transition': None}
                return band_gap, None, None

    band_structure = vasprun.get_band_structure()
    vbm = band_structure.get_vbm()
    cbm = band_structure.get_cbm()

    band_gap = band_structure.get_band_gap()

    if band_gap["energy"] == 0.0:
        return band_gap, None, None

    vbm_band = {str(k): v for k, v in vbm['band_index'].items()}
    cbm_band = {str(k): v for k, v in cbm['band_index'].items()}

    vbm_info = {'energy': vbm['energy'],
                'band_index': vbm_band,
                'kpoints':
                    [vasprun.actual_kpoints[x] for x in vbm["kpoint_index"]]}
    cbm_info = {'energy': cbm['energy'],
                'band_index': cbm_band,
                'kpoints':
                    [vasprun.actual_kpoints[x] for x in cbm["kpoint_index"]]}

    return band_gap, vbm_info, cbm_info

