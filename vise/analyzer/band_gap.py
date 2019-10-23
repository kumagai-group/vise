#!/usr/bin/env python

from typing import Tuple, Optional

from pymatgen.io.vasp.outputs import Vasprun


def band_gap_properties(vasprun: Vasprun,
                        digit: int = 1) -> Optional[Tuple[str, dict, dict]]:
    """
    Args:
        vasprun (Vasprun):
        digit (int):
    """

    for s in vasprun.eigenvalues:
        occupation = [round(sum([i[1] for i in k]), digit)
                      for k in vasprun.eigenvalues[s]]

        if len(set(occupation)) != 1:
            return None

    band_structure = vasprun.get_band_structure()
    vbm = band_structure.get_vbm()
    cbm = band_structure.get_cbm()

    band_gap = band_structure.get_band_gap()

    if band_gap["energy"] == 0.0:
        return None

    vbm_info = {'energy': vbm['energy'],
                'band_index': dict(vbm['band_index']),
                'kpoints':
                    [vasprun.actual_kpoints[x] for x in vbm["kpoint_index"]]}
    cbm_info = {'energy': cbm['energy'],
                'band_index': dict(cbm['band_index']),
                'kpoints':
                    [vasprun.actual_kpoints[x] for x in cbm["kpoint_index"]]}

    return band_gap, vbm_info, cbm_info

