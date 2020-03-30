# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from typing import Tuple, Optional

from pymatgen.io.vasp.outputs import Vasprun


def band_gap_properties(vasprun: Vasprun,
                        frac_threshold: float = 0.1
                        ) -> Tuple[dict, Optional[dict], Optional[dict]]:
    """Evaluate the band gap properties from vasprun.xml and OUTCAR files

    Args:
        vasprun (Vasprun):
            Pymatgen Vasprun class instance with eigenvalue info.
        frac_threshold:
            Threshold to judge if the number is integer or not, which must be
            less than 0.5.

    Return:
        Tuple of band_gap, vbm, and cbm-related info. For metals,
        {'energy': 0.0, 'direct': None, 'transition': None}, None, None
        is returned.
    """
    kpoints = vasprun.actual_kpoints
    metal_info = {'energy': 0.0, 'direct': None, 'transition': None}, None, None

    vbm_info = {'energy': -float("inf")}
    cbm_info = {'energy': float("inf")}

    for s, eigval in vasprun.eigenvalues.items():

        occupation = np.sum(eigval[:, :, 1], axis=1)
        if occupation.max() - occupation.min() > frac_threshold:
            return metal_info

        hob_band_index = int(round(occupation.mean())) - 1

        hob_eigenvalue = np.amax(eigval[:, hob_band_index, 0])
        hob_k_index = np.where(
            eigval[:, hob_band_index, 0] == hob_eigenvalue)[0][0]

        lub_eigenvalue = np.amin(eigval[:, hob_band_index + 1, 0])
        lub_k_index = np.where(
            eigval[:, hob_band_index + 1, 0] == lub_eigenvalue)[0][0]

        spin = None if len(vasprun.eigenvalues) == 1 else int(s)

        if hob_eigenvalue > vbm_info["energy"]:
            vbm_info = {'energy': hob_eigenvalue,
                        'spin': spin,
                        'band_index': hob_band_index,
                        'kpoints': kpoints[hob_k_index]}

        if lub_eigenvalue < cbm_info["energy"]:
            cbm_info = {'energy': lub_eigenvalue,
                        'spin': spin,
                        'band_index': hob_band_index + 1,
                        'kpoints': kpoints[lub_k_index]}

    if vbm_info["energy"] > cbm_info["energy"]:
        return metal_info

    is_direct = (vbm_info["spin"] == cbm_info["spin"]
                 and vbm_info["kpoints"] == cbm_info["kpoints"])

    band_gap = {"energy": round(cbm_info["energy"] - vbm_info["energy"], 4),
                "direct": is_direct}

    return band_gap, vbm_info, cbm_info


