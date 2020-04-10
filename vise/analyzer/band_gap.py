# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from typing import Tuple, Optional, List, Dict

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun, Outcar


def band_gap_from_vasp(vasprun: Vasprun,
                       outcar: Outcar,
                       ) -> Tuple[dict, Optional[dict], Optional[dict]]:
    """Evaluate the band gap properties from vasprun.xml

    Args:
        vasprun (Vasprun):
            Pymatgen Vasprun class instance with eigenvalue info.
        outcar:
    """
    eigenvalues = {spin: e[:, :, 1] for spin, e in vasprun.eigenvalues.items()}
    return band_gap_properties(eigenvalues,
                               outcar.nelect,
                               outcar.magnetization,
                               vasprun.actual_kpoints)


def band_gap_properties(eigenvalues: Dict[Spin, np.ndarray],
                        nelect: float,
                        magnetization: float,
                        kpoints: List[List[float]],
                        frac_threshold: float = 0.1
                        ) -> Tuple[dict, Optional[dict], Optional[dict]]:
    """Evaluate the band gap properties from vasprun.xml and OUTCAR files

    Args:
        eigenvalues:
            Eigenvalues[Spin][band_idx][k-point_idx] = eigenvalue
        nelect:
            Number of electrons.
        magnetization:
            Magnetization in Bohr magneton.
        kpoints (list):
            List of actual k-point fractional coordinates.
        frac_threshold (float):
            Threshold to judge if the number is integer or not, must be between
            0 and 0.5.

    Return:
        Tuple of band_gap, vbm, and cbm-related info. For metals,
        {'energy': 0.0, 'direct': None, 'transition': None}, None, None
        is returned.
    """
    if not 0 < frac_threshold < 0.5:
        raise ValueError("frac_threshold must be between 0 and 0.5.")

    metal_info = {'energy': 0.0, 'direct': None, 'transition': None}, None, None
    if max(nelect % 1, magnetization % 1) > frac_threshold:
        return metal_info

    vbm_info = {'energy': -float("inf")}
    cbm_info = {'energy': float("inf")}

    for spin, eigs_per_spin in eigenvalues.items():

        spin_int = None if len(eigenvalues) == 1 else int(spin)
        hob_band_index = int(round(nelect + magnetization * int(spin)) / 2) - 1
        hob_eigenvalue = np.amax(eigs_per_spin[:, hob_band_index])
        lub_eigenvalue = np.amin(eigs_per_spin[:, hob_band_index + 1])

        if hob_eigenvalue > vbm_info["energy"]:
            hob_k_index = np.where(
                eigs_per_spin[:, hob_band_index] == hob_eigenvalue)[0][0]
            vbm_info = {'energy': round(hob_eigenvalue, 4),
                        'spin': spin_int,
                        'band_index': hob_band_index,
                        'kpoints': kpoints[hob_k_index]}

        if lub_eigenvalue < cbm_info["energy"]:
            lub_k_index = np.where(
                eigs_per_spin[:, hob_band_index + 1] == lub_eigenvalue)[0][0]
            cbm_info = {'energy': round(lub_eigenvalue, 4),
                        'spin': spin_int,
                        'band_index': hob_band_index + 1,
                        'kpoints': kpoints[lub_k_index]}

    if vbm_info["energy"] > cbm_info["energy"]:
        return metal_info

    is_direct = (vbm_info["spin"] == cbm_info["spin"]
                 and vbm_info["kpoints"] == cbm_info["kpoints"])

    band_gap = {"energy": round(cbm_info["energy"] - vbm_info["energy"], 4),
                "direct": is_direct}

    return band_gap, vbm_info, cbm_info


