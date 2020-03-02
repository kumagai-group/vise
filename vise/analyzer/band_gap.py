# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from typing import Tuple, Optional, Union
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.electronic_structure.core import Spin


def band_gap_properties(vasprun: Union[Vasprun, str],
                        outcar: Union[Outcar, str],
                        int_m_threshold: float = 0.1
                        ) -> Tuple[dict, Optional[dict], Optional[dict]]:
    """Evaluate the band gap properties from vasprun.xml and OUTCAR files

    Args:
        vasprun (Vasprun/ str):
            Vasprun file.
        outcar (Outcar/ str):
           OUTCAR file.
        int_m_threshold:
            The threshold to judge if the magnetization is not integer.

    Return:
        Tuple of band_gap, vbm, and cbm-related info. For metals,
        {'energy': 0.0, 'direct': False, 'transition': None}, None, None
        is returned.
    """
    if isinstance(vasprun, str):
        vasprun = Vasprun(vasprun)
    if isinstance(outcar, str):
        outcar = Outcar(outcar)

    eigenvalues = vasprun.eigenvalues
    kpts = vasprun.actual_kpoints
    metal = {'energy': 0.0, 'direct': False, 'transition': None}, None, None
    mag = outcar.total_mag

    if mag is None:
        if abs(outcar.nelect - round(outcar.nelect / 2) * 2) > int_m_threshold:
            return metal
        hob_index = round((outcar.nelect / 2) - 1)

        hob_eig, hob_k_index, lub_eig, lub_k_index = \
            edge_info(eigenvalues, hob_index)

        vbm_info = {'energy': hob_eig,
                    'spin': None,
                    'band_index': hob_index,
                    'kpoints': kpts[hob_k_index]}
        cbm_info = {'energy': lub_eig,
                    'spin': None,
                    'band_index': hob_index + 1,
                    'kpoints': kpts[lub_k_index]}

    else:
        frac_occupation = mag - round(mag)
        if abs(frac_occupation) > int_m_threshold:
            return metal

        up_hob_index = round((outcar.nelect + (mag + 0.1)) / 2) - 1
        up_hob_eig, up_hob_k_index, up_lub_eig, up_lub_k_index = \
            edge_info(eigenvalues, up_hob_index)

        down_hob_index = round((outcar.nelect - (mag + 0.1)) / 2) - 1
        down_hob_eig, down_hob_k_index, down_lub_eig, down_lub_k_index = \
            edge_info(eigenvalues, down_hob_index, Spin.down)

        if up_hob_eig >= down_hob_eig:
            vbm_info = {'energy': up_hob_eig,
                        'spin': int(Spin.up),
                        'band_index': up_hob_index,
                        'kpoints': kpts[up_hob_k_index]}
        else:
            vbm_info = {'energy': down_hob_eig,
                        'spin': int(Spin.down),
                        'band_index': down_hob_index,
                        'kpoints': kpts[down_hob_k_index]}

        if up_lub_eig <= down_lub_eig:
            cbm_info = {'energy': up_lub_eig,
                        'spin': int(Spin.up),
                        'band_index': up_hob_index + 1,
                        'kpoints': kpts[up_lub_k_index]}
        else:
            cbm_info = {'energy': down_lub_eig,
                        'spin': int(Spin.down),
                        'band_index': down_hob_index + 1,
                        'kpoints': kpts[down_lub_k_index]}

    if vbm_info["energy"] > cbm_info["energy"]:
        return metal

    is_direct = (vbm_info["spin"] == cbm_info["spin"]
                 and vbm_info["kpoints"] == cbm_info["kpoints"])

    band_gap = {"energy": cbm_info["energy"] - vbm_info["energy"],
                "direct": is_direct}

    return band_gap, vbm_info, cbm_info


def edge_info(eigenvalues: dict,
              hob_index: int,
              spin=Spin.up
              ) -> Tuple[float, int, float, int]:
    """Evaluate the band edge info from eigenvalues and given band index.

    Args:
        eigenvalues (dict):
            Dict of eigenvalues.
            eigenvalue, occupation =
            Eigenvalues[spin][k-point index][band index][0:2]
        hob_index (int):
            Highest-occupied band (HOB) index starting from0.
        spin:
            Pymatgen Spin object.

    Return:
        Tuple of HOB eigenvalue, its k index, lowest-unoccupied band (LUB)
        eigenvalue, and its k index.
    """

    hob_eig = np.amax(eigenvalues[spin][:, hob_index, 0])
    hob_k_index = int(np.where(
        eigenvalues[spin][:, hob_index, 0] == hob_eig)[0][0])

    lub_eig = np.amin(eigenvalues[spin][:, hob_index + 1, 0])
    lub_k_index = int(np.where(
        eigenvalues[spin][:, hob_index + 1, 0] == lub_eig)[0][0])

    return hob_eig, hob_k_index, lub_eig, lub_k_index

