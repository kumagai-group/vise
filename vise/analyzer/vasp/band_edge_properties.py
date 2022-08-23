# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections import defaultdict

import numpy
from pymatgen.io.vasp import Vasprun, Outcar, Procar
from vise.analyzer.band_edge_properties import BandEdgeProperties


class VaspBandEdgeProperties(BandEdgeProperties):
    def __init__(self,
                 vasprun: Vasprun,
                 outcar: Outcar,
                 integer_criterion: float = 0.1):

        total_mag = outcar.total_mag if outcar.total_mag else 0.0

        self.structure = vasprun.final_structure

        super().__init__(eigenvalues=eigenvalues_from_vasprun(vasprun),
                         nelect=outcar.nelect,
                         magnetization=total_mag,
                         kpoint_coords=vasprun.actual_kpoints,
                         integer_criterion=integer_criterion)


def edge_orbital_contributions(procar: Procar, structure, vbm_info, cbm_info):
    elements = structure.composition.elements
    ele_indices = []
    for e in elements:
        ele_indices.append(
            [i for i, s in enumerate(structure) if s.specie == e])

    # Kpoint index first and then band index
    vbm = procar.data[vbm_info.spin][vbm_info.kpoint_index][vbm_info.band_index]
    cbm = procar.data[cbm_info.spin][cbm_info.kpoint_index][cbm_info.band_index]

    # orbital type
    if len(vbm[0]) == 3:
        orb_indices = [[0], [1], [2]]
    elif len(vbm[0]) == 4:
        orb_indices = [[0], [1], [2], [3]]
    elif len(vbm[0]) == 9:
        orb_indices = [[0], [1, 2, 3], [4, 5, 6, 7, 8]]
    elif len(vbm[0]) == 16:
        orb_indices = [[0], [1, 2, 3], [4, 5, 6, 7, 8],
                       [9, 10, 11, 12, 13, 14, 15]]
    else:
        raise ValueError("")

    vbm_comp, cbm_comp = defaultdict(dict), defaultdict(dict)

    for e, e_idxs in zip(elements, ele_indices):
        for o, o_idxs in zip(["s", "p", "d", "f"], orb_indices):
            vbm_comp[str(e)][o] = vbm[numpy.ix_(e_idxs, o_idxs)].sum()
            cbm_comp[str(e)][o] = cbm[numpy.ix_(e_idxs, o_idxs)].sum()

    return dict(vbm_comp), dict(cbm_comp)


def eigenvalues_from_vasprun(vasprun):
    return {spin: e[:, :, 0] for spin, e in vasprun.eigenvalues.items()}