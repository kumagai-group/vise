# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.io.vasp import Vasprun, Outcar

from vise.analyzer.band_edge_properties import BandEdgeProperties


class VaspBandEdgeProperties(BandEdgeProperties):
    def __init__(self,
                 vasprun: Vasprun,
                 outcar: Outcar,
                 integer_criterion: float = 0.1):

        super().__init__(eigenvalues=eigenvalues_from_vasprun(vasprun),
                         nelect=outcar.nelect,
                         magnetization=outcar.total_mag,
                         kpoints=vasprun.actual_kpoints,
                         integer_criterion=integer_criterion)


def eigenvalues_from_vasprun(vasprun):
    return {spin: e[:, :, 0] for spin, e in vasprun.eigenvalues.items()}