# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from typing import Tuple, Optional, List, Dict

from dataclasses import dataclass

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun, Outcar


@dataclass()
class BandEdge:
    energy: float
    spin: Spin = None
    band_index: int = None
    kpoint_coords: List[float] = None

    def is_direct(self, other: "BandEdge"):
        return (self.spin == other.spin and
                self.kpoint_coords == other.kpoint_coords)


class BandEdgeProperties:
    def __init__(self,
                 eigenvalues: Dict[Spin, np.ndarray],
                 nelect: float,
                 magnetization: float,
                 kpoints: List[List[float]],
                 integer_criterion: float = 0.1):

        assert 0 < integer_criterion < 0.5

        # [Spin][k-point_idx][band_idx] = eigenvalue
        self._eigenvalues = eigenvalues
        self._nelect = nelect
        #  In Bohr magneton.
        self._magnetization = magnetization
        self._kpoints = kpoints
        self._integer_criterion = integer_criterion

        self._calculate_vbm_cbm()

        if self._is_metal():
            self.vbm = None
            self.cbm = None

    def _calculate_vbm_cbm(self):
        self.vbm = BandEdge(float("-inf"))
        self.cbm = BandEdge(float("inf"))
        for spin, eigenvalues in self._eigenvalues.items():
            # ho = highest occupied
            ho_band_index = self._hob_band_index(spin)
            ho_eigenvalue = np.amax(eigenvalues[:, ho_band_index])
            if ho_eigenvalue > self.vbm.energy:
                self.vbm = self.band_edge(
                    eigenvalues, ho_band_index, ho_eigenvalue, spin)

            # lu = lowest unoccupied
            lu_band_index = ho_band_index + 1
            lu_eigenvalue = np.amin(eigenvalues[:, lu_band_index])
            if lu_eigenvalue < self.cbm.energy:
                self.cbm = self.band_edge(
                    eigenvalues, lu_band_index, lu_eigenvalue, spin)

    def band_edge(self, eigenvalues, band_index, eigenvalue, spin):
        k_index = np.where(eigenvalues[:, band_index] == eigenvalue)[0][0]
        return BandEdge(eigenvalue, spin, band_index, self._kpoints[k_index])

    def _hob_band_index(self, spin):
        if spin == Spin.up:
            num_occupied_band = (self._nelect + self._magnetization) / 2
        else:
            num_occupied_band = (self._nelect - self._magnetization) / 2

        return round(num_occupied_band) - 1

    def _is_metal(self):
        nelect_frac = abs(round(self._nelect) - self._nelect)
        is_nelect_frac = nelect_frac > self._integer_criterion

        mag_frac = abs(round(self._magnetization) - self._magnetization)
        is_mag_frac = mag_frac > self._integer_criterion

        is_vbm_higher_than_cbm = self.vbm.energy > self.cbm.energy

        return is_nelect_frac or is_mag_frac or is_vbm_higher_than_cbm

    @property
    def is_direct(self):
        return self.vbm.is_direct(self.cbm) if self.vbm else None

    @property
    def band_gap(self):
        return self.cbm.energy - self.vbm.energy if self.vbm else None


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

