# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List
import numpy as np

from monty.json import MSONable
from numpy import linalg
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class EffectiveMass(MSONable, ToJsonFileMixIn):
    p: List[List[List]]  # [temperature][carrier concentration]
    n: List[List[List]]
    concentrations: List[float]

    def effective_mass(self, carrier_type, temp, concentration):
        i_c = self.concentrations.index(concentration)
        return self.__getattribute__(carrier_type)[i_c]


def eigvals_and_vecs(matrix: np.ndarray):
    assert matrix.shape == (3, 3)
    eigvals, eigvecs = linalg.eig(matrix)
    idx = eigvals.argsort()
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx].T  # Need to transverse
    return eigvals, eigvecs


def lowest_eigval_and_vecs(matrix: np.ndarray):
    eigvals, eigvecs = eigvals_and_vecs(matrix)
    assert eigvals[0] > 0.0
    eigval = eigvals[0]
    lowest_eigvals = np.where(eigvals < eigval * 1.0001)
    eigvecs = np.round(eigvecs[lowest_eigvals], 2)
    return eigval, eigvecs
