# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List
import numpy as np

from monty.json import MSONable
from numpy import linalg
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class EffectiveMass(MSONable, ToJsonFileMixIn):
    p: List[List[List]]
    n: List[List[List]]
    temperature: float
    concentrations: List[float]

    def average_mass(self, carrier_type, concentration):
        em = self.effective_mass(carrier_type, concentration)
        return sum(em[i][i] for i in range(3)) / 3

    def minimum_mass(self, carrier_type, concentration):
        em = self.effective_mass(carrier_type, concentration)
        result = lowest_eigval_and_vecs(np.array(em))[0]
        if isinstance(result, complex):
            if abs(result.imag) < 1e-3:
                result = result.real
            else:
                raise ValueError("The minimum effective mass shows complex "
                                 f"value of f{result}.")
        return result

    def effective_mass(self, carrier_type, concentration):
        i_c = self.concentrations.index(concentration)
        return self.__getattribute__(carrier_type)[i_c]

    def __str__(self):
        lines = [f"temperature: {self.temperature}"]
        for c, pp, nn in zip(self.concentrations, self.p, self.n):
            lines.append("-"*30)
            lines.append(f"concentration: {c:g}")
            lines.append(f"p:")
            lines.append(tabulate(pp))
            lines.append(f"n:")
            lines.append(tabulate(nn))
        return "\n".join(lines)


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
