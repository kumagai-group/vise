# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from math import sqrt, pi
from typing import List

import numpy as np
from monty.json import MSONable
from tqdm import tqdm
from vise.util.mix_in import ToJsonFileMixIn
from scipy.constants import physical_constants as pc

eV_to_inv_cm = pc["electron volt-inverse meter relationship"][0] / 100


def diele_func_to_coeff(freq, real, imag):
    return (2 * sqrt(2) * pi * sqrt(sqrt(real ** 2 + imag ** 2) - real)
            * freq * eV_to_inv_cm)


@dataclass
class DieleFuncData(MSONable, ToJsonFileMixIn):
    energies: List[float]  # in eV
    diele_func_real: List[List[float]]  # [xx, yy, zz, xy, yz, xz]
    diele_func_imag: List[List[float]]  # [xx, yy, zz, xy, yz, xz]
    band_gap: float  # in eV

    @property
    def ave_absorption_coeff(self):
        reals = [sum(self.diele_func_real[i][:3]) / 3
                 for i in range(len(self.energies))]
        imags = [sum(self.diele_func_imag[i][:3]) / 3
                 for i in range(len(self.energies))]
        return [diele_func_to_coeff(freq, real, imag)
                for freq, real, imag in zip(self.energies, reals, imags)]

    def target_coeff_min_e(self, target_coeff: float = 10**4):
        for e, coeff in zip(self.energies, self.ave_absorption_coeff):
            if coeff > target_coeff:
                return e
        return None


def make_shifted_diele_func(diele_func_data: DieleFuncData,
                            original_band_gap: float,
                            shift: float) -> DieleFuncData:
    imag = imag_shift(diele_func_data.diele_func_imag,
                      diele_func_data.energies,
                      original_band_gap + shift, shift)
    real = kramers_kronig_trans(imag, diele_func_data.energies)
    return DieleFuncData(diele_func_data.energies,
                         real.tolist(),
                         imag.tolist(),
                         original_band_gap + shift)


def imag_shift(diele_func_imag: List[List[float]],
               energies: List[float],
               band_gap: float,
               shift: float) -> np.ndarray:
    energies = np.array(energies)
    assert shift > 0
    result = []
    for energy_grid in energies:
        old_e = energy_grid - shift
        right_idx = np.argwhere(energies > old_e)[0][0]
        left_e, right_e = energies[right_idx - 1], energies[right_idx]
        # linear interpolation
        left_ratio = (right_e - old_e) / (right_e - left_e)

        inner_result = []
        for imag_idx in range(6):
            if energy_grid < band_gap:
                inner_result.append(0.0)
            else:
                old_diele = \
                    diele_func_imag[right_idx - 1][imag_idx] * left_ratio + \
                    diele_func_imag[right_idx][imag_idx] * (1 - left_ratio)
                inner_result.append(
                    old_diele * (energy_grid - shift) / energy_grid)

        result.append(inner_result)

    return np.array(result)


def kramers_kronig_trans(diele_func_imag: np.array,
                         energies: List[float],
                         ita: float = 0.01) -> np.ndarray:
    mesh = energies[1] - energies[0]
    result = []
    ee2ss = [[e ** 2 - energy_grid ** 2 for e in energies]
             for energy_grid in energies]
    for imag_idx in tqdm(range(6)):
        imags = diele_func_imag[:, imag_idx]
        if imag_idx == 0 or \
                (imag_idx > 0
                 and np.allclose(
                            imags, diele_func_imag[:, imag_idx - 1]) is False):
            if np.count_nonzero(imags) == 0:
                inner_result = [0.0] * len(energies)
            else:
                inner_result = []
                for ee2s in ee2ss:
                    integrals = [e * imag * ee2 / (ee2 ** 2 + ita ** 2)
                                 for e, ee2, imag in zip(energies, ee2s, imags)]
                    integral = sum(integrals) * mesh * 2 / pi
                    if imag_idx < 3:
                        integral += 1
                    inner_result.append(integral)

        result.append(inner_result)

    return np.array(result).T