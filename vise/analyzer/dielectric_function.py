# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from math import sqrt, pi
from typing import List, Optional

import numpy as np
import pandas as pd
from monty.json import MSONable
from tqdm import tqdm
from vise.util.mix_in import ToJsonFileMixIn, ToCsvFileMixIn
from vise.util.logger import get_logger
from scipy.constants import physical_constants as pc

logger = get_logger(__name__)

eV_to_inv_cm = pc["electron volt-inverse meter relationship"][0] / 100


def diele_func_to_coeff(energy: float, real: float, imag: float) -> float:
    return (2 * sqrt(2) * pi * sqrt(sqrt(real ** 2 + imag ** 2) - real)
            * energy * eV_to_inv_cm)


def refractive_idx_real(e_real: float, e_imag: float) -> float:
    return sqrt(e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)


def refractive_idx_imag(e_real: float, e_imag: float) -> float:
    return sqrt(-e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)


@dataclass
class DieleFuncData(MSONable, ToJsonFileMixIn, ToCsvFileMixIn):
    energies: List[float]  # in eV
    diele_func_real: List[List[float]]  # [xx, yy, zz, xy, yz, xz]
    diele_func_imag: List[List[float]]  # [xx, yy, zz, xy, yz, xz]
    band_gap: Optional[float] = None  # in eV

    @classmethod
    def real_columns(cls):
        return [f"real_{d}" for d in ["xx", "yy", "zz", "xy", "yz", "xz"]]

    @classmethod
    def imag_columns(cls):
        return [f"imag_{d}" for d in ["xx", "yy", "zz", "xy", "yz", "xz"]]

    @property
    def csv_column_names(self):
        result = ["energies(eV)"]
        result.extend(self.real_columns())
        result.extend(self.imag_columns())
        return result

    @property
    def csv_data(self):
        result = []
        for i, j, k in zip(self.energies, self.diele_func_real, self.diele_func_imag):
            result.append([i] + j + k)
        return result

    @classmethod
    def from_csv(cls, filename: str):
        df = pd.read_csv(filename)
        real = df.loc[:, cls.real_columns()].values
        imag = df.loc[:, cls.imag_columns()].values
        return cls(energies=df["energies(eV)"].tolist(),
                   diele_func_real=real, diele_func_imag=imag)

    @property
    def absorption_coeff(self):
        return [[diele_func_to_coeff(e, r, i) for r, i in zip(reals, imags)]
                for e, reals, imags in
                zip(self.energies, self.diele_func_real, self.diele_func_imag)]

    @property
    def refractive_idx_real(self):
        return [[refractive_idx_real(r, i) for r, i in zip(reals, imags)]
                for reals, imags in
                zip(self.diele_func_real, self.diele_func_imag)]

    @property
    def refractive_idx_imag(self):
        return [[refractive_idx_imag(r, i) for r, i in zip(reals, imags)]
                for reals, imags in
                zip(self.diele_func_real, self.diele_func_imag)]

    @property
    def reflectivity(self):
        return [[((n - 1) ** 2 + k ** 2) / ((n + 1) ** 2 + k ** 2)
                 for n, k in zip(ns, ks)]
                for ns, ks in zip(self.refractive_idx_real,
                                  self.refractive_idx_imag)]


def min_e_w_target_coeff(energies, quantities, target_value):
    for e, quantity in zip(energies, quantities):
        if quantity > target_value:
            return e
    logger.warning(
        f"The target value {target_value} is not reached in the entire range.")
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