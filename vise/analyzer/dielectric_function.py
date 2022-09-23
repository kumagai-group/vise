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


def absorption_coeff(energy: float, real: float, imag: float) -> float:
    return (2 * sqrt(2) * pi * sqrt(sqrt(real ** 2 + imag ** 2) - real)
            * energy * eV_to_inv_cm)


def refractive_idx_real(e_real: float, e_imag: float) -> float:
    return sqrt(e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)


def refractive_idx_imag(e_real: float, e_imag: float) -> float:
    return sqrt(-e_real + sqrt(e_real ** 2 + e_imag ** 2)) / sqrt(2)


@dataclass
class DieleFuncData(MSONable, ToJsonFileMixIn, ToCsvFileMixIn):
    energies: List[float]  # in eV
    directions: List[str]  # ["xx", "yy", "ave"]
    diele_func_real: List[List[float]]  # [list_of_xx, list_of_xy, list_of_ave]
    diele_func_imag: List[List[float]]  # [list_of_xx, list_of_xy, list_of_ave]
    band_gap: Optional[float] = None  # in eV

    def __post_init__(self):
        try:
            assert len(self.directions) == len(self.diele_func_real)
        except AssertionError:
            print(f"{len(self.directions)} vs {len(self.diele_func_real)}")
            raise

    @property
    def real_columns(self):
        return [f"real_{d}" for d in self.directions]

    @property
    def imag_columns(self):
        return [f"imag_{d}" for d in self.directions]

    @property
    def to_dataframe(self) -> pd.DataFrame:
        d = {"energies(eV)": self.energies}
        for x, y in zip(self.real_columns, self.diele_func_real):
            d[x] = y
        for x, y in zip(self.imag_columns, self.diele_func_imag):
            d[x] = y
        d["band_gap"] = [None] * len(self.energies)
        d["band_gap"][0] = self.band_gap

        return pd.DataFrame.from_dict(d)

    @classmethod
    def from_dataframe(cls, df):
        real_T, imag_T, directions = [], [], []
        for column_name, item in df.iteritems():
            if column_name in ["energies(eV)", "band_gap"]:
                continue
            elif "real" in column_name:
                real_T.append(item.tolist())
                directions.append(column_name.split("_")[-1])
            elif "imag" in column_name:
                imag_T.append(item.tolist())
            else:
                raise KeyError("The input CSV does not have proper format.")

        return cls(energies=df["energies(eV)"].tolist(),
                   diele_func_real=real_T,
                   diele_func_imag=imag_T,
                   directions=directions,
                   band_gap=float(df.loc[0, "band_gap"]))

    @property
    def absorption_coeff(self):
        return [[absorption_coeff(e, r, i)
                for e, r, i in zip(self.energies, reals, imags)]
                for reals, imags in
                zip(self.diele_func_real, self.diele_func_imag)]

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
                         diele_func_data.directions,
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
        for imag_idx in range(len(diele_func_imag)):
            if energy_grid < band_gap:
                inner_result.append(0.0)
            else:
                old_diele = \
                    diele_func_imag[imag_idx][right_idx - 1] * left_ratio + \
                    diele_func_imag[imag_idx][right_idx] * (1 - left_ratio)
                inner_result.append(
                    old_diele * (energy_grid - shift) / energy_grid)

        result.append(inner_result)

    return np.array(result).T


def kramers_kronig_trans(diele_func_imag: np.array,
                         energies: List[float],
                         ita: float = 0.01) -> np.ndarray:
    mesh = energies[1] - energies[0]
    result = []
    ee2ss = [[e ** 2 - energy_grid ** 2 for e in energies]
             for energy_grid in energies]
    for imag_idx in tqdm(range(len(diele_func_imag))):
        imags = diele_func_imag[imag_idx, :]
        if imag_idx == 0 or \
                (imag_idx > 0
                 and np.allclose(
                            imags, diele_func_imag[imag_idx - 1, :]) is False):
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

    return np.array(result)