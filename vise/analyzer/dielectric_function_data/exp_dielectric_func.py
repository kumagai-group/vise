# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path
from pandas import DataFrame
import pandas as pd
from vise.util.logger import get_logger

logger = get_logger(__name__)

parent_dir = Path(__file__).parent


band_gap = {"GaAs": 1.41,
            "Si": 1.17}


class ExpDieleFunc:
    def __init__(self, material):
        self.material = material
        try:
            self.dataframe = pd.read_csv(parent_dir / f"{material}.csv", sep='\s+')
        except FileNotFoundError:
            logger.warning(f"Csv file for {material} does not exist.")
            raise

#        self.dataframe.absorption = self.dataframe.absorption.astype(float)

    @property
    def energies(self):
        return self.dataframe.energy

    @property
    def dielectric_real(self):
        return self.dataframe.real_part

    @property
    def dielectric_imag(self):
        return self.dataframe.imaginary_part

    @property
    def absorption_coeff(self):
        return self.dataframe.absorption

    @property
    def band_gap(self):
        return band_gap[self.material]

    @property
    def reference(self):
        return {"title": "THE HANDBOOK ON OPTICAL CONSTANTS OF SEMICONDUCTORS",
                "author": "Sadao Adachi",
                "publisher": "World Scientific",
                "doi": "doi.org/10.1142/8480",
                "year": 2012}
