# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.

from dataclasses import dataclass
from typing import List, Dict

from monty.json import MSONable


@dataclass
class BZPlotInfo(MSONable):
    faces: List[List[List[float]]]
    # {"X": {"cart": [0.5, 0, 0], "frac": [0.7514, 0, 0]"}}
    labels: Dict[str, Dict[str, List[float]]]
    band_paths: List[List[List[float]]] = None
    rec_lat_vec: List[List[float]] = None
