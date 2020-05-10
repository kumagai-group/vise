# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional

from pymatgen.io.vasp.sets import Potcar

from vise.error import ViseError
from vise.input_set.datasets.dataset_util import PotcarSet
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)


def generate_potcar(symbol_list: list,
                    xc: Xc,
                    potcar_set: PotcarSet = PotcarSet.normal,
                    overridden_potcar: Optional[dict] = None):

    potcar_dict = potcar_set.overridden_potcar_dict(overridden_potcar)
    try:
        potcar_symbols = [potcar_dict[el] for el in symbol_list]
    except KeyError as e:
        raise ViseNoPotcarError(e)

    return Potcar(potcar_symbols, functional=xc.potcar_functional)


class ViseNoPotcarError(ViseError):
    pass

