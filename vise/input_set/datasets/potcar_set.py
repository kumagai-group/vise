# -*- coding: utf-8 -*-

#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
from pathlib import Path
from typing import Optional, Dict

from monty.serialization import loadfn

from vise.util.enum import ExtendedEnum


class PotcarSet(ExtendedEnum):
    normal = "normal"
    mp_relax_set = "mp"
    gw = "gw"

    def overridden_potcar_dict(self,
                               override_potcar_set: Optional[dict] = None
                               ) -> Dict[str, str]:
        result = deepcopy(self.potcar_dict())
        if override_potcar_set:
            result.update(override_potcar_set)

        return result  # e.g. {"Zr": "Zr_pv", ...}

    def potcar_dict(self) -> Dict[str, str]:
        _potcar_list = loadfn(Path(__file__).parent / "potcar_set.yaml")
        set_names = _potcar_list.pop("set_names")
        result = {}
        for name in set_names:
            result[name] = {}  # normal, mp_relax_set, gw

        for element, potcar_string in _potcar_list.items():
            potcars = potcar_string.split()

            def sanitize(val: str):
                if val == "---":
                    return potcars[0]
                elif val == "None":
                    return None
                else:
                    return val

            for index, potcar_single in enumerate(potcars):
                set_name = set_names[index]
                result[set_name][element] = sanitize(potcar_single)

        return result[self.value]