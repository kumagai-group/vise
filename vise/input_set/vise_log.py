# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from copy import copy
from dataclasses import dataclass
from typing import Optional, Dict

import yaml
from monty.json import MSONable
from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.mix_in import ToYamlFileMixIn


@dataclass
class ViseLog(MSONable, ToYamlFileMixIn):
    version: str
    task: Task
    xc: Xc
    input_options: dict
    user_incar_settings: Optional[dict] = None
    ldauu: Optional[Dict[str, float]] = None
    ldaul: Optional[Dict[str, float]] = None

    def to_yaml(self):
        return yaml.safe_dump(self.as_dict(), sort_keys=False)

    def as_dict(self) -> dict:
        d = {"version": self.version, "task": str(self.task),
             "xc": str(self.xc), "input_options": self.input_options}
        if self.user_incar_settings is not None:
            d["user_incar_settings"] = self.user_incar_settings
        if self.ldauu is not None:
            d["ldauu"] = self.ldauu
        if self.ldaul is not None:
            d["ldaul"] = self.ldaul
        if "kpt_mode" in self.input_options:
            self.input_options["kpt_mode"] = str(self.input_options["kpt_mode"])
        if "potcar_set" in self.input_options:
            self.input_options["potcar_set"] = \
                str(self.input_options["potcar_set"])
        return d

    @classmethod
    def from_dict(cls, d):
        dd = copy(d)
        dd["task"] = Task.from_string(d["task"])
        dd["xc"] = Xc.from_string(d["xc"])
        if "kpt_mode" in d["input_options"]:
            dd["input_options"]["kpt_mode"] = \
                KpointsMode.from_string(dd["input_options"]["kpt_mode"])
        if "potcar_set" in d["input_options"]:
            dd["input_options"]["potcar_set"] = \
                PotcarSet.from_string(dd["input_options"]["potcar_set"])

        return super().from_dict(dd)
