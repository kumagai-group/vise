# -*- coding: utf-8 -*-

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun, Outcar, Potcar
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.defaults import defaults


@dataclass
class PriorInfo(MSONable):
    """Used to control the input parameters"""
    structure: Structure = None
    energy_per_atom: float = None
    band_gap: float = None
    vbm_cbm: list = field(default_factory=list)
    total_magnetization: float = None
    data_source: str = None
    is_cluster: bool = None
    charge: int = None
    icsd_ids: List[int] = None
    incar: dict = field(default_factory=dict)

    def dump_yaml(self, filename: str = "prior_info.yaml") -> None:
        with open(filename, "w") as f:
            f.write(yaml.dump(self.as_dict()))

    @classmethod
    def load_yaml(cls, filename: str = "prior_info.yaml"):
        with open(filename, "r") as f:
            d = yaml.load(f, Loader=yaml.SafeLoader)

        return cls.from_dict(d)

    def dump_json(self, filename: str = "prior_info.json") -> None:
        with open(filename, "w") as fw:
            json.dump(self.as_dict(), fw, indent=2)

    @classmethod
    def load_json(cls, filename: str = "prior_info.json"):
        return loadfn(filename)

    @property
    def is_magnetic(self) -> Optional[bool]:
        try:
            return self.total_magnetization > defaults.integer_criterion
        except TypeError:
            return

    @property
    def has_band_gap(self) -> bool:
        return self.band_gap > defaults.band_gap_criterion

    @property
    def is_metal(self) -> bool:
        return not self.has_band_gap

    @property
    def input_options_kwargs(self):
        result = {}
        if self.vbm_cbm:
            result["vbm_cbm"] = self.vbm_cbm
        if isinstance(self.is_magnetic, bool):
            result["is_magnetization"] = self.is_magnetic
        if self.band_gap:
            result["band_gap"] = self.band_gap
        if self.charge:
            result["charge"] = self.charge
        return result


def prior_info_from_calc_dir(prev_dir_path: Path,
                             vasprun: str = "vasprun.xml",
                             outcar: str = "OUTCAR",
                             potcar: str = "POTCAR"):

    vasprun = Vasprun(str(prev_dir_path / vasprun))
    outcar = Outcar(str(prev_dir_path / outcar))
    potcar = Potcar.from_file(str(prev_dir_path / potcar))

    charge = get_net_charge_from_vasp(vasprun.final_structure,
                                      vasprun.parameters["NELECT"],
                                      potcar)
    structure = vasprun.final_structure.copy()
    energy_per_atom = outcar.final_energy / len(structure)
    band_edge_property = VaspBandEdgeProperties(vasprun, outcar)
    total_magnetization = outcar.total_mag

    return PriorInfo(structure=structure,
                     charge=charge,
                     energy_per_atom=energy_per_atom,
                     band_gap=band_edge_property.band_gap,
                     vbm_cbm=band_edge_property.vbm_cbm,
                     total_magnetization=total_magnetization)


def get_net_charge_from_vasp(structure: Structure,
                             nelect: int,
                             potcar: Potcar):
    """
    Returns the defect charge by comparing nion, number of electrons in POTCAR,
    and NELECT in INCAR.
    """
    nuclei_charge = 0

    for elem, potcar in zip(structure.composition, potcar):
        if potcar.element != str(elem):
            raise ValueError("The sequence of elements in POTCAR and Structure "
                             "is different.")
        nuclei_charge += potcar.nelectrons * structure.composition[elem]

    # charge is minus of difference of the electrons
    return int(nuclei_charge - nelect)
