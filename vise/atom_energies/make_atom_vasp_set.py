# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import warnings
from pathlib import Path
from typing import Union, List

import fire
import yaml
from monty.serialization import loadfn
from pymatgen.core import Structure, Lattice, Element
from pymatgen.io.vasp import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet, BadInputSetWarning
from vise.input_set.datasets.potcar_set import PotcarSet

mags = loadfn(Path(__file__).parent / "mp_atom_mag.yaml")


warnings.simplefilter('ignore', BadInputSetWarning)


def is_target_element(elem: Union[Element, str]):
    elem = Element(elem) if type(elem) is str else elem
    return elem.Z <= 58 or 71 <= elem.Z <= 83


def make_atom_poscar_dirs(path: Path, elems: List[Element] = None):
    elems = [str(e) for e in elems] if elems else mags.keys()

    for element in elems:
        d = path / element
        Path(d).mkdir()
        structure = Structure(Lattice.cubic(10), [element], [[0.5]*3])
        structure.to(fmt="POSCAR", filename=str(d / "POSCAR"))
        nupdown = mags[element]
        prior_info = {"incar": {"ISPIN": 2, "NUPDOWN": float(nupdown),
                                "NELM": 300},
                      "is_cluster": True}
        (d / "prior_info.yaml").write_text(yaml.dump(prior_info))


def make_atom_mp_relax_set():
    for element, potcar in PotcarSet.mp_relax_set.potcar_dict().items():
        if potcar is None:
            continue

        if is_target_element(element) is False:
            continue

        structure = Structure(Lattice.cubic(10),
                              coords=[[0.5]*3], species=[element])

        mp_set = MPRelaxSet(structure,
                            user_kpoints_settings=Kpoints(kpts=((1, 1, 1),)),
                            user_incar_settings={"ALGO": "D",
                                                 "ISIF": 2,
                                                 "ISMEAR": 0,
                                                 "MAGMOM": {"H": 1.0},
                                                 "NELM": 300})
        Path(element).mkdir()
        mp_set.write_input(element)


if __name__ == '__main__':
    def make_mp_set():
        make_atom_mp_relax_set()

    fire.Fire(make_mp_set)
