# -*- coding: utf-8 -*-

import json
import shutil
from pathlib import Path
from typing import List

from pymatgen import Element, MPRester, Composition
from vise.chempotdiag.free_energy_entries import logger
from vise.chempotdiag.gas import Gas


def get_mp_materials(elements: List[str],
                     properties: List[str],
                     e_above_hull: float = 1e-4,
                     api_key=None):
    """
    """
    exclude_z = list(i for i in range(1, 100))
    excluded_elements = [str(Element.from_Z(z)) for z in exclude_z]
    for e in elements:
        excluded_elements.remove(e)
    with MPRester(api_key) as m:
        materials = \
            m.query(criteria={"elements": {"$in": elements,
                                           "$nin": excluded_elements},
                              "e_above_hull": {"$lte": e_above_hull}},
                    properties=properties)

    return materials


def make_poscars_from_mp(elements,
                         path: Path = Path.cwd(),
                         e_above_hull=0.01,
                         api_key=None,
                         add_molecules=True,
                         only_molecules=True) -> None:
    """

    Args:
        elements(list): like ["Cu", "O"]
        path(str):
        e_above_hull(float):
        api_key(str):
        add_molecules(bool):
        only_molecules:

    Returns:
        None
    """
    if not path.is_dir:
        raise NotADirectoryError(f"{path} is not directory.")

    mol_dir = Path(__file__).parent / "molecules"

    molecules_formula_list = []
    if add_molecules:
        for g in Gas:
            comp = Composition(str(g))
            if set([str(e) for e in comp.elements]) < set(elements):
                molecules_formula_list.append(comp.reduced_formula)
                dirname = path / f"mol_{str(comp)}"
                if dirname.exists():
                    logger.critical(f"{dirname} exists! So, skip creating it.")
                else:
                    dirname.mkdir()
                    shutil.copyfile(mol_dir / str(comp) / "POSCAR",
                                    dirname / "POSCAR")

    properties = ["task_id",
                  "full_formula",
                  "final_energy",
                  "structure",
                  "spacegroup",
                  "band_gap",
                  "total_magnetization",
                  "magnetic_type"]
    materials = get_mp_materials(elements, properties, e_above_hull, api_key)

    for m in materials:
        comp = Composition(m.pop("full_formula")).reduced_formula
        if only_molecules and comp in molecules_formula_list:
            continue
        m_path = path / f"{m['task_id']}_{comp}"
        m_path.mkdir()
        m.pop("structure").to(filename=m_path / "POSCAR")
        json_path = m_path / "prior_info.json"
        with open(str(json_path), "w") as fw:
            json.dump(m, fw)