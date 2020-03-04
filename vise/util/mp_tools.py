# -*- coding: utf-8 -*-

import json
import shutil
from pathlib import Path
from typing import List, Optional
import yaml

from pymatgen import Element, MPRester, Composition

from vise.util.logger import get_logger
from vise.chempotdiag.gas import Gas

logger = get_logger(__name__)


def get_mp_materials(elements: List[str],
                     properties: List[str],
                     e_above_hull: float = 1e-4,
                     api_key: str = None) -> list:
    """Get Materials Project materials as list

    Args:
        elements (list):
            Element list like ["Cu", "O"]
        properties:
            Returned properties.
        e_above_hull:
            Energy criterion in eV/atom.
        api_key:
            Materials Project REST api key.

    Returns:
        List of materials composed of properties.
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


def make_poscars_from_mp(elements: list,
                         path: Path = Path.cwd(),
                         e_above_hull: float = 0.01,
                         api_key: str = None,
                         molecules: bool = True,
                         properties: Optional[list] = None,
                         only_unary=False) -> None:
    """

    Args:
        elements(list):
        e_above_hull(float):
        api_key(str):
        properties (list):
            See docstrings of get_mp_materials above.
        path(str):
            Path to create directories.
        molecules(bool):
            Whether to replace molecular entries to molecules in a large box.
        only_unary:
            Create only unary simple substances.

    Returns:
        None
    """
    if not path.is_dir:
        raise NotADirectoryError(f"{path} is not directory.")

    mol_dir = Path(__file__).parent / ".." / "chempotdiag" / "molecules"

    molecules_formula_list = []
    if molecules:
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
                    shutil.copyfile(mol_dir / str(comp) / "prior_info.yaml",
                                    dirname / "prior_info.yaml")

    properties = properties or ["task_id",
                                "full_formula",
                                "final_energy",
                                "structure",
                                "spacegroup",
                                "band_gap",
                                "total_magnetization",
                                "magnetic_type"]
    materials = get_mp_materials(elements, properties, e_above_hull, api_key)

    for m in materials:
        comp = Composition(m.pop("full_formula"))
        formula = comp.reduced_formula
        if molecules and formula in molecules_formula_list:
            continue
        if only_unary and len(comp) > 1:
            continue

        m_path = path / f"{formula}_{m['task_id']}"
        m_path.mkdir()
        m.pop("structure").to(filename=m_path / "POSCAR")
        yaml_path = m_path / "prior_info.yaml"
        with open(str(yaml_path), "w") as fw:
            fw.write(yaml.dump(m))
