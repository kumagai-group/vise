# -*- coding: utf-8 -*-
from math import ceil
from pathlib import Path
from typing import Optional
from typing import Tuple

from monty.serialization import loadfn
from pymatgen import Composition
from pymatgen.io.vasp import Potcar
from vise.input_set.datasets.element_parameters import unoccupied_bands
from vise.util.logger import get_logger

logger = get_logger(__name__)
SET_DIR = Path(__file__).parent / "datasets"


def load_potcar_yaml(set_name: Optional[str] = "normal",
                     override_potcar_set: Optional[dict] = None) -> dict:
    """Load the yaml setting files for config and POTCAR list.

    Args:
        set_name (str):
            Potcar set name.
        override_potcar_set (dict):
            User specifying POTCAR set

    Returns:
          Dictionary of potcar_set, like {"Zr": "Zr_pv", ...}
    """
    potcar_set = loadfn(SET_DIR / "potcar_set.yaml")
    try:
        potcar = potcar_set[set_name]
    except KeyError:
        logger.critical(f"The accepted potcar set name is {potcar_set.keys()}")
        raise

    if override_potcar_set:
        potcar.update(override_potcar_set)

    return potcar


def load_default_incar_settings(yaml_filename: str,
                                required_flags: set,
                                optional_flags: set,
                                key_name: str) -> dict:
    """Load the yaml setting files for config and POTCAR list.

    Args:
        yaml_filename (str):
            Yaml filename with xxx_incar_set.yaml.
        required_flags (set):
            INCAR flags that must exist in the dict as a value of key_name key
            in yaml file.
        optional_flags (set):
            INCAR flags that that could exist in the dict as a value of key_name
            key in yaml file.
        key_name (str):
            Key name such as "structure_opt" or "hse".

    Returns:
          settings (dict):
    """
    incar_set = loadfn(SET_DIR / yaml_filename)
    settings = {}
    for f in required_flags:
        settings[f] = incar_set[f][key_name]

    for f in optional_flags:
        try:
            settings[f] = incar_set[f][key_name]
        except KeyError:
            continue

    return settings


def check_keys(d: dict, required: set, optional: set) -> bool:
    """Check if the required keys exist and other keys are included in optional

    Args:
        d (dict):
            Dictionary to be checked.
        required (set):
            Required key set.
        optional (set):
            Optional key set.

    Returns:
        Simply raise KeyError if the condition is not satisfied.
        If succeed, True is returned.
    """
    check_required = set(required)
    for key in d:
        check_required.discard(key)
        if key not in required | optional:
            raise KeyError(f"{key} does not belong.")

    if check_required:
        raise KeyError(f"{check_required} must exist in {d}.")

    return True


def nelect(composition: Composition, potcar: Potcar, charge: int = 0) -> int:
    """Gets the default number of electrons for a given composition. """
    num_nelect = 0
    for pt in potcar:
        num_nelect += composition.element_composition[pt.element] * pt.ZVAL

    return int(num_nelect) - charge


def nbands(composition: Composition, potcar: Potcar) -> int:
    """Calculate total number of bands for unoccupied related properties

    Useful for optical absorption, band structure, and DOS.

    Args:
        composition (Composition): Input composition
        potcar (Potcar): Input Potcar object

    Returns:
        Int of number of bands
    """
    num_bands = \
        sum([composition[c] * (p.nelectrons / 2 + unoccupied_bands[str(c)])
             for c, p in zip(composition, potcar)])

    return ceil(num_bands)


def calc_npar_kpar(num_kpoints: int, num_nodes: int) -> Tuple[int, int]:
    """Recommend suitable NPAR and KPAR params in INCAR from num of kpts.

    The test is performed using 36 core node.

    Args:
        num_kpoints (int): Number of irreducible k-points.
        num_nodes (int): Number of nodes.

    Returns:
        Tuple of recommended NPAR and KPAR.
    """

    kpar_set = {
        1:    [1, 1, 1],
        2:    [2, 2, 2],
        3:    [3, 3, 3],
        4:    [2, 4, 4],
        5:    [1, 2, 2],
        6:    [3, 6, 6],
        7:    [2, 2, 4],
        8:    [4, 8, 8],
        9:    [3, 3, 3],
        10:   [2, 2, 2],
        11:   [2, 2, 2],
        12:   [3, 4, 4],
        13:   [2, 2, 2],
        14:   [2, 2, 2],
        15:   [3, 2, 2],
        16:   [4, 4, 8],
        17:   [2, 4, 4],
        18:   [3, 6, 4],
        19:   [2, 4, 4],
        20:   [2, 4, 4],
        21:   [3, 3, 3],
        22:   [2, 4, 4],
        23:   [2, 4, 4],
        24:   [3, 6, 8],
        25:   [2, 4, 4],
        26:   [2, 4, 4],
        27:   [3, 4, 4],
        28:   [2, 4, 4],
        29:   [2, 4, 4],
        30:   [3, 6, 4],
        31:   [3, 4, 4],
        32:   [3, 4, 8],
        33:   [3, 4, 4],
        34:   [3, 4, 4],
        35:   [3, 4, 4],
        36:   [4, 8, 16],
        48:   [4, 8, 16],
        60:   [4, 8, 16],
        72:   [4, 8, 16],
        None: [4, 8, 16]}

    num_kpt_key = num_kpoints if num_kpoints in kpar_set else None

    if num_nodes == 2:
        kpar = kpar_set[num_kpt_key][1]
    elif num_nodes == 4:
        kpar = kpar_set[num_kpt_key][2]
    else:
        kpar = kpar_set[num_kpt_key][0]

    if kpar == 1:
        npar = 2
    else:
        npar = 1

    return kpar, npar
