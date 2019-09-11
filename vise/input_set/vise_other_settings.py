from math import ceil
from pathlib import Path
from typing import Optional

import numpy as np
from monty.serialization import loadfn
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Potcar, Kpoints
from vise.input_set.datasets.element_parameters import unoccupied_bands
from vise.input_set.make_kpoints import make_kpoints
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.structure_handler import find_spglib_primitive

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

    Return:
          Dictionary of potcar_set, like {"Zr": "Zr_pv", ...}
    """
    potcar_set = loadfn(SET_DIR / "potcar_set.yaml")
    try:
        potcar = potcar_set[set_name]
    except KeyError:
        logger.warning(f"The accepted potcar set name is {potcar_set.keys()}")
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
            Required INCAR flags.
        optional_flags (set):
            Optional INCAR flags.
        key_name (str):
            Key name such as "structure_opt" or "hse".

    Return:
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

    Return:
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
    """Gets the default number of electrons for a given structure.


     """
    # if structure is not sorted this can cause problems, so must take
    # care to remove redundant symbols when counting electrons
    num_nelect = 0
    for pt in potcar:
        num_nelect += composition.element_composition[pt.element] * pt.ZVAL

    return int(num_nelect) - charge


def nbands(composition: Composition, potcar: Potcar) -> int:
    """
    Calculate the total number of bands required for the unoccupied related
    properties such as optical absorption, band structure, and DOS.

    Args:
        composition:
        potcar (Potcar):
    """
    num_bands = \
        sum([composition[c] * (p.nelectrons / 2 + unoccupied_bands[str(c)])
             for c, p in zip(composition, potcar)])

    return ceil(num_bands)


def calc_npar_kpar(num_kpoints, num_cores_per_node, num_nodes):
    """
    :param num_kpoints:
    :param num_cores_per_node:
    :param num_nodes:
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

    return kpar, npar, num_kpoints


class TaskStructureKpoints:

    def __init__(self,
                 structure: Structure,
                 kpoints: Kpoints,
                 is_structure_changed: bool,
                 sg: int,
                 num_kpts: int,
                 factor: int):

        self.structure = structure
        self.kpoints = kpoints
        self.is_structure_changed = is_structure_changed
        self.sg = sg
        self.num_kpts = num_kpts
        self.factor = factor

    @classmethod
    def from_options(cls,
                     task: Task,
                     original_structure: Structure,
                     standardize_structure: bool,
                     sort_structure: bool,
                     is_magnetization: bool,
                     kpt_mode: str,
                     kpt_density: float,
                     kpt_shift: list,
                     only_even: bool,
                     band_ref_dist: float,
                     factor: Optional[int],
                     symprec: float,
                     angle_tolerance: float):
        """

        See docstrings in the make_input method in ViseInputSet class
        :return:
        """
        if sort_structure:
            structure = original_structure.get_sorted_structure()
            if structure.symbol_set != original_structure.symbol_set:
                logger.warning(
                    "CAUTION: The sequence of the species is changed."
                    f"Symbol set in the original structure " 
                    f"{original_structure.symbol_set} "
                    f"Symbol set in the generated structure "
                    f"{structure.symbol_set}")
        else:
            structure = original_structure.copy()

        is_structure_changed = False
        if task == Task.defect:
            kpt_mode = "manual_set"
        elif task == Task.cluster_opt:
            kpt_mode = "manual_set"
            kpt_density = 1e-5
            only_even = False
            kpt_shift = [0, 0, 0]
        elif task == Task.band:
            kpt_mode = "band"
        elif task == Task.phonon_force:
            kpt_mode = "manual_set"
            if kpt_shift != [0, 0, 0]:
                logger.warning("For phonon force calculaitons, Gamma centering "
                               "is forced for k-point sampling.")
            kpt_shift = [0, 0, 0]
        else:
            primitive_structure, is_structure_changed = \
                find_spglib_primitive(structure, symprec, angle_tolerance)

            if standardize_structure:
                org = original_structure.lattice.matrix
                primitive = primitive_structure.lattice.matrix
                if is_structure_changed:
                    with np.printoptions(precision=3, suppress=True):
                        logger.warning(
                            "CAUTION: The structure is changed.\n"
                            f"Original lattice\n {org} \n"
                            f"Generated lattice\n {primitive} \n")
                structure = primitive_structure
            else:
                if is_structure_changed and kpt_mode != "manual_set":
                    logger.warning(
                        "Standardizaion is set to False and the given "
                        "structure is not a primitive cell. Thus, the "
                        "kpoint set is switched to manual_set.")
                    kpt_mode = "manual_set"
                else:
                    logger.info("The structure is a standardized primitive.")

        if factor is None:
            if task == Task.dielectric_function:
                factor = 3
            elif task in (Task.dos, Task.dielectric_dfpt,
                          Task.dielectric_finite_field):
                factor = 2
            else:
                factor = 1
        print(kpt_density)
        # - KPOINTS construction
        kpoints, structure, sg, num_kpts = \
            make_kpoints(mode=kpt_mode,
                         structure=structure,
                         kpt_density=kpt_density,
                         only_even=only_even,
                         ref_distance=band_ref_dist,
                         kpt_shift=kpt_shift,
                         factor=factor,
                         symprec=symprec,
                         angle_tolerance=angle_tolerance,
                         is_magnetization=is_magnetization)

        return cls(structure, kpoints, is_structure_changed, sg, num_kpts,
                   factor)


class XcTaskPotcar:

    def __init__(self,
                 potcar: Potcar,
                 max_enmax: int):

        self.potcar = potcar
        self.max_enmax = max_enmax

    @classmethod
    def from_options(cls,
                     xc: Xc,
                     symbol_set: tuple,
                     potcar_set_name: str = None,
                     override_potcar_set: dict = None):

        potcar_functional = "LDA" if xc == Xc.lda else "PBE"
        potcar_set_name = potcar_set_name or "normal"

        potcar_list = load_potcar_yaml(potcar_set_name, override_potcar_set)
        potcar_symbols = [potcar_list.get(el, el) for el in symbol_set]
        potcar = Potcar(potcar_symbols, functional=potcar_functional)

        max_enmax = max([p.enmax for p in potcar])

        return cls(potcar, max_enmax)


