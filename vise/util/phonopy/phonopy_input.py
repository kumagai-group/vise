# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Optional, Dict, Any

import numpy as np
from monty.json import MSONable

try:
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp_from_strings, \
        parse_force_constants, get_born_vasprunxml
    from phonopy.structure.atoms import PhonopyAtoms
    from pymatgen.electronic_structure.boltztrap2 import VasprunBSLoader, \
        BztInterpolator, BztTransportProperties
except ImportError:
    raise ImportError('Calculating effective mass requires BoltzTrap2')

from pymatgen.core import IStructure, Structure
from pymatgen.io.vasp import Vasprun
from vise.util.bravais_lattice import BravaisLattice
from vise.util.centering import Centering
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import StructureSymmetrizer


logger = get_logger(__name__)


@dataclass
class PhonopyInput(MSONable, ToJsonFileMixIn):
    """Parameters for non-analytical term correction

    dict
        Parameters used for non-analytical term correction
        'born': ndarray
            Born effective charges
            shape=(primitive cell atoms, 3, 3), dtype='double', order='C'
        'factor': float
            Unit conversion factor
        'dielectric': ndarray
            Dielectric constant tensor
            shape=(3, 3), dtype='double', order='C'

    """
    unitcell: Structure
    supercell: Structure
    supercell_matrix: List[List[float]]
    force_constants: Optional[np.ndarray] = None
    nac_params: Optional[Dict[str, Any]] = None

    @property
    def to_phonopy(self) -> Phonopy:
        if self.force_constants is None:
            raise ValueError("Set force_constants first.")
        _nac = "" if self.nac_params else "*not* "
        logger.info(f"Parameters for a non-analytical term correction are "
                    f"{_nac}set.")

        result = Phonopy(unitcell=structure_to_phonopy_atoms(self.unitcell),
                         supercell_matrix=self.supercell_matrix,
                         nac_params=self.nac_params)
        result.force_constants = self.force_constants
        return result

    def set_force_constants_from_vasprun(self, vasprun_name: str) -> None:
        self.force_constants, _ = parse_force_constants(vasprun_name)
        vasprun = Vasprun(vasprun_name)
        if vasprun.vasp_version[0] == "6":
            self.force_constants /= 15.633302 ** 2

    def set_nac_params_from_vasprun(self, vasprun_name: str) -> None:
        born, epsilon, _, = get_born_vasprunxml(
            vasprun_name, primitive_matrix=np.eye(3),
            supercell_matrix=self.supercell_matrix)
        self.nac_params = {"born": born, "dielectric": epsilon}
        # XX: need factor?


def structure_to_phonopy_atoms(structure: IStructure) -> PhonopyAtoms:
    return read_vasp_from_strings(structure.to(fmt="poscar"), None)


def make_phonopy_input(unitcell: Structure,
                       supercell_matrix: List[List[float]] = None,
                       conventional_base: bool = True) -> PhonopyInput:
    symmetrizer = StructureSymmetrizer(unitcell)
#    assert symmetrizer.is_primitive_lattice_changed
    supercell_matrix = (supercell_matrix
                        or default_supercell_matrix(symmetrizer.bravais))

    if conventional_base:
        centering = Centering(symmetrizer.centering)
        p_to_c_matrix = centering.primitive_to_conv
        supercell_matrix = np.dot(supercell_matrix, p_to_c_matrix)

    supercell = symmetrizer.primitive * supercell_matrix
    return PhonopyInput(unitcell=symmetrizer.primitive, supercell=supercell,
                        supercell_matrix=supercell_matrix.tolist())


def default_supercell_matrix(bravais: BravaisLattice) -> np.ndarray:
    if bravais in [BravaisLattice.hR, BravaisLattice.hP]:
        return np.array([[2, 1, 0], [1, 2, 0], [0, 0, 2]])
    else:
        return np.eye(3) * 2

