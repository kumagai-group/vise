# -*- coding: utf-8 -*-
import math
import warnings

import numpy as np
import seekpath
import spglib

from pymatgen import Structure
from pymatgen.io.vasp import Poscar

from obadb.database.atom import symbols_to_atom

from vise.core.config import ANGLE_TOL, SYMPREC

from atomate.utils.utils import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_symmetry_dataset(structure: Structure,
                         symprec: float = SYMPREC,
                         angle_tolerance: float = ANGLE_TOL) -> dict:
    """
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.
    """
    cell = structure_to_spglib_cell(structure)
    return spglib.get_symmetry_dataset(cell, symprec=symprec,
                                       angle_tolerance=angle_tolerance)


# def get_point_group_from_dataset(sym_dataset: dict,
#                                  coords: list,
#                                  lattice: np.ndarray,
#                                  symprec: float = SYMPREC) -> tuple:
#     """
#     Args:
#         sym_dataset (dict):
#             spglib get_symmetry_dataset.
#         coords (list):
#             Fractional coordinates.
#         lattice (np.ndarray):
#             3x3 numpy ndarray
#         symprec (float):
#             Distance tolerance in cartesian coordinates Unit is compatible with
#             the cell.
#     """
#     full_rotations = sym_dataset["rotations"]
#     translations = sym_dataset["translations"]
#     rotations = get_rotations(coords, lattice, full_rotations, translations,
#                               symprec)
#     return get_point_group_from_rotations(rotations)


def get_point_group_from_rotations(rotations):
    ptg = spglib.get_pointgroup(rotations)
    return ptg[0].strip(), ptg[2]


def get_rotations(coords, lattice, rotations, translations, symprec=SYMPREC):
    """
    Args:
        coords (list):
            Cartesian coordinates.
        lattice (numpy.array):
            3x3 numpy array
        rotations (dict):
            list of 3x3 rotation matrix.
        translations (dict):
            list of 3 translation column.
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
    """
    site_symmetries = []

    for r, t in zip(rotations, translations):
        rot_pos = np.dot(coords, r.T) + t
        diff = coords - rot_pos
        diff -= np.rint(diff)
        diff = np.dot(diff, lattice)
        if np.linalg.norm(diff) < symprec:
            site_symmetries.append(r)

    return np.array(site_symmetries, dtype='intc')


def structure_to_spglib_cell(structure):
    """
    Return a *cell* tuple parsed by spglib that is composed of lattice
    parameters, atomic positions in fractional coordinates, and corresponding
    atomic numbers.
    Args:
        structure (Structure):
            Pymatgen Structure class object
    """
    lattice = list(structure.lattice.matrix)
    positions = structure.frac_coords.tolist()
    atomic_numbers = [i.specie.number for i in structure.sites]
    return lattice, positions, atomic_numbers


def spglib_cell_to_structure(cell):
    """
    Return a pymatgen Structure class object from spglib cell tuple.
    Args:
        cell (3 tuple):
            Lattice parameters, atomic positions in fractional coordinates,
            and corresponding atom numbers
    """
    species = [symbols_to_atom[i] for i in cell[2]]
    return Structure(cell[0], species, cell[1])


def find_spglib_standard_conventional(structure, symprec=SYMPREC,
                                      angle_tolerance=ANGLE_TOL):
    """
    Return a standard conventional unit cell.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.
    """
    cell = structure_to_spglib_cell(structure)
    return spglib_cell_to_structure(
        spglib.standardize_cell(cell, to_primitive=False, no_idealize=False,
                                symprec=symprec,
                                angle_tolerance=angle_tolerance))


def find_spglib_standard_primitive(structure: Structure,
                                   symprec: float = SYMPREC,
                                   angle_tolerance: float = ANGLE_TOL):
    """
    Return a primitive unit cell.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.

    """
    cell = structure_to_spglib_cell(structure)
    primitive_structure = \
        spglib_cell_to_structure(
            spglib.find_primitive(cell, symprec=symprec,
                                  angle_tolerance=angle_tolerance))
    is_structure_changed = structure.lattice != primitive_structure.lattice

    return primitive_structure, is_structure_changed


def find_hpkot_primitive(structure, symprec=SYMPREC,
                         angle_tolerance=ANGLE_TOL):
    """
    Return a hpkot primitive unit cell.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.
    """
    cell = structure_to_spglib_cell(structure)
    res = seekpath.get_explicit_k_path(structure=cell, symprec=symprec,
                                       angle_tolerance=angle_tolerance)

    return seekpath_to_hpkot_structure(res)


def structure_to_seekpath(structure, time_reversal=True, ref_distance=0.025,
                          recipe='hpkot', threshold=1.e-7, symprec=SYMPREC,
                          angle_tolerance=ANGLE_TOL):
    """
    Return the full information for the band path of the given Structure class
    object generated by seekpath.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        time_reversal (bool):
            If the time reversal symmetry exists
        ref_distance (float):
            Distance for the k-point mesh.
        threshold (float):
            To use to verify if we are in edge case (see seekpath).
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.
    """
    cell = structure_to_spglib_cell(structure)
    res = seekpath.get_explicit_k_path(cell,
                                       with_time_reversal=time_reversal,
                                       reference_distance=ref_distance,
                                       recipe=recipe,
                                       threshold=threshold,
                                       symprec=symprec,
                                       angle_tolerance=angle_tolerance)

    # If numpy.allclose is too strict in pymatgen.core.lattice __eq__,
    # make almost_equal
    if structure.lattice != seekpath_to_hpkot_structure(res).lattice:
        warnings.warn("Given structure is modified to be compatible with HPKOT "
                      "k-path.")

    return res

def seekpath_to_hpkot_structure(res):
    """
    Return a pymatgen Structure class object from seekpath res dictionary.
    Args:
        res (dict):
            seekpath res dictionary.
    """
    lattice = res["primitive_lattice"]
    element_types = res["primitive_types"]
    species = [symbols_to_atom[i] for i in element_types]
    positions = res["primitive_positions"]
    return Structure(lattice, species, positions)


def fold_float(x):
    """ Return the folded float number, e.g., 3.5 -> 0.5 and -0.7 -> 0.3. """
    return x - math.floor(x)


def fold_position_structure(structure):
    """ Modify positions out of box (x<0 or x>=1) into box (0 <= x < 1)

    For example, coords of site changes from [-0.3, 1.9, 0.5] to [0.7, 0.9, 0.5]

    Args:
        structure(Structure):
            Input Structure.

    Returns:
        structure(Structure):
            Structure with folded fractional coordinates
    """
    for i, site in enumerate(structure):
        modification_vector = [-math.floor(v) for v in site.frac_coords]
        structure.translate_sites(i, modification_vector)
    return structure


def fold_position_poscar(poscar):
    """ Same as fold_position_structure but for POSCAR.

    Args:
        poscar(Poscar):
            Input POSCAR.

    Returns:
        poscar(Poscar):
            POSCAR with folded fractional coordinates
    """
    s = poscar.structure
    fold_position_structure(s)
    return Poscar(s)



