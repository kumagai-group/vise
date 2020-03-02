# -*- coding: utf-8 -*-

from itertools import groupby
import operator
from typing import Tuple

import seekpath
import spglib

from pymatgen import Structure
from pymatgen.core.periodic_table import Element

from vise.config import ANGLE_TOL, SYMMETRY_TOLERANCE, BAND_REF_DIST
from vise.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_symmetry_dataset(structure: Structure,
                         symprec: float = SYMMETRY_TOLERANCE,
                         angle_tolerance: float = ANGLE_TOL) -> dict:
    """ Get spglib symmetry dataset from a Structure.

    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.

    Return:
        spglib symmetry dataset. See docstrings of spglib for details.
    """
    cell = structure_to_spglib_cell(structure)
    return spglib.get_symmetry_dataset(cell=cell,
                                       symprec=symprec,
                                       angle_tolerance=angle_tolerance)


def structure_to_spglib_cell(structure: Structure) -> Tuple[list, list, list]:
    """Convert structure to spglib cell tuple.

    Return a cell tuple used by spglib that is composed of lattice parameters,
    atomic positions in fractional coordinates, and atomic numbers.

    Args:
        structure (Structure):
            Pymatgen Structure class object

    Returns:
        Tuple of (lattice matrix, fractional coordinates, atomic numbers
    """
    lattice = list(structure.lattice.matrix)
    positions = structure.frac_coords.tolist()
    atomic_numbers = [i.specie.number for i in structure.sites]

    return lattice, positions, atomic_numbers


def spglib_cell_to_structure(cell: tuple) -> Structure:
    """
    Args:
        cell (3 tuple):
            Lattice parameters, atomic positions in fractional coordinates,
            and corresponding atom numbers
    Returns:
         pymatgen Structure class object from spglib cell tuple.
    """
    species = [Element.from_Z(i) for i in cell[2]]
    return Structure(cell[0], species, cell[1])


def find_spglib_standard_conventional(structure: Structure,
                                      symprec: float = SYMMETRY_TOLERANCE,
                                      angle_tolerance: float = ANGLE_TOL
                                      ) -> Structure:
    """

    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.
    Returns:
        Structure of a standard conventional unit cell.
    """
    cell = structure_to_spglib_cell(structure)
    return spglib_cell_to_structure(
        spglib.standardize_cell(cell=cell,
                                to_primitive=False,
                                no_idealize=False,
                                symprec=symprec,
                                angle_tolerance=angle_tolerance))


def find_spglib_primitive(structure: Structure,
                          symprec: float = SYMMETRY_TOLERANCE,
                          angle_tolerance: float = ANGLE_TOL
                          ) -> Tuple[Structure, bool]:
    """
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.

    Returns:
        Structure of a primitive unit cell.
    """
    cell = structure_to_spglib_cell(structure)
    primitive_cell = spglib.find_primitive(cell=cell,
                                           symprec=symprec,
                                           angle_tolerance=angle_tolerance)
    primitive_structure = spglib_cell_to_structure(primitive_cell)
    is_structure_changed = structure.lattice != primitive_structure.lattice

    return primitive_structure, is_structure_changed


def find_hpkot_primitive(structure: Structure,
                         symprec: float = SYMMETRY_TOLERANCE,
                         angle_tolerance: float = ANGLE_TOL
                         ) -> Structure:
    """
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.

    Returns:
         Structure of a hpkot primitive unit cell.
    """
    cell = structure_to_spglib_cell(structure)
    res = seekpath.get_explicit_k_path(structure=cell,
                                       symprec=symprec,
                                       angle_tolerance=angle_tolerance)

    return seekpath_to_hpkot_structure(res)


def structure_to_seekpath(structure: Structure,
                          time_reversal: bool = True,
                          ref_distance: float = BAND_REF_DIST,
                          recipe: str = 'hpkot',
                          threshold: float = 1e-7,
                          symprec: float = SYMMETRY_TOLERANCE,
                          angle_tolerance: float = ANGLE_TOL) -> dict:
    """Get full information of seekpath band path.

    Args:
        structure (Structure):
            Pymatgen Structure class object
        time_reversal (bool):
            If the time reversal symmetry exists
        ref_distance (float):
            Distance for the k-point mesh.
        recipe (str):
            See docstrings of seekpath.
        threshold (float):
            To use to verify if we are in edge case (see seekpath).
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analyzer.

    Return:
        Dict with some properties. See docstrings of seekpath.
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
        logger.warning(
            "Given structure is modified to be compatible with HPKOT k-path.")

    return res


def seekpath_to_hpkot_structure(res: dict) -> Structure:
    """
    Args:
        res (dict):
            seekpath res dictionary.

    Returns:
         pymatgen Structure obtained from seekpath res dictionary.
    """
    lattice = res["primitive_lattice"]
    element_types = res["primitive_types"]
    species = [Element.from_Z(i) for i in element_types]
    positions = res["primitive_positions"]

    return Structure(lattice, species, positions)


def get_symbol_list(structure: Structure):
    species = [str(s) for s in structure.species]
    # unique_justseen https://docs.python.org/ja/3/library/itertools.html
    # ["H", "H", "O", "O", "H"] -> ['H', 'O', 'H']
    return list(map(next, map(operator.itemgetter(1), groupby(species, None))))
