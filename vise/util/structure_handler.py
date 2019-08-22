# -*- coding: utf-8 -*-
import math
from collections import defaultdict
import warnings

import numpy as np
import seekpath
import spglib

from pymatgen import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.core.periodic_table import DummySpecie

from obadb.database.atom import symbols_to_atom
from obadb.database.atom import charge as charge_list
from obadb.atomate.vasp.config import SYMMETRY_TOLERANCE as SYMPREC
from obadb.atomate.vasp.config import ST_MATCHER_ANGLE_TOL

from atomate.utils.utils import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_symmetry_dataset(structure, symprec=SYMPREC,
                         angle_tolerance=ST_MATCHER_ANGLE_TOL):
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


def get_point_group_from_dataset(sym_dataset, coords, lattice, symprec=SYMPREC):
    """
    Args:
        sym_dataset (dict):
            spglib get_symmetry_dataset.
        coords (list):
            Fractional coordinates.
        lattice (numpy.array):
            3x3 numpy array
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
    """
    full_rotations = sym_dataset["rotations"]
    translations = sym_dataset["translations"]
    rotations = get_rotations(coords, lattice, full_rotations, translations,
                              symprec)
    return get_point_group_from_rotations(rotations)


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
                                      angle_tolerance=ST_MATCHER_ANGLE_TOL):
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


def find_spglib_standard_primitive(structure, symprec=SYMPREC,
                                   angle_tolerance=ST_MATCHER_ANGLE_TOL):
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
                         angle_tolerance=ST_MATCHER_ANGLE_TOL):
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
                          angle_tolerance=ST_MATCHER_ANGLE_TOL):
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


def get_coordination_environment(structure, index, factor=1.3):
    """

    Since the coordination is determined by the sum of the ionic radii,
    the coordination may be different around cation and anion.
    For example, in ZnS, when factor
        Zn:  S: 2.36 2.36 2.36 2.36
        S: Zn: 2.36 2.36 2.36 2.36 S: 3.85 3.85 3.85 3.85 3.85 3.85 3.85 3.85
                                      3.85 3.85 3.85 3.85

    Args:
        structure (Structure):
            pmg Structure class object
        index (int):
            The atomic index
        factor (float):
            Multiple number of the distance of the sum of the averaged ionic
            radii that is used to detect the coordination.

    Return:
        coords (dict):
            values are tuples of Element object and distance.
            example  {'O': [(PeriodicSite: O (0.0000, 0.0000, -2.1234)
            [-0.5000, -0.5000, 0.5000], 2.123447), ...
    """
    ionic_radii = {}
    for e in structure.types_of_specie:
        if isinstance(e, DummySpecie):
            ionic_radii[e] = 1.5
            logger.warning("Use atomic radius for {} of {}".
                           format(str(e), ionic_radii[e]))
        else:
            try:
                charge = e.oxi_state
            except AttributeError:
                charge = charge_list[str(e)]
            try:
                ionic_radii[e] = e.ionic_radii[charge]
            except KeyError or TypeError:
                ionic_radii[e] = e.atomic_radius * 1.2
                logger.warning("20% increased atomic radius {} is used for "
                               "ionic radius of {}".
                               format(ionic_radii[e], str(e)))

    coords = []
    first_ionic_radius = ionic_radii[structure.sites[index].specie]
    for e in structure.types_of_specie:
        second_ionic_radius = ionic_radii[e]
        cutoff = (first_ionic_radius + second_ionic_radius) * factor
        neighbors = structure.get_neighbors(structure.sites[index], cutoff)
        for site in neighbors:
            if site[0].specie == e:
                coords.append(site)
    
    return coords


def get_coordination_distances(structure, index, factor=1.3):
    """
    Return:
        coords (dict):
            example {"Mg": [1.92, 1.73], "Al": [2.01, 2.11]}
    """
    coordination_environment = \
        get_coordination_environment(structure, index, factor)
    coordination_distances = defaultdict(list)

    for ce in coordination_environment:
        # reomve oxidation state from species_string
        specie = ''.join([i for i in ce[0].species_string if not (i.isdigit() or i == "+" or i == "-")])
        coordination_distances[specie].append(round(float(ce[1]), 2))

    for k, v in coordination_distances.items():
        coordination_distances[k] = sorted(v) 

    return dict(coordination_distances)


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



