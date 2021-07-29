# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np
import seekpath
import spglib
from monty.json import MSONable
from more_itertools import consecutive_groups
from pymatgen.core import Structure, Element
from collections import defaultdict
from itertools import groupby

from tabulate import tabulate
from vise.defaults import defaults
from vise.util.centering import Centering
from vise.util.logger import get_logger
from vise.error import ViseError
from vise.util.bravais_lattice import BravaisLattice

logger = get_logger(__name__)


def cell_to_structure(
        cell: Tuple[List[List[float]], List[List[float]], List[int]]
) -> Structure:
    """
    cell: (Lattice parameters
           [[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]],
          Fractional atomic coordinates in an Nx3 array,
          Z numbers of species in a length N array)
    """
    species = [Element.from_Z(i) for i in cell[2]]
    return Structure(lattice=cell[0], coords=cell[1], species=species)


class StructureSymmetrizer:
    def __init__(self,
                 structure: Structure,
                 symprec: float = defaults.symmetry_length_tolerance,
                 angle_tolerance: float = defaults.symmetry_angle_tolerance,
                 time_reversal: bool = True,
                 band_mesh_distance: float = defaults.band_mesh_distance):
        """Get full information of seekpath band path.

        Note: site properties such as magmom are removed.

        The structures of aP (SG:1, 2), mC (5, 8, 9, 12, 15) and
        oA (38, 39, 40, 41) can be different between spglib and seekpath.
        see Y. Hinuma et al. Comput. Mater. Sci. 128 (2017) 140–184
        -- spglib mC
         6.048759 -3.479491 0.000000
         6.048759  3.479491 0.000000
        -4.030758  0.000000 6.044512
        -- seekpath mC
         6.048759  3.479491  0.000000
        -6.048759  3.479491  0.000000
        -4.030758  0.000000  6.044512
        """
        self.structure = structure.copy()
        if structure.site_properties:
            logger.warning(f"Site property {structure.site_properties.keys()} "
                           f"removed in primitive and conventional structures.")
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.time_reversal = time_reversal
        self.ref_distance = band_mesh_distance
        lattice_matrix = structure.lattice.matrix
        positions = structure.frac_coords.tolist()
        atomic_numbers = [i.specie.number for i in structure.sites]
        self.cell = (lattice_matrix, positions, atomic_numbers)
        # evaluated lazily
        self._spglib_sym_data = None
        self._primitive = None
        self._second_primitive = None
        self._seekpath_data = None
        self._band_primitive = None

    def __repr__(self):
        sym_data = self.spglib_sym_data
        is_primitive = self.structure == self.primitive
        lines = [f"Symprec: {self.symprec}",
                 f"Angle tolerance: {self.angle_tolerance}",
                 f"Space group: {sym_data['international']}",
                 f"Is primitive: {is_primitive}"]
        site_info_header = ["site", "wyckoff", "site sym", "equiv sites"]
        site_infos = []
        for name, site in self.sites.items():
            site_infos.append([f"{name}",
                               f"{site.wyckoff_letter}",
                               f"{site.site_symmetry}",
                               f"{site.pprint_equiv_atoms}"])
        lines.append(tabulate(site_infos, headers=site_info_header))

        return "\n".join(lines)

    @property
    def spglib_sym_data(self) -> dict:
        if not self._spglib_sym_data:
            self._spglib_sym_data = spglib.get_symmetry_dataset(
                self.cell, self.symprec, self.angle_tolerance)
        return self._spglib_sym_data

    @property
    def conventional(self) -> Structure:
        center = Centering.from_string(self.centering)
        return self.primitive * center.primitive_to_conv

        # # I don't know if this is fine for spglib cyclic behavior.
        # if self._conventional is None:
        #     conventional = \
        #         spglib.standardize_cell(self.cell,
        #                                 symprec=self.symprec,
        #                                 angle_tolerance=self.angle_tolerance)
        #     if conventional is None:
        #         raise ViseSymmetryError(
        #             "Spglib couldn't find the conventional cell. Change the "
        #             "symprec and/or angle_tolerance.")
        #     else:
        #         self._conventional = \
        #             cell_to_structure(conventional).get_sorted_structure()
        # return self._conventional

    @property
    def primitive(self) -> Structure:
        if self._primitive is None:
            primitive = \
                spglib.find_primitive(self.cell,
                                      symprec=self.symprec,
                                      angle_tolerance=self.angle_tolerance)
            if primitive is None:
                raise ViseSymmetryError(
                    "Spglib couldn't find the primitive cell. "
                    "Change the symprec and/or angle_tolerance.")
            else:
                # To manage spglib cyclic behavior, we need to run this again.
                second_primitive = cell_to_structure(spglib.find_primitive(
                    primitive, symprec=self.symprec,
                    angle_tolerance=self.angle_tolerance)).get_sorted_structure()
                primitive = cell_to_structure(primitive)

                if primitive != second_primitive:
                    if first_structure_is_primitive(primitive, second_primitive):
                        self._primitive = primitive
                        self._second_primitive = second_primitive

                    else:
                        self._primitive = second_primitive
                        self._second_primitive = primitive
                else:
                    self._primitive = primitive

        return self._primitive

    @property
    def second_primitive(self) -> Structure:
        return self._second_primitive

    def find_seekpath_data(self) -> None:
        """Get full information of seekpath band path. """
        self._seekpath_data = \
            seekpath.get_explicit_k_path(structure=self.cell,
                                         symprec=self.symprec,
                                         angle_tolerance=self.angle_tolerance,
                                         with_time_reversal=self.time_reversal,
                                         reference_distance=self.ref_distance)
        lattice = self._seekpath_data["primitive_lattice"]
        element_types = self._seekpath_data["primitive_types"]
        species = [Element.from_Z(i) for i in element_types]
        positions = self._seekpath_data["primitive_positions"]
        self._band_primitive = Structure(lattice, species, positions)

    @property
    def sg_number(self):
        return self.spglib_sym_data["number"]

    @property
    def space_group(self):
        return self.spglib_sym_data["international"]

    @property
    def point_group(self):
        return self.spglib_sym_data["pointgroup"]

    @property
    def seekpath_data(self):
        if self._seekpath_data is None:
            self.find_seekpath_data()
        return self._seekpath_data

    @property
    def band_primitive(self) -> Structure:
        if self._band_primitive is None:
            self.find_seekpath_data()
        return self._band_primitive

    @property
    def is_primitive_lattice_changed(self) -> bool:
        # np.allclose is used for lattice comparison.
        # def allclose(a, b, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.structure.lattice != self.primitive.lattice

    @property
    def is_band_primitive_lattice_changed(self) -> bool:
        return self.structure.lattice != self.band_primitive.lattice

    @property
    def band_primitive_differ_primitive(self) -> bool:
        return self.primitive.lattice != self.band_primitive.lattice

    def irreducible_kpoints(self,
                            num_kpt_list: List[int],
                            kpt_shift: List[float]  # vasp convention
                            ) -> List[Tuple[List[float], int]]:
        # evaluated at a given Brillouin zone.
        mapping, integer_grid_points = spglib.get_ir_reciprocal_mesh(
            mesh=num_kpt_list,
            cell=self.cell,
            is_shift=np.array(kpt_shift) * 2,
            is_time_reversal=self.time_reversal,
            symprec=self.symprec)

        results = []
        shift_in_integer_grid = np.array(kpt_shift)
        for repr_index, count in zip(*np.unique(mapping, return_counts=True)):
            integer_grid_point = \
                integer_grid_points[repr_index] + shift_in_integer_grid
            normalized_grid_point = integer_grid_point / num_kpt_list
            results.append((normalized_grid_point.tolist(), count))

        return results  # [(irreducible kpoint in fractional coords, weight)]

    def grouped_atom_indices(self):
        result = {}
        for name, site in self.sites.items():
            key = "_".join([name, site.wyckoff_letter])
            result[key] = site.equivalent_atoms
        return result

    @property
    def sites(self) -> Dict[str, "Site"]:
        wyckoffs = self.spglib_sym_data["wyckoffs"]
        equivalent_atoms = self.spglib_sym_data["equivalent_atoms"]
        site_symmetries = self.spglib_sym_data["site_symmetry_symbols"]
        equiv_indices = sorted(enumerate(equivalent_atoms), key=lambda x: x[1])
        result = {}
        element_idx_dict = defaultdict(int)
        for _, equiv_sites in groupby(equiv_indices, lambda x: x[1]):
            equiv_site_list = list(equiv_sites)
            repr_idx = equiv_site_list[0][0]
            element = self.structure[repr_idx].specie.name
            element_idx_dict[element] += 1
            index = str(element_idx_dict[str(element)])
            name = element + index
            result[name] = Site(element=element,
                                wyckoff_letter=wyckoffs[repr_idx],
                                site_symmetry=site_symmetries[repr_idx],
                                equivalent_atoms=[s[0] for s in equiv_site_list])
        return result

    @property
    def bravais(self):
        return BravaisLattice.from_sg_num(self.spglib_sym_data["number"])

    @property
    def centering(self):
        return self.spglib_sym_data["international"][0]


def first_structure_is_primitive(structure1: Structure, structure2: Structure):
    for s1, s2 in zip(structure1, structure2):
        for c1, c2 in zip(s1.frac_coords, s2.frac_coords):
            if c1 < c2 - 1e-5:
                return True
            elif c2 < c1 - 1e-5:
                return False
    raise ViseSymmetryError("Primitive cannot be unique.")


@dataclass(frozen=True)
class Site(MSONable):
    element: str
    wyckoff_letter: str
    site_symmetry: str
    equivalent_atoms: List[int]

    @property
    def pprint_equiv_atoms(self):
        str_list = []
        for consecutive_ints in consecutive_groups(self.equivalent_atoms):
            ints = list(consecutive_ints)
            if len(ints) >= 3:
                str_list.append("..".join([str(ints[0]), str(ints[-1])]))
            else:
                str_list.append(" ".join([str(j) for j in ints]))
        return " ".join(str_list)


num_sym_op = {"1":     1,
              "-1":    2,
              "2":     2,
              "m":     2,
              "2/m":   4,
              "222":   4,
              "2mm":   4,
              "m2m":   4,
              "mm2":   4,
              "mmm":   8,
              "4":     4,
              "-4":    4,
              "4/m":   8,
              "422":   8,
              "4mm":   8,
              "-4m2":  8,
              "-42m":  8,
              "4/mmm": 16,
              "3":     3,
              "-3":    6,
              "32":    6,
              "3m":    6,
              "-3m":   12,
              "6":     6,
              "-6":    6,
              "6/m":   12,
              "622":   12,
              "6mm":   12,
              "-6m2":  12,
              "6/mmm": 24,
              "23":    12,
              "m3":    24,
              "m-3":   24,
              "432":   24,
              "-43m":  24,
              "m-3m":  48}


def num_symmetry_operation(point_group: str) -> int:
    """ Return number of symmetry operations from Hermann–Mauguin notation.

    Args:
         point_group (str):
            Point group in Hermann–Mauguin notation.
            "." will be removed, e.g., ..6 -> 6.

    Returns:
         Number of symmetry operation.
    """
    # remove "." from the given point_group
    point_group = "".join([s for s in point_group if s != "."])
    return num_sym_op[point_group]


class ViseSymmetryError(ViseError):
    pass
