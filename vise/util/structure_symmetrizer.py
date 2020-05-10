# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List, Tuple

import numpy as np
import seekpath
import spglib
from pymatgen import Structure, Element
from collections import defaultdict
from itertools import groupby

from vise.defaults import defaults
from vise.util.logger import get_logger
from vise.error import ViseError

logger = get_logger(__name__)


def cell_to_structure(
        cell: Tuple[List[List[float]], List[List[float]], List[int]]
) -> Structure:
    # cell: (Lattice parameters
    #        [[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]],
    #       Fractional atomic coordinates in an Nx3 array,
    #       Z numbers of species in a length N array)
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
        self._conventional = None
        self._primitive = None
        self._seekpath_data = None
        self._band_primitive = None
        self._irreducible_kpoints = None

    @property
    def spglib_sym_data(self) -> dict:
        if not self._spglib_sym_data:
            self._spglib_sym_data = spglib.get_symmetry_dataset(
                self.cell, self.symprec, self.angle_tolerance)
        return self._spglib_sym_data

    @property
    def conventional(self) -> Structure:
        if self._conventional is None:
            conventional = \
                spglib.standardize_cell(self.cell,
                                        symprec=self.symprec,
                                        angle_tolerance=self.angle_tolerance)
            if conventional is None:
                raise ViseSymmetryError(
                    "Spglib couldn't find the conventional cell. Change the "
                    "symprec and/or angle_tolerance.")
            else:
                self._conventional = cell_to_structure(conventional)
        return self._conventional

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
                self._primitive = cell_to_structure(primitive)
        return self._primitive

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
        wyckoffs = self.spglib_sym_data["wyckoffs"]
        site_symmetries = self.spglib_sym_data["site_symmetry_symbols"]
        equivalent_atom_set = \
            enumerate(self.spglib_sym_data["equivalent_atoms"].tolist())
        elem_wyckoff_idx = defaultdict(int)

        for index, same_sites_index_list in groupby(equivalent_atom_set,
                                                    key=lambda x: x[1]):
            elem = str(self.structure[index].specie)
            wyckoff = wyckoffs[index]
            elem_wyckoff = "_".join([elem, wyckoff])
            elem_wyckoff_idx[elem_wyckoff] += 1
            name = elem_wyckoff + str(elem_wyckoff_idx[elem_wyckoff])
            result[name] = [i[0] for i in list(same_sites_index_list)]

        return result


class ViseSymmetryError(ViseError):
    pass
