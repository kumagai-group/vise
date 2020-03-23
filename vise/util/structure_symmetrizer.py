# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen import Structure, Element

import seekpath

import spglib

from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL, BAND_MESH_DISTANCE
from vise.util.logger import get_logger


logger = get_logger(__name__)


def cell_to_structure(cell: tuple) -> Structure:
    """Structure is returned from spglib cell """
    return Structure(cell[0], [Element.from_Z(i) for i in cell[2]], cell[1])


class StructureSymmetrizer:

    def __init__(self,
                 structure: Structure,
                 symprec: float = SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL):
        """Get full information of seekpath band path.

        Note: site properties such as magmom are removed.

        Args:
            structure (Structure):
                Pymatgen Structure class object
            symprec (float):
                Distance tolerance in cartesian coordinates Unit is compatible
                with the structure.
            angle_tolerance (float):
                Angle tolerance used for symmetry analyzer.
        """
        self.structure = structure
        props = self.structure.site_properties
        if props:
            logger.warning(f"The site properties {props.keys()} are removed"
                           f"in the primitive and conventional structures.")
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        lattice = list(structure.lattice.matrix)
        positions = structure.frac_coords.tolist()
        atomic_numbers = [i.specie.number for i in structure.sites]
        self.cell = lattice, positions, atomic_numbers
        self.kwargs = {"cell": self.cell,
                       "symprec": symprec,
                       "angle_tolerance": angle_tolerance}
        self._spglib_sym_data = None
        self._conventional = None
        self._primitive = None
        self._seekpath_data = None
        self._band_primitive = None

    @property
    def spglib_sym_data(self):
        if not self._spglib_sym_data:
            self._spglib_sym_data = spglib.get_symmetry_dataset(**self.kwargs)
        return self._spglib_sym_data

    @property
    def conventional(self):
        if self._conventional is None:
            conventional = spglib.standardize_cell(**self.kwargs)
            if conventional is None:
                raise ViseSymmetryError(
                    "Spglib couldn't find the conventional cell. Change the "
                    "symprec and/or angle_tolerance.")
            else:
                self._conventional = cell_to_structure(conventional)
        return self._conventional

    @property
    def primitive(self):
        if self._primitive is None:
            primitive = spglib.find_primitive(**self.kwargs)
            if primitive is None:
                raise ViseSymmetryError(
                    "Spglib couldn't find the primitive cell. Change the "
                    "symprec and/or angle_tolerance.")
            else:
                self._primitive = cell_to_structure(primitive)
        return self._primitive

    def find_seekpath_data(self,
                           time_reversal: bool = True,
                           ref_distance: float = BAND_MESH_DISTANCE):
        """Get full information of seekpath band path.

        Args:
            time_reversal (bool):
                If the time reversal symmetry exists
            ref_distance (float):
                Mesh distance for the k-point mesh.

        Return:
            Dict with some properties. See docstrings of seekpath.
        """
        self._seekpath_data = \
            seekpath.get_explicit_k_path(structure=self.cell,
                                         symprec=self.symprec,
                                         angle_tolerance=self.angle_tolerance,
                                         with_time_reversal=time_reversal,
                                         reference_distance=ref_distance)
        lattice = self._seekpath_data["primitive_lattice"]
        element_types = self._seekpath_data["primitive_types"]
        species = [Element.from_Z(i) for i in element_types]
        positions = self._seekpath_data["primitive_positions"]
        self._band_primitive = Structure(lattice, species, positions)

    @property
    def seekpath_data(self):
        if self._seekpath_data is None:
            raise ViseSymmetryError("seekpath_data is not set.")
        return self._seekpath_data

    @property
    def band_primitive(self):
        if self._band_primitive is None:
            raise ViseSymmetryError("seekpath_primitive is not set.")
        return self._band_primitive

    @property
    def is_lattice_changed(self) -> bool:
        if not hasattr(self, "primitive"):
            raise ViseSymmetryError("Primitive cell is not searched for.")
        # For Lattice comparison, np.allclose is used.
        # def allclose(a, b, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.structure.lattice != self.primitive.lattice

# if self.structure.lattice != lattice:
#     logger.warning(
#         "Given structure is modified consistently with HPKOT k-path.")


class ViseSymmetryError(Exception):
    """Raised when the spglib return is inadequate."""
    pass
