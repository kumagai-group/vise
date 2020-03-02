# -*- coding: utf-8 -*-

from enum import Enum, unique
from math import ceil, modf, pow
from typing import List

import numpy as np
from pymatgen import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from vise.config import BAND_REF_DIST, KPT_DENSITY
from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from vise.input_set.datasets.kpt_centering import kpt_centering
from vise.util.logger import get_logger
from vise.util.structure_handler import (
    structure_to_seekpath, find_spglib_primitive, get_symmetry_dataset)

logger = get_logger(__name__)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


@unique
class KpointsMode(Enum):
    """K-point generation type

       Supporting modes are:
       "band":
           Kpoints with the band path will be returned based on the
           seekpath program. The space group is analyzed and primitive
           unitcell that must be used for the band structure calculation is
           returned as well.
       "primitive_uniform":
           Kpoints with uniform k-point sampling. The k-point sampling mesh
           and centering are determined based on the standardized primitive
           unitcell. Structure is also changed if not primitive.
       "manual_set":
           Kpoints with uniform k-point sampling. The k-point sampling mesh
           and centering are determined based on the given lattice. Note
           that only when the angles are 90 degrees, the centering is
           shifted along the perpendicular direction.
           This mode is useful when calculating the supercells.
    """

    band = "band"
    primitive_uniform = "primitive_uniform"
    manual_set = "manual_set"

    def __repr__(self):
        return self.value

    # Don't remove __str__. This is a must.
    def __str__(self):
        return self.value

    @classmethod
    def from_string(cls, s: str):
        for m in cls:
            if m.value == s or m.name == s:
                return m
        raise AttributeError(f"k-point mode {str(s)} is not proper.\n",
                             f"Supported modes:\n {cls.name_list()}")

    @classmethod
    def name_list(cls):
        return ', '.join([e.name for e in cls])


class MakeKpoints:
    """Make Kpoint based on default settings depending on the task.

        # The structures of aP (SG:1, 2), mC (5, 8, 9, 12, 15) and
        # oA (38, 39, 40, 41) are different between spglib and seekpath.
        # see Y. Hinuma et al. Comput. Mater. Sci. 128 (2017) 140–184
        # -- spglib mC
        #  6.048759 -3.479491 0.000000
        #  6.048759  3.479491 0.000000
        # -4.030758  0.000000 6.044512
        # -- seekpath mC
        #  6.048759  3.479491  0.000000
        # -6.048759  3.479491  0.000000
        # -4.030758  0.000000  6.044512
        # -- spglib oA
        #  6.373362  0.000000  0.000000
        #  0.000000  3.200419  5.726024
        #  0.000000 -3.200419  5.726024
        # -- seekpath oA
        #  0.000000  3.200419 -5.726024
        #  0.000000  3.200419  5.726024
        #  6.373362  0.000000  0.000000
    """

    def __init__(self,
                 mode: str,
                 structure: Structure,
                 kpt_density: float = KPT_DENSITY,
                 only_even: bool = False,
                 manual_kpts: list = None,
                 ref_distance: float = BAND_REF_DIST,
                 kpt_shift: list = None,
                 factor: int = 1,
                 symprec: float = SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL,
                 is_magnetization: bool = False):
        """Kpoint object based on default settings depending on the task.

        Note that it does not check if primitive cell is standardized or not.

        Args:
            mode (str):
                KpointsMode in string.
            structure (Structure/IStructure):
                An input Structure object.
            kpt_density (float):
                Density of k-point mesh along each direction.
            only_even (bool):
                Numbered k points are ceiled to be even numbers.
            manual_kpts (3x1 list):
                Manual set of the numbers of k-points.
            ref_distance (float):
                Distance for the k-point mesh.
            kpt_shift (1x3 list):
                K-point shift in the definition of the vasp setting.
            factor (int):
                Multiplier factor. This is useful for the use of NKRED = factor
                for hybrid xc calculations.
            symprec (float):
                Precision in Angstrom used for the symmetry search.
            angle_tolerance (float):
                Angle tolerance used for symmetry analyzer.
            is_magnetization (bool):
                Whether the magnetization is considered or not.
                This modifies the band structure path of systems w/o inversion.
        """
        self.mode = KpointsMode.from_string(mode)
        self.initial_structure = structure
        self.kpt_density = kpt_density
        self.only_even = only_even
        self.manual_kpts = manual_kpts
        self.ref_distance = ref_distance
        self.kpt_shift = kpt_shift
        self.factor = factor
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.is_magnetization = is_magnetization

        self.comment = None
        self.kpt_mesh = None
        self.kpoints = None
        self.num_kpts = None
        self.corresponding_structure = None
        self.sg = None
        self.sg_symbol = None

        if self.mode == KpointsMode.band:
            self.seekpath_info = \
                structure_to_seekpath(structure=structure,
                                      ref_distance=ref_distance,
                                      time_reversal=(not self.is_magnetization),
                                      symprec=symprec,
                                      angle_tolerance=angle_tolerance)
        else:
            self.seekpath_info = None

    def make_kpoints(self):
        self._set_structure()
        self._set_kmesh()
        if self.kpt_shift is None:
            self._set_centering()

        self.kpoints = Kpoints(comment=self.comment,
                               kpts=(self.kpt_mesh,),
                               kpts_shift=self.kpt_shift)
        if self.mode in (KpointsMode.primitive_uniform, KpointsMode.manual_set):
            self.num_kpts = \
                num_irreducible_kpoints(kpoints=self.kpoints,
                                        structure=self.corresponding_structure,
                                        symprec=self.symprec,
                                        angle_tolerance=self.angle_tolerance)
        else:
            irreducible_kpts = \
                irreducible_kpoints(structure=self.corresponding_structure,
                                    kpoints=self.kpoints,
                                    symprec=self.symprec,
                                    angle_tolerance=self.angle_tolerance)
            kpts = [k[0] for k in irreducible_kpts]
            weights = [k[1] for k in irreducible_kpts]
            labels = [None] * len(irreducible_kpts)
            self.kpoints = Kpoints(comment=self.comment,
                                   style=Kpoints.supported_modes.Reciprocal,
                                   num_kpts=len(kpts),
                                   kpts=kpts,
                                   kpts_weights=weights,
                                   labels=labels)
            self._add_band_kpts()
            self.num_kpts = len(self.kpoints.kpts)

    def _set_structure(self):

        if self.mode == KpointsMode.manual_set:
            self.corresponding_structure = self.initial_structure

        elif self.mode == KpointsMode.primitive_uniform:
            d = {"structure": self.initial_structure,
                 "symprec": self.symprec,
                 "angle_tolerance": self.angle_tolerance}

            self.corresponding_structure, _ = find_spglib_primitive(**d)
            sym_dataset = get_symmetry_dataset(**d)
            self.sg = sym_dataset["number"]
            self.sg_symbol = sym_dataset["international"]

        elif self.mode == KpointsMode.band:
            lattice = self.seekpath_info["primitive_lattice"]
            element_types = self.seekpath_info["primitive_types"]
            species = [Element.from_Z(i) for i in element_types]
            positions = self.seekpath_info["primitive_positions"]

            self.corresponding_structure = \
                Structure(lattice, species, positions)
            self.sg = self.seekpath_info["spacegroup_number"]
            self.sg_symbol = self.seekpath_info["spacegroup_international"]

        else:
            raise ValueError(f"{self.mode} mode is not supported yet.")

    def _set_kmesh(self):
        self.comment = f"Generated by vise. "
        if self.manual_kpts:
            self.comment += f"Mode: {self.mode}, manual kpt setting. "
            self.kpt_mesh = self.manual_kpts

        else:
            reciprocal_lat = \
                self.corresponding_structure.lattice.reciprocal_lattice
            reciprocal_abc = reciprocal_lat.abc
            if self.mode in (KpointsMode.band, KpointsMode.primitive_uniform):

                body_centered_ortho = {23, 24, 44, 45, 46, 71, 72, 73, 74}
                body_centered_tetra = {79, 80, 81, 82, 87, 88, 97, 98, 107, 108,
                                       109, 110, 119, 120, 121, 122, 139, 140,
                                       141, 142}
                if self.sg in body_centered_ortho | body_centered_tetra:
                    logger.warning("To keep the space group symmetry, the "
                                   "number of k-points along three directions "
                                   "are kept the same for oI and tI Bravais "
                                   "lattice.")
                    average_abc = pow(np.prod(reciprocal_lat.abc), 1 / 3)
                    reciprocal_abc = (average_abc, average_abc, average_abc)

            if self.only_even:
                self.kpt_mesh = [int(ceil(self.kpt_density * r / 2) * 2)
                                 for r in reciprocal_abc]
            else:
                self.kpt_mesh = [int(ceil(self.kpt_density * r))
                                 for r in reciprocal_abc]

            self.kpt_mesh = [i * self.factor for i in self.kpt_mesh]

            self.comment += f"Mode: {self.mode}, " \
                            f"kpt density: {self.kpt_density}, " \
                            f"factor: {self.factor}. "

    def _set_centering(self):
        self.kpt_shift = []
        if self.mode is KpointsMode.manual_set:
            angle = self.corresponding_structure.lattice.angles

            for i in range(3):
                # shift kpt mesh center only for the lattice vector being normal
                # to a lattice plane and even number of k-points.
                if max([abs(angle[i - 2] - 90), abs(angle[i - 1] - 90)]) < 1e-5:
                    self.kpt_shift.append(0.5)
                else:
                    self.kpt_shift.append(0.0)
        else:
            self.kpt_shift = kpt_centering[self.sg]

        # If number of k-point is odd, off centering has no meaning.
        for i in range(3):
            if self.kpt_mesh[i] % 2 == 1:
                self.kpt_shift[i] = 0.0

    def _add_band_kpts(self):
        k_path = self.seekpath_info["explicit_kpoints_rel"]
        k_path_label = self.seekpath_info["explicit_kpoints_labels"]

        formula = self.corresponding_structure.composition.reduced_formula
        self.kpoints.comment += \
            f"k-path added by seekpath. Formula: {formula} SG: {self.sg}"
        self.kpoints.num_kpts += len(k_path)
        self.kpoints.kpts += list(k_path)
        self.kpoints.labels += k_path_label
        self.kpoints.kpts_weights += [0] * len(k_path)

    @property
    def is_structure_changed(self):
        return not np.allclose(self.initial_structure.lattice.matrix,
                               self.corresponding_structure.lattice.matrix,
                               atol=self.symprec)


def irreducible_kpoints(structure: Structure,
                        kpoints: Kpoints,
                        symprec: float = SYMMETRY_TOLERANCE,
                        angle_tolerance: float = ANGLE_TOL) -> List[tuple]:
    """Estimate the irreducible k-points at the given k-mesh written in KPOINTS

    Args:
        kpoints (Kpoints):
            Kpoints file. Supported modes are only Gamma and Monkhorst.
        structure (Structure):
            Input structure corresponding to kpoints.
        symprec (float):
        angle_tolerance (float):
            See docstrings in spglib.

    Returns:
        List of irreducible kpoints and their weights tuples
        [(ir_kpoint, weight)], with ir_kpoint given in fractional coordinates.
    """
    if kpoints.style == Kpoints.supported_modes.Monkhorst:
        kpts_shift = [0.5, 0.5, 0.5]
    elif kpoints.style == Kpoints.supported_modes.Gamma:
        kpts_shift = [0.0, 0.0, 0.0]
    else:
        raise ValueError("Only Gamma or Monkhorst k-points modes supported.")

    kpts_shift = [kpts_shift[x] + kpoints.kpts_shift[x] for x in range(3)]
    # modf(x) returns (fraction part, integer part)
    shift = [1 if modf(i)[0] == 0.5 else 0 for i in kpts_shift]

    sga = SpacegroupAnalyzer(structure=structure,
                             symprec=symprec,
                             angle_tolerance=angle_tolerance)

    return sga.get_ir_reciprocal_mesh(mesh=kpoints.kpts[0], is_shift=shift)


def num_irreducible_kpoints(kpoints: Kpoints,
                            structure: Structure = None,
                            symprec: float = SYMMETRY_TOLERANCE,
                            angle_tolerance: float = ANGLE_TOL) -> int:
    """Calculate number of k-points

    Args:
        kpoints (Kpoints):
            Kpoints file. Supported modes are only Reciprocal, Gamma and
            Monkhorst.
        structure (Structure):
            Input structure corresponding to kpoints.
        symprec (float):
        angle_tolerance (float):
            See docstrings in spglib.

    Return:
        Int of Number of k-points
    """
    if kpoints.style == Kpoints.supported_modes.Reciprocal:
        return len(kpoints.kpts)
    else:
        return len(irreducible_kpoints(structure=structure,
                                       kpoints=kpoints,
                                       symprec=symprec,
                                       angle_tolerance=angle_tolerance))

