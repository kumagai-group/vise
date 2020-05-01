# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.util.enum import ExtendedEnum


class KpointsMode(ExtendedEnum):
    """K-point generation type

       Supporting modes are:
       "band":
           Kpoints with the band path will be returned based on the
           seekpath program. The space group is analyzed and primitive
           unitcell that must be used for the band structure calculation is
           returned as well.
       "primitive":
           Kpoints with uniform k-point sampling. The k-point sampling mesh
           and centering are determined based on the standardized primitive
           unitcell. Structure is also changed if not primitive.
       "uniform":
           Kpoints with uniform k-point sampling. The k-point sampling mesh
           and centering are determined based on the given lattice. Note
           that only when the angles are 90 degrees, the centering is
           shifted along the perpendicular direction.
           This mode is useful when calculating the supercells.
    """
    band = "band"
    primitive = "primitive"
    uniform = "uniform"

    @property
    def band_or_primitive(self):
        return True if self in [self.band, self.primitive] else False
