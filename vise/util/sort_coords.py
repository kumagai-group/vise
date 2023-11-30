# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.
import numpy as np
from numpy.linalg import det
from numpy import concatenate, clip, dot, arctan2


def sort_coords(coords: np.ndarray) -> np.ndarray:
    """Sort coordinates based on the angle with first coord from the center.

    Args:
        coords (np.ndarray):
            Coordinates to be sorted. The format of coords is as follows.
            np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]

    Returns:
        np.ndarray for sorted coordinates.
    """
    if len(coords[0]) != 3:
        raise ValueError("Only valid for 3D vector")

    center = np.average(coords, axis=0)
    relative_coords = coords - center
    external_prod = np.cross(relative_coords[0], relative_coords[1])

    # Skip parallel vectors.
    if abs(np.linalg.norm(external_prod)) < 1e-8 and len(relative_coords) > 2:
        external_prod = np.cross(relative_coords[0], relative_coords[2])
    normal_to_12_plane = external_prod / np.linalg.norm(external_prod)

    v0 = relative_coords[0] / np.linalg.norm(relative_coords[0])

    def angle_between_v0(index: int) -> float:
        """
        Args:
            index (int): index of coords.

        Returns (float):
            Angle between rays from the center to rel_coords[0] and
            rel_coords[int].
        """
        v = relative_coords[index] / np.linalg.norm(relative_coords[index])
        matrix = concatenate(([v0], [v], [normal_to_12_plane]), axis=0)
        determinant = det(matrix)
        angle = arctan2(clip(dot(v0, v), -1.0, 1.0), determinant)
        return angle

    indices = [i for i in range(len(coords))]
    indices.sort(key=angle_between_v0)
    return coords[indices]
