# -*- coding: utf-8 -*-
import numpy as np

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def normalized_random_3d_vector():
    """
    Generates a random 3d unit vector with a uniform spherical distribution.
    stackoverflow.com/questions/5408276/python-uniform-spherical-distribution

    Return:
        normed random vector
    """
    phi = np.random.uniform(0, np.pi * 2)
    cos_theta = np.random.uniform(-1, 1)
    theta = np.arccos(cos_theta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def random_vector(normed_vector, distance):
    """
    Returns a vector scaled by distance * x, where 0<x<1.

    Args:
        normed_vector (3x1 array): Normed 3d vector.
        distance (float): distance
    """
    return normed_vector * distance * np.random.random()