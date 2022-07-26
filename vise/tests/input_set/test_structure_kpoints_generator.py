# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure

from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.structure_kpoints_generator import \
    StructureKpointsGenerator
from vise.input_set.task import Task
from vise.util.structure_symmetrizer import StructureSymmetrizer


def test_constructor_num_kpt_multiplication_factor(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.dielectric_function,
                                          kpt_density=1.0)
    assert generator._num_kpt_factor == 3


def test_constructor_only_even_num_kpts(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          only_even_num_kpts=True,
                                          task=Task.cluster_opt,
                                          kpt_density=1.0)
    assert generator._only_even_num_kpts is False


def test_constructor_kpt_mode(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.cluster_opt,
                                          kpt_density=1.0)
    assert generator._kpt_mode is KpointsMode.uniform


def test_constructor_kpt_shift_structure_opt(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.structure_opt,
                                          kpt_density=1.0)
    assert generator._gamma_centered is None


def test_constructor_kpt_shift_dos(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.dos, kpt_density=1.0)
    assert generator._gamma_centered is True


def test_constructor_symmetrizer(sc_structure, mocker):
    mock = mocker.patch(
        "vise.input_set.structure_kpoints_generator.StructureSymmetrizer")
    StructureKpointsGenerator(sc_structure,
                              task=Task.dos,
                              symprec=10.0,
                              angle_tolerance=20.0,
                              band_ref_dist=30.0,
                              is_magnetization=True,
                              kpt_density=1.0)
    mock.assert_called_once_with(
        structure=sc_structure, symprec=10.0, angle_tolerance=20.0,
        band_mesh_distance=30.0, time_reversal=False)

# Structure related unittests.


def test_mc_structure_primitive(mc_structure):
    generator = StructureKpointsGenerator(mc_structure, task=Task.structure_opt,
                                          kpt_density=1.0)
    generator.generate_input()
    actual = generator.structure.lattice.matrix
    expected = [[ 6.,-3., 0.],
                [ 6., 3., 0.],
                [-4., 0., 7.]]
    np.testing.assert_array_almost_equal(actual, expected)


def test_mc_structure_band(mc_structure):
    generator = StructureKpointsGenerator(mc_structure, task=Task.band,
                                          kpt_density=1.0)
    generator.generate_input()
    actual = generator.structure.lattice.matrix
    expected = [[ 6., 3., 0.],
                [-6., 3., 0.],
                [-4., 0., 7.]]
    np.testing.assert_array_almost_equal(actual, expected)


def test_mc_structure_dos_uniform(mc_structure):
    generator = StructureKpointsGenerator(mc_structure, task=Task.dos,
                                          kpt_mode=KpointsMode.uniform,
                                          kpt_density=1.0)
    generator.generate_input()
    actual = generator.structure.lattice.matrix
    expected = [[ 6., 3., 0.],
                [ 6.,-3., 0.],
                [-4., 0., 7.]]
    np.testing.assert_array_almost_equal(actual, expected)


# KPOINTS related unittests.

def test_manual_kpts():
    lattice =[[10.0, 0.0, 0.0],
              [ 0.0, 8.0, 0.0],
              [-2.0, 0.0, 8.0]]
    coords = [[0.0, 0.0, 0.0]]
    # lattice angles are (90.0, 104.0362434679265, 90.0)
    # so only kpoint shift along b-direction must be 0.5.
    structure = Structure(lattice=lattice, species=["H"], coords=coords)
    generator = StructureKpointsGenerator(structure,
                                          kpt_density=2.5,
                                          task=Task.structure_opt,
                                          kpt_mode=KpointsMode.uniform)
    generator.generate_input()
    assert generator.kpoints.kpts[0] == [2, 2, 2]
    assert generator.kpoints.kpts_shift == [0.0, 0.5, 0.0]
    generator = StructureKpointsGenerator(structure,
                                          task=Task.structure_opt,
                                          kpt_mode=KpointsMode.uniform,
                                          gamma_centered=True,
                                          kpt_density=1.0)
    generator.generate_input()
    assert generator.kpoints.kpts_shift == [0.0, 0.0, 0.0]


def test_only_even_num_kpts(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.structure_opt,
                                          kpt_density=2.)
    generator.generate_input()
    num_kpt_list = generator.kpoints.kpts[0]
    assert num_kpt_list == [13, 13, 13]

    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.structure_opt,
                                          kpt_density=2.,
                                          only_even_num_kpts=True)
    generator.generate_input()
    num_kpt_list = generator.kpoints.kpts[0]
    assert num_kpt_list == [14, 14, 14]


def test_kpt_factor(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.structure_opt,
                                          kpt_density=2.,
                                          num_kpt_factor=2)
    generator.generate_input()
    num_kpt_list = generator.kpoints.kpts[0]
    assert num_kpt_list == [26, 26, 26]


def test_cluster_opt_kpoints(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.cluster_opt,
                                          kpt_density=2.)
    generator.generate_input()
    num_kpt_list = generator.kpoints.kpts[0]
    assert num_kpt_list == [1, 1, 1]
    kpt_shift = generator.kpoints.kpts_shift
    assert kpt_shift == [0, 0, 0]


def test_band_path(sc_structure):
    generator = StructureKpointsGenerator(sc_structure,
                                          task=Task.band,
                                          kpt_density=0.5,
                                          band_ref_dist=0.5)
    # reciprocal length = 2 * pi / 1 = 6.28
    # so num_kpt_list = [12, 12, 12] and kpt_shift = [0.5, 0.5, 0.5] are set.
    generator.generate_input()
    assert generator.kpoints.num_kpts == 44
    assert generator.kpoints.kpts[0] == [0.125, 0.125, 0.125]


def test_oi_ti_bravais():
    structure = Structure(Lattice.orthorhombic(20, 8, 6),
                          species=["H"]*2,
                          coords=[[0.0]*3, [0.5]*3])
    generator = StructureKpointsGenerator(structure,
                                          task=Task.structure_opt,
                                          kpt_density=5)
    generator.generate_input()

    primitive_structure = StructureSymmetrizer(structure).primitive

    expected_generator = StructureKpointsGenerator(primitive_structure,
                                                   task=Task.structure_opt,
                                                   kpt_density=5)
    expected_generator.generate_input()

    assert generator.num_kpts == expected_generator.num_kpts
    assert generator.kpoints.kpts_shift == expected_generator.kpoints.kpts_shift


def test_hexagonal():
    hexagonal_lattice =[[10.0, 0.0,      0.0],
                        [-5.0, 8.660254, 0.0],
                        [ 0.0, 0.0,     10.0]]
    coords = [[0.0, 0.0, 0.0]]
    # recipro_lat_abc = (0.85, 1.09, 1.31)
    structure = Structure(lattice=hexagonal_lattice,
                          species=["H"], coords=coords)
    generator = StructureKpointsGenerator(structure,
                                          task=Task.structure_opt,
                                          kpt_density=5)
    generator.generate_input()
    assert generator.kpoints.kpts_shift == [0.0, 0.0, 0.5]


def test_conventional_input(tmpdir):
    structure = Structure(Lattice.cubic(3.183372), species=["H"] * 2,
                          coords=[[0] * 3, [0.5] * 3])
    generator = StructureKpointsGenerator(initial_structure=structure,
                                          task=Task.structure_opt,
                                          kpt_density=2.5,
                                          num_kpt_factor=2)
    generator.generate_input()

    expected = Structure([[-1.591686, 1.591686, 1.591686],
                          [1.591686, -1.591686, 1.591686],
                          [1.591686, 1.591686, -1.591686]],
                         species=["H"], coords=[[0]*3])
    assert generator.structure == expected
    assert generator.num_kpts == 104
