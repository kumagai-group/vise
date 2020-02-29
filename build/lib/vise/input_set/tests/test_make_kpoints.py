# -*- coding: utf-8 -*-

import numpy as np

from pymatgen.io.vasp import Kpoints
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from vise.input_set.make_kpoints import (
    KpointsMode, MakeKpoints, irreducible_kpoints, num_irreducible_kpoints)
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class KpointsModeTest(ViseTest):
    def test(self):
        band = KpointsMode.from_string("band")
        self.assertEqual(KpointsMode.band, band)

    def test2(self):
        with self.assertRaises(AttributeError):
            KpointsMode.from_string("fail_string")


class MakeKpointsTest(ViseTest):

    def setUp(self) -> None:
        self.sg1 = self.get_structure_by_sg(1)  # test for hinuma reduced cell
        self.sg23 = self.get_structure_by_sg(23)  # I222, body-centered ortho
#        self.sg38 = self.get_structure_by_sg(38)  # Amm2, base-centered ortho

    def test_band_sg1(self):
        k = MakeKpoints(mode="band",
                        only_even=True,
                        structure=self.sg1)
        k.make_kpoints()
        expected_kpt_mesh = [2, 4, 4]
        expected_lattice = Lattice([[1.702429, 1.411278, - 11.284670],
                                    [1.502067, -6.074800, 0.000000],
                                    [-6.169335, 0.000000, 0.000000]])
        expected_shift = [0.5, 0.5, 0.5]
        self.assertEqual(expected_kpt_mesh, k.kpt_mesh)
        self.assertEqual(expected_lattice, k.corresponding_structure.lattice)
        self.assertEqual(expected_shift, k.kpt_shift)
        self.assertEqual(161, k.num_kpts)

    def test_band_sg23(self):
        k = MakeKpoints(mode="band",
                        only_even=True,
                        structure=self.sg23)
        k.make_kpoints()
        expected_kpt_mesh = [4, 4, 4]
        expected_lattice = Lattice([
            [-3.5435465044428671, 3.7004010644681178, 4.1207036549338687],
            [3.5435465044428671, -3.7004010644681178, 4.1207036549338687],
            [3.5435465044428671, 3.7004010644681178, -4.1207036549338687]])
        expected_shift = [0.5, 0.5, 0.5]
        self.assertEqual(expected_kpt_mesh, k.kpt_mesh)
        self.assertEqual(expected_lattice, k.corresponding_structure.lattice)
        self.assertEqual(expected_shift, k.kpt_shift)
        self.assertEqual(224, k.num_kpts)

    def test_primitive_uniform_sg1(self):
        k = MakeKpoints(mode="primitive_uniform",
                        structure=self.sg1)
        k.make_kpoints()
        expected_kpt_mesh = [3, 3, 2]
        expected_lattice = Lattice([[6.169335, 0.000000, 0.000000],
                                    [-1.502067, 6.074800, 0.000000],
                                    [-1.702429, -1.411278, 11.284670]])

        expected_shift = [0.0, 0.0, 0.5]
        self.assertEqual(expected_kpt_mesh, k.kpt_mesh)
        self.assertEqual(expected_lattice, k.corresponding_structure.lattice)
        self.assertEqual(expected_shift, k.kpt_shift)
        self.assertEqual(9, k.num_kpts)

    def test_manual_set_sg1(self):
        k = MakeKpoints(mode="manual_set",
                        structure=self.sg1)
        k.make_kpoints()
        expected_kpt_mesh = [3, 3, 2]
        expected_lattice = Lattice([[6.16933470, 0.00090247, -0.00277192],
                                    [-1.5029601, 6.07456979, -0.01071425],
                                    [-1.6971538, -1.3903666, 11.28805959]])

        expected_shift = [0.0, 0.0, 0.0]
        self.assertEqual(expected_kpt_mesh, k.kpt_mesh)
        self.assertEqual(expected_lattice, k.corresponding_structure.lattice)
        self.assertEqual(expected_shift, k.kpt_shift)
        self.assertEqual(10, k.num_kpts)


class IrreducibleKpointsTest(ViseTest):
    def setUp(self):
        lattice = [[10, 0, 0], [-5, 8.660254, 0], [0, 0, 5]]
        self.hexagonal = Structure(lattice=lattice,
                                   species=["H"],
                                   coords=[[0, 0, 0]])
        self.kpoints_g = Kpoints.gamma_automatic(kpts=(2, 2, 2))

        lattice2 = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        self.ortho = Structure(lattice=lattice2,
                               species=["H"],
                               coords=[[0, 0, 0]])
        self.kpoints_mp = Kpoints.monkhorst_automatic(kpts=(2, 2, 2))

        self.kpoints_reciprocal = \
            Kpoints(style=Kpoints.supported_modes.Reciprocal,
                    num_kpts=2,
                    kpts=((0, 0, 0), (0.5, 0.5, 0.5)),
                    kpts_weights=[1, 1])

    def test_gamma(self):
        actual = irreducible_kpoints(structure=self.hexagonal,
                                     kpoints=self.kpoints_g)
        expected = [(np.array([0.0,  0.0,  0.0]), 1),
                    (np.array([0.5,  0.0,  0.0]), 3),
                    (np.array([0.0,  0.0,  0.5]), 1),
                    (np.array([0.5,  0.0,  0.5]), 3)]
        for i in range(4):
            self.assertEqual(expected[i][0].tolist(), actual[i][0].tolist())
            self.assertEqual(expected[i][1], actual[i][1])

    def test_mp(self):
        actual = irreducible_kpoints(structure=self.ortho,
                                     kpoints=self.kpoints_mp)
        expected = [(np.array([0.25,  0.25,  0.25]), 8)]
        self.assertEqual(expected[0][0].tolist(), actual[0][0].tolist())
        self.assertEqual(expected[0][1], actual[0][1])

    def test_num_kpts(self):
        self.assertEqual(2, num_irreducible_kpoints(self.kpoints_reciprocal))
        self.assertEqual(4, num_irreducible_kpoints(structure=self.hexagonal,
                                                    kpoints=self.kpoints_g))
        self.assertEqual(1, num_irreducible_kpoints(structure=self.ortho,
                                                    kpoints=self.kpoints_mp))
