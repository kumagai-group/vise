# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.input_set.settings_structure_kpoints import TaskStructureKpoints
from vise.input_set.task import Task


class TaskStructureKpointsTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO")
        self.structure_opt = \
            TaskStructureKpoints.from_options(task=Task.band,
                                              original_structure=mgo,
                                              standardize_structure=True,
                                              sort_structure=True,
                                              is_magnetization=False,
                                              kpt_mode="primitive",
                                              kpt_density=3,
                                              kpt_shift=[0, 0, 0],
                                              only_even=True,
                                              band_ref_dist=0.025,
                                              factor=1,
                                              symprec=0.01,
                                              angle_tolerance=3)

    def test(self):
        expected = 287
        actual = self.structure_opt.kpoints.num_kpts
        self.assertEqual(expected, actual)

        expected = ['GAMMA', 'X', 'U', 'K', 'GAMMA', 'L', 'W', 'X']
        actual = [l for l in self.structure_opt.kpoints.labels if l]
        self.assertEqual(expected, actual)
