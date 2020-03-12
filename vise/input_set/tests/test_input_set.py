# -*- coding: utf-8 -*-

import tempfile
import os
from pathlib import Path

from vise.input_set.input_set import ViseInputSet
from vise.input_set.xc import Xc
from vise.util.testing import ViseTest


class InputSetTest(ViseTest):

    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO")
        self.input_set = ViseInputSet.make_input(structure=mgo, xc=Xc.hse)

    def test_write(self):
        with tempfile.TemporaryDirectory() as tmp_dirname:
            os.chdir(tmp_dirname)
            self.input_set.write_input(output_dir=tmp_dirname)
            os.remove("INCAR")
            os.remove("POSCAR")
            os.remove("POTCAR")
            os.remove("KPOINTS")
            os.remove("vise.json")
            Path.cwd()  # may be safer to go back to cwd

    def test_dict(self):
        expected = self.input_set.as_dict()
        actual = ViseInputSet.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)

    def test_msonable(self):
        self.assertMSONable(self.input_set)
