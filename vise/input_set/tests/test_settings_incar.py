# -*- coding: utf-8 -*-
from copy import deepcopy

from pymatgen.io.vasp import Potcar

from vise.input_set.settings_incar import TaskIncarSettings
from vise.input_set.task import Task

from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class TaskIncarSettingsTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO")
        self.default_kwargs = {"task": Task.structure_opt,
                               "structure": mgo,
                               "potcar": Potcar(["Mg", "O"]),
                               "num_kpoints": 8,
                               "max_enmax": 400.0,
                               "is_magnetization": False,
                               "vbm_cbm": [3.0, 8.0],
                               "npar_kpar": True,
                               "num_nodes": 1,
                               "encut": None,
                               "structure_opt_encut_factor": 1.3,
                               }

    def test_structure_opt(self):
        setting = TaskIncarSettings.from_options(**self.default_kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 3,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-07,
                    'NSW': 50,
                    'EDIFFG': -0.005,
                    'KPAR': 4,
                    'ENCUT': 520.0}
        self.assertEqual(expected, setting.settings)

    def test_structure_opt_rough(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.structure_opt_rough
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 3,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-04,
                    'NSW': 50,
                    'EDIFFG': -0.2,
                    'KPAR': 4,
                    'POTIM': 0.1,
                    'ENCUT': 520.0}
        self.assertEqual(expected, setting.settings)

    def test_structure_opt_tight(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.structure_opt_tight
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'ADDGRID': True,
                    'IBRION': 2,
                    'PREC': 'A',
                    'ISIF': 3,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-08,
                    'NSW': 50,
                    'EDIFFG': -0.001,
                    'KPAR': 4,
                    'ENCUT': 520.0}
        self.assertEqual(expected, setting.settings)

    def test_cluster_opt(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.cluster_opt
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 2,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-07,
                    'NSW': 50,
                    'EDIFFG': -0.005,
                    'KPAR': 4,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_phonon_force(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.phonon_force
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'ADDGRID': True,
                    'IBRION': 2,
                    'PREC': 'A',
                    'ISIF': 2,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-08,
                    'KPAR': 4,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_defect(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.defect
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 2,
                    'ISMEAR': 0,
                    'ISPIN': 2,
                    'LREAL': "A",
                    'EDIFF': 1e-05,
                    'NSW': 50,
                    'EDIFFG': -0.05,
                    'KPAR': 4,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_band(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.band
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-05,
                    'KPAR': 4,
                    'NBANDS': 12,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_dos(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.dos
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': -5,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-05,
                    'KPAR': 4,
                    'NBANDS': 12,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_dielectric_dfpt(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.dielectric_dfpt
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 8,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LREAL': False,
                    'LEPSILON': True,
                    'EDIFF': 1e-06,
                    'KPAR': 4,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_dielectric_finite_field(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.dielectric_finite_field
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 6,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': 0,
                    'ISPIN': 1,
                    'LCALCEPS': True,
                    'POTIM': 0.015,
                    'LREAL': False,
                    'EDIFF': 1e-06,
                    'KPAR': 4,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_dielectric_function(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.dielectric_function
        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': -5,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-05,
                    'KPAR': 4,
                    'NBANDS': 12,
                    'CSHIFT': 0.01,
                    'LOPTICS': True,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

    def test_args(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.dos
        kwargs["num_kpoints"] = 1
        kwargs["is_magnetization"] = True
        kwargs["vbm_cbm"] = None
        kwargs["encut"] = 800.0

        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': 0,
                    'ISPIN': 2,
                    'LREAL': False,
                    'EDIFF': 1e-05,
                    'KPAR': 1,
                    'NBANDS': 12,
                    'ENCUT': 800.0}
        self.assertEqual(expected, setting.settings)


    def test_args(self):
        kwargs = deepcopy(self.default_kwargs)
        kwargs["task"] = Task.dos
        kwargs["num_kpoints"] = 18
        kwargs["num_nodes"] = 2
        kwargs["structure_opt_encut_factor"] = 1.5

        setting = TaskIncarSettings.from_options(**kwargs)
        expected = {'IBRION': 2,
                    'PREC': 'N',
                    'ISIF': 0,
                    'ISMEAR': -5,
                    'ISPIN': 1,
                    'LREAL': False,
                    'EDIFF': 1e-05,
                    'KPAR': 6,
                    'NBANDS': 12,
                    'ENCUT': 400.0}
        self.assertEqual(expected, setting.settings)

