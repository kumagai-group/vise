# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.input_set.incar import incar_flags


def test_incar_flags():
    expected = {'ADDGRID': 'task_optional',
                'EDIFF': 'task_required',
                'ENAUG': None,
                'ENCUT': 'task_optional',
                'LASPH': 'common_required',
                'LREAL': 'task_required',
                'NELM': 'common_required',
                'NELMDL': None,
                'NELMIN': None,
                'PREC': 'task_required'}
    actual = incar_flags["accuracy"]
    assert expected == actual



