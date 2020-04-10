# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.input_set.incar import incar_flags, ViseIncar


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
    assert actual == expected


def test_from_string():
    incar_string = """
    ENCUT = 500; LREAD = True
     ALGO = D     ; PREC = A  
      MAGMOM = 5 5 5
    """
    actual = ViseIncar.from_string(incar_string)
    expected = {'ALGO': 'D', 'ENCUT': 500, 'LREAD': True,
                'MAGMOM': [5, 5, 5], 'PREC': 'A'}
    assert actual == expected


def test_get_string():
    incar = ViseIncar.from_dict({"PREC": "Normal",
                                 "LREAL": False,
                                 "EDIFF": 1.0,
                                 "ISYM": 0})
    actual = incar.get_string()
    expected = """# accuracy
PREC   =  Normal
LREAL  =  False
EDIFF  =  1.0

# symmetry
ISYM  =  0"""
    assert actual == expected


def test_is_ncs_calc():
    incar = ViseIncar.from_dict({"LNONCOLLINEAR": "True"})
    assert incar.is_ncl_calc
    incar = ViseIncar.from_dict({})
    assert incar.is_ncl_calc is False
