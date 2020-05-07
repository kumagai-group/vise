# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.input_set.incar import ViseIncar, tabulated_string


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
                                 "ISYM": 0,
                                 "MAGMOM": [3, 3, 3]})
    actual = incar.get_string()
    expected = """# accuracy
PREC   =  Normal
LREAL  =  False
EDIFF  =  1.0

# symmetry
ISYM  =  0

# spin
MAGMOM  =  3 3 3"""
    assert actual == expected


def test_setting_to_str():
    incar = ViseIncar.from_dict({"EDIFF": 1.0,
                                 "MAGMOM": [3, 3, 3]})
    assert incar.setting_to_str("EDIFF") == str(1.0)
    assert incar.setting_to_str("MAGMOM") == "3 3 3"


def test_is_ncs_calc():
    incar = ViseIncar.from_dict({"LNONCOLLINEAR": "True"})
    assert incar.is_ncl_calc
    incar = ViseIncar.from_dict({})
    assert incar.is_ncl_calc is False


def test_tabulated_string():
    actual = tabulated_string([["short", "short"], ["long_str", "long_str"]])
    expected = """short     =  short
long_str  =  long_str"""
    assert actual == expected
