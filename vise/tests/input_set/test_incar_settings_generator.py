# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen.core import Composition
from pymatgen.io.vasp import Potcar

from vise.input_set.incar_settings_generator import (
    IncarSettingsGenerator)
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.defaults import defaults

@pytest.fixture()
def default_dict():
    return {"composition": Composition("HO2"),
            "symbol_list": ["H", "O"],
            "num_kpts": 5,
            "num_kpt_multiplication_factor": 1,
            "potcar": Potcar(["H", "O"]),
            "xc": Xc.pbe,
            "task": Task.structure_opt}


def test_pbe_structure_opt(default_dict):
    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    expected = {
        "ALGO": "Normal",
        "PREC": "Normal",
        "LREAL": False,
        "EDIFF": 1e-07,
        "ENCUT": 520.0,
        "LASPH": True,
        "NELM": 100,
        "ISIF": 3,
        "IBRION": 2,
        "EDIFFG": -0.005,
        "NSW": 50,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "LWAVE": False,
        "LCHARG": False,
        "LORBIT": 10,
        "KPAR": 1,
    }
    assert generator.incar_settings == expected


def test_scan_structure_opt(default_dict):
    default_dict["xc"] = Xc.scan
    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    assert generator.incar_settings["METAGGA"] == "SCAN"


def test_hse_structure_opt(default_dict):
    default_dict.update({"xc": Xc.hse,
                         "composition": Composition("UO2"),
                         "symbol_list": ["U", "O"],
                         "potcar": Potcar(["U", "O"]),
                         "exchange_ratio": 0.5,
                         "set_hubbard_u": True})

    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    expected = {
        "ALGO": "Damped",
        "PREC": "Normal",
        "LREAL": False,
        "EDIFF": 1e-07,
        "ENCUT": 520.0,
        "LASPH": True,
        "NELM": 100,
        "ISIF": 3,
        "IBRION": 2,
        "EDIFFG": -0.005,
        "NSW": 50,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "LWAVE": True,
        "LCHARG": False,
        "LORBIT": 10,
        "LHFCALC": True,
        "PRECFOCK": "Fast",
        "TIME": 0.4,
        "AEXX": 0.5,
        "HFSCREEN": 0.208,
        "LDAU": True,
        "LDAUTYPE": 2,
        "LMAXMIX": 6,
        "LDAUPRINT": 1,
        "LDAUU": [5, 0],
        "LDAUL": [3, -1],
        "KPAR": 1,
    }
    assert generator.incar_settings == expected


def test_dos(default_dict):
    default_dict.update({"vbm_cbm": [0.0, 1.0], "task": Task.dos})
    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    expected = {
        "ALGO": "Normal",
        "PREC": "Normal",
        "LREAL": False,
        "EDIFF": 1e-05,
        "ENCUT": 400.0,
        "LASPH": True,
        "NELM": 100,
        "ISIF": 0,
        "IBRION": 2,
        "NSW": 1,
        "ISMEAR": -5,
        "SIGMA": 0.1,
        "LWAVE": False,
        "LCHARG": False,
        "NBANDS": 19,
        "LORBIT": 11,
        "EMIN": -15.01,
        "EMAX": 16,
        "NEDOS": 3102,
        "KPAR": 1,
    }
    assert generator.incar_settings == expected


def test_defect(default_dict):
    default_dict.update({"task": Task.defect})
    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    expected = {
        "ALGO": "Normal",
        "PREC": "Normal",
        "LREAL": "Auto",
        "EDIFF": 1e-05,
        "ENCUT": 400.0,
        "LASPH": True,
        "NELM": 100,
        "ISIF": 2,
        "IBRION": 2,
        "EDIFFG": -0.03,
        "NSW": 50,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "ISPIN": 2,
        "LWAVE": False,
        "LCHARG": False,
        "LORBIT": 10,
        "KPAR": 1,
    }
    assert generator.incar_settings == expected


def test_phonon_force(default_dict):
    default_dict.update({"task": Task.phonon_force})
    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    expected = {
        "ADDGRID": True,
        "ALGO": "Normal",
        "PREC": "Accurate",
        "LREAL": False,
        "EDIFF": 1e-08,
        "ENCUT": 400.0,
        "LASPH": True,
        "NELM": 100,
        "ISIF": 2,
        "IBRION": 2,
        "NSW": 1,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "LWAVE": False,
        "LCHARG": False,
        "KPAR": 1,
    }
    assert generator.incar_settings == expected


def test_dielectric_function(default_dict):
    default_dict.update({"task": Task.dielectric_function})
    generator = IncarSettingsGenerator(**default_dict)
    # for k, v in generator.incar_settings.items():
    #     print(f'\"{k}\": {v}, ')
    expected = {
        "ALGO": "Normal",
        "PREC": "Normal",
        "LREAL": False,
        "EDIFF": 1e-05,
        "ENCUT": 400.0,
        "LASPH": True,
        "NELM": 100,
        "ISIF": 0,
        "IBRION": 2,
        "NSW": 1,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "LWAVE": False,
        "LCHARG": False,
        "NBANDS": 19,
        "LORBIT": 10,
        "EMIN": -20.01,
        "EMAX": 20,
        "NEDOS": 4002,
        "LOPTICS": True,
        "CSHIFT": 0.0,
        "KPAR": 1,
    }
    assert generator.incar_settings == expected


def test_with_band_gap(default_dict):
    default_dict.update({"band_gap": defaults.band_gap_criterion + 1e-5})
    generator = IncarSettingsGenerator(**default_dict)
    assert generator.incar_settings["ISMEAR"] == -5


def test_with_band_gap_band(default_dict):
    default_dict.update({"band_gap": defaults.band_gap_criterion + 1e-5,
                         "task": Task.band})
    generator = IncarSettingsGenerator(**default_dict)
    assert generator.incar_settings["ISMEAR"] == 0


def test_ldau_option(default_dict):
    generator = IncarSettingsGenerator(composition=Composition("Zn"),
                                       symbol_list=["Zn"],
                                       potcar=Potcar(["Zn"]),
                                       num_kpts=5,
                                       num_kpt_multiplication_factor=1,
                                       xc=Xc.pbe,
                                       task=Task.structure_opt)
    assert "LDAU" not in generator.incar_settings

