# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from vise.analyzer.irrep import Irrep, Character
from vise.analyzer.vasp.make_irreps import make_irreps_from_wavecar, find_irrep, \
    ViseNoIrrepError, special_points_from_kpoints
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


def test_special_points_from_kpoints(test_data_files: Path):
    kpoints_file = str(test_data_files / "H_band_KPOINTS")
    actual = special_points_from_kpoints(kpoints_file)
    assert actual == (["GM", "H", "N", "P"], [9, 14, 17, 24])


def test_make_irreps_with_actual_vasp_files(test_data_files: Path):
    wavecar_file = str(test_data_files / "H_band_WAVECAR")
    poscar_file = str(test_data_files / "H_band_POSCAR")
    actual = make_irreps_from_wavecar(special_point_characters=["GM", "H"],
                                      kpt_indices=[9, 14],
                                      wavecar_filename=wavecar_file,
                                      poscar_filename=poscar_file,
                                      )
    assert actual.sg_num == 229

    assert_dataclass_almost_equal(actual.irreps["H"],
                                  Irrep((0.5, -0.5, 0.5),
                                      [Character("H1+", -5.93323, 1),
                                       Character("H4-", 2.71352, 3),
                                       Character("H3+", 3.80069, 2)]), digit=5)

    assert_dataclass_almost_equal(actual.irreps["GM"],
                                  Irrep((0.0, 0.0, 0.0),
                                        [Character("GM1+", -6.5536, 1),
                                         Character("GM1+", 0.9245, 1),
                                         Character("GM4-", 6.7080, 3),
                                         Character("Unknown", 8.6659, 1)]),
                                  digit=2)


def test_find_irrep():
    d = {"A": (0.991, 0.0), "B": (0.009, 0.0)}
    assert find_irrep(d) == "A"
    with pytest.raises(ViseNoIrrepError):
        find_irrep(d, threshold=0.992)

