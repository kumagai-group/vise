# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from vise.analyzer.plot_band import Irrep
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
    # K does not exist in irreptables.
    actual = make_irreps_from_wavecar(special_point_symbols=["GM", "K", "H"],
                                      kpt_indices=[9, 15, 14],
                                      wavecar_filename=wavecar_file,
                                      poscar_filename=poscar_file,
                                      )
    assert actual.sg_num == 229

    assert_dataclass_almost_equal(actual.irreps["H"],
                                  Irrep([0.5, -0.5, 0.5],
                                        ["H1+", "H4-", "H3+"],
                                        [-5.93323, 2.71352, 3.80069],
                                        [1, 3, 2]), digit=3)
    assert_dataclass_almost_equal(actual.irreps["Γ"],
                                  Irrep([0.0, 0.0, 0.0],
                                        ["Γ1+", "Γ1+", "Γ4-", "Unknown"],
                                        [-6.5536, 0.9245, 6.7080, 8.6659],
                                        [1, 1, 3, 1]), digit=2)


def test_find_irrep():
    d = {"A": (0.991, 0.0), "B": (0.009, 0.0)}
    assert find_irrep(d) == "A"
    with pytest.raises(ViseNoIrrepError):
        find_irrep(d, threshold=0.992)

