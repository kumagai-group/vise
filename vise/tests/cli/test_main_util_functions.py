# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pymatgen.core import Element
from vise.cli.main_util_functions import make_atom_poscars, \
    make_spin_decomposed_volumetric_files, make_light_weight_vol_data

_file_path = "vise.cli.main_util_functions"


# def test_make_phonon_poscars():
#     args = Namespace(unitcell=)


def test_make_atom_poscars(mocker):
    args = Namespace(dirname=Path("a"), elements=[Element.H, Element.He])
    mock = mocker.patch(f"{_file_path}.make_atom_poscar_dirs")
    make_atom_poscars(args)
    mock.assert_called_once_with(Path("a"), [Element.H, Element.He])


def test_make_spin_decomposed_volumetric_files(mocker):
    args = Namespace(chgcar="CHGCAR")
    mock_chgcar = mocker.patch(f"{_file_path}.Chgcar")
    mock_make_spin_charges = mocker.patch(f"{_file_path}.make_spin_charges")
    chgcar = mock_chgcar.from_file.return_value
    make_spin_decomposed_volumetric_files(args)
    mock_chgcar.from_file.assert_called_once_with("CHGCAR")
    mock_make_spin_charges.assert_called_once_with(chgcar)


def test_make_light_weight_vol_data(mocker, tmpdir):
    tmpdir.chdir()
    Path("original.vesta").write_text("a")

    args = Namespace(volumetric_file="CHGCAR",
                     output_lw_volumetric_filename=Path("CHGCAR_lw"),
                     border_fractions=[0.1],
                     output_vesta_filename=Path("test_w_chg.vesta"),
                     original_vesta_file=Path("original.vesta"))
    mock_chgcar = mocker.patch(f"{_file_path}.Chgcar")
    mock_vol_text = mocker.patch(f"{_file_path}.light_weight_vol_text")
    mock_calc_isurfs = mocker.patch(f"{_file_path}.calc_isurfs")
    mock_add_density = mocker.patch(f"{_file_path}.add_density")

    chgcar = mock_chgcar.from_file.return_value
    mock_add_density.return_value = "add_density_return"
    mock_vol_text.return_value = "light_weight_vol_text"
    mock_calc_isurfs.return_value = [0.0001]

    make_light_weight_vol_data(args)
    mock_vol_text.assert_called_once_with(chgcar, [0.1])
    mock_calc_isurfs.assert_called_once_with([0.1], True,
                                             chgcar.structure.volume)
    mock_add_density.assert_called_once_with("a", [0.0001], "CHGCAR_lw")
    assert Path("test_w_chg.vesta").read_text() == "add_density_return"
    assert Path("CHGCAR_lw").read_text() == "light_weight_vol_text"
