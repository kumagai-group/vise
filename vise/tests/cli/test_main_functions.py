# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os
from argparse import Namespace
from copy import deepcopy
from pathlib import Path

import pytest
from pymatgen.core import Structure, Element
from vise.analyzer.atom_grouping_type import AtomGroupingType
from vise.analyzer.dielectric_function import DieleFuncData
from vise.analyzer.plot_diele_func_data import DieleFuncPlotType
from vise.cli.main_functions import get_poscar_from_mp, VaspSet, plot_band, \
    plot_dos, band_edge_properties, plot_diele_func, \
    structure_info, get_most_stable_mp_id_from_formula, \
    get_most_stable_mp_id_from_formula_w_new_mprester
from vise.defaults import defaults
from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task
from vise.input_set.xc import Xc

default_option_args = {"poscar": Path("POSCAR"),
                       "task": Task.structure_opt,
                       "xc": Xc.pbe,
                       "kpt_density": 1.0,
                       "overridden_potcar": ["Mn_pv"]}

default_args = deepcopy(default_option_args)
default_args.update({"dirs": [Path.cwd()],
                     "user_incar_settings": None,
                     "prev_dir": None,
                     "options": None,
                     "file_transfer_type": None,
                     "uniform_kpt_mode": False,
                     "vasprun": Path("vasprun.xml"),
                     "outcar": Path("outcar.xml"),
                     })


def test_structure_info(mocker):
    args = Namespace(poscar=Path("POSCAR"), symprec=0.1, angle_tolerance=5,
                     show_conventional=False, show_primitive=False)
    lattice = [[10.0,  0.0,  0.0], [0.0, 10.0,  0.0], [-2.0,  0.0, 10.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]

    mock = mocker.patch("vise.cli.main_functions.Structure")
    mock.from_file.return_value = Structure(lattice=lattice, species=["H"] * 2,
                                            coords=coords)
    structure_info(args)
    args = Namespace(poscar=Path("POSCAR"), symprec=0.1, angle_tolerance=5,
                     show_conventional=True, show_primitive=False)
    structure_info(args)
    args = Namespace(poscar=Path("POSCAR"), symprec=0.1, angle_tolerance=5,
                     show_conventional=True, show_primitive=True)
    structure_info(args)


def test_get_most_stable_mp_id_from_formula():
    actual = get_most_stable_mp_id_from_formula_w_new_mprester("MgO")
    assert actual == "mp-1265"
#    actual = get_most_stable_mp_id_from_formula("MgO")
#    assert actual == "mp-1265"


def test_get_most_stable_mp_id_from_formula_not_present_formula():
    actual = get_most_stable_mp_id_from_formula_w_new_mprester("MgO100")
    assert actual is None


def test_get_poscar_from_mp(tmpdir):
    args = Namespace(mpid="mp-110", formula="He")
    tmpdir.chdir()
    get_poscar_from_mp(args)
    expected = """Mg1
1.0
   2.8949120000000002    0.0000000000000000   -1.0235060000000000
  -1.4474560000000001    2.5070679999999999   -1.0235060000000000
   0.0000000000000000    0.0000000000000000    3.0705179999999999
Mg
1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Mg"""
    assert Structure.from_file("POSCAR") == Structure.from_str(expected,
                                                               fmt="POSCAR")
    assert Path("prior_info.yaml").read_text() == """band_gap: 0.0
data_source: new MPRester mp-110
total_magnetization: 0.0
"""
#     expected = """Mg1
# 1.0
# -1.789645  1.789645  1.789645
#  1.789645 -1.789645  1.789645
#  1.789645  1.789645 -1.789645
# 1
# direct
# 0.000000 0.000000 0.000000 Mg
# """
#     assert Structure.from_file("POSCAR") == Structure.from_str(expected,
#                                                                fmt="POSCAR")
#     assert Path("prior_info.yaml").read_text() == """band_gap: 0.0
# data_source: mp-110
# icsd_ids:
# - 180455
# - 642652
# total_magnetization: 0.0001585
# """
    # Need to remove file to avoid the side effect for other unittests.
    os.remove("prior_info.yaml")


def test_get_poscar_from_mp_by_formula(tmpdir):
    args = Namespace(mpid=None, formula="Mg")
    print(tmpdir)
    tmpdir.chdir()
    get_poscar_from_mp(args)
    print(Structure.from_file("POSCAR").to(fmt="poscar"))
    expected = """Mg1
1.0
  -2.2442560000000000   -2.2442560000000000    0.0000000000000000
  -2.2442560000000000    0.0000000000000000   -2.2442560000000000
   0.0000000000000000   -2.2442560000000000   -2.2442560000000000
Mg
1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Mg"""
    assert Structure.from_file("POSCAR") == Structure.from_str(expected,
                                                               fmt="POSCAR")
    assert Path("prior_info.yaml").read_text() == """band_gap: 0.0
data_source: new MPRester mp-1056702
total_magnetization: 8.6e-06
"""
    # expected = """Mg3
# 1.0
# 7.698262 -1.605623 0.000000
# 7.698262 1.605623 0.000000
# 7.363377 0.000000 2.760785
# Mg
# 3
# direct
# 0.000000 0.000000 0.000000 Mg
# 0.222208 0.222208 0.222208 Mg
# 0.777792 0.777792 0.777792 Mg"""
#     assert Structure.from_file("POSCAR") == Structure.from_str(expected,
#                                                                fmt="POSCAR")
#     assert Path("prior_info.yaml").read_text() == """band_gap: 0.0
# data_source: mp-1094122
# icsd_ids: []
# total_magnetization: 0.00010333333333333333
# """
    # Need to remove file to avoid the side effect for other unittests.
    os.remove("prior_info.yaml")


def test_get_poscar_from_mp_wo_arguments():
    args = Namespace(mpid=None, formula=None)
    with pytest.raises(ValueError):
        get_poscar_from_mp(args)


test_data = [
    ({}, {}, {}, {}),

    ({"user_incar_settings": ["ISPIN", "2"]},
     {"ISPIN": 2},
     {},
     {}),

    ({"options": ["only_even_num_kpts", "True"]},
     {},
     {"only_even_num_kpts": True},
     {}),

    ({"uniform_kpt_mode": True},
     {},
     {"kpt_mode": KpointsMode.uniform},
     {}),

    ({"prev_dir": Path("a"), "file_transfer_type": ["file", "c"]},
     {},
     {},
     {"x": "y"}),

    ({"options": ["kpt_density", "3.5"], "kpt_density": 2.5},
     {},
     {"kpt_density": 2.5},
     {})
]


# Be careful about literal; all the words are connected (see position of commas)
@pytest.mark.parametrize("modified_settings,"
                         "overridden_incar_settings,"
                         "overridden_options_args,"
                         "prior_info", test_data)
def test_user_incar_settings(mocker,
                             modified_settings,
                             overridden_incar_settings,
                             overridden_options_args,
                             prior_info):
    args = deepcopy(default_args)
    args.update(modified_settings)

    structure = mocker.patch("vise.cli.main_functions.Structure")
    prior_info_mock = mocker.patch(
        "vise.cli.main_functions.prior_info_from_calc_dir")
    options = mocker.patch("vise.cli.main_functions.CategorizedInputOptions")
    vif = mocker.patch("vise.cli.main_functions.VaspInputFiles")

    prior_info_mock.return_value.input_options_kwargs = prior_info

    name_space = Namespace(**args)
    VaspSet(name_space)

    option_args = deepcopy(default_option_args)
    option_args.update(overridden_options_args)
    option_args.update(prior_info)
    option_args["overridden_potcar"] = {"Mn": "Mn_pv"}
    option_args.pop("poscar")
    option_args["structure"] = structure.from_file.return_value

    options.assert_called_once_with(**option_args)

    incar_settings = defaults.user_incar_settings
    incar_settings.update(overridden_incar_settings)

    vif.assert_called_once_with(options.return_value, incar_settings)


# def test_user_incar_settings_with_kpt_density_in_options():
#     args = deepcopy(default_args)
#     args.update(dict(kpt_density=2.5, options={"kpt_density": 3.5}))
#     name_space = Namespace(**args)
#     VaspSet(name_space)


def test_plot_band(tmpdir, test_data_files):
    tmpdir.chdir()
    args = Namespace(vasprun=test_data_files / "H_band_vasprun.xml",
                     kpoints_filename=str(test_data_files / "H_band_KPOINTS"),
                     y_range=[-10, 10],
                     wavecar_filename=None,
                     poscar=str(test_data_files / "H_band_POSCAR"),
                     plotly=False,
                     filename="test.pdf")

    plot_band(args)


def test_plot_band_irrep_and_plotly(test_data_files):
    args = Namespace(vasprun=test_data_files / "H_band_vasprun.xml",
                     kpoints_filename=str(test_data_files / "H_band_KPOINTS"),
                     y_range=[-10, 10],
                     wavecar_filename=str(test_data_files / "H_band_WAVECAR"),
                     poscar=str(test_data_files / "H_band_POSCAR"),
                     plotly=True,
                     filename="test.pdf")
    plot_band(args)


def test_plot_dos(tmpdir, test_data_files):
    tmpdir.chdir()  # comment out when one wants to see the figure
    args = Namespace(vasprun=test_data_files / "MnO_uniform_vasprun.xml",
                     outcar=test_data_files / "MnO_uniform_OUTCAR",
                     type=AtomGroupingType.non_equiv_sites,
                     legend=True,
                     base_energy=None,
                     crop_first_value=True,
                     x_range=[-5, 5],
                     y_max_ranges=[10, 5, 7],
                     target=["Mn", "O"],
                     filename="test.pdf",
                     title="test")
    plot_dos(args)


def test_band_edge_info(test_data_files):
    args = Namespace(vasprun=test_data_files / "MnO_uniform_vasprun.xml",
                     outcar=test_data_files / "MnO_uniform_OUTCAR")
    band_edge_properties(args)


def test_plot_diele_func(tmpdir, test_data_files):
    print(tmpdir)
    tmpdir.chdir()
    args = Namespace(vasprun=test_data_files / "MgSe_absorption_vasprun.xml",
                     outcar=test_data_files / "MgSe_absorption_OUTCAR",
                     input_csv_name=None,
                     y_range=[2.0, 8.0],
                     calc_kk=False,
                     ita=0.1,
                     filename=None,
                     directions=["ave"],
                     plot_type=DieleFuncPlotType.absorption_coeff,
                     to_csv=True,
                     title=None)
    plot_diele_func(args)

    args = Namespace(vasprun=test_data_files / "MgSe_absorption_vasprun.xml",
                     outcar=test_data_files / "MgSe_absorption_OUTCAR",
                     input_csv_name="diele_func_data.csv",
                     y_range=None,
                     calc_kk=False,
                     ita=0.1,
                     filename=None,
                     directions=["xx"],
                     plot_type=DieleFuncPlotType.refraction,
                     to_csv=False,
                     title="main func test")
    plot_diele_func(args)
