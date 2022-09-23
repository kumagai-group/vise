# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from pathlib import Path

from vise.analyzer.atom_grouping_type import AtomGroupingType
from vise.analyzer.plot_diele_func_data import DieleFuncPlotType
from vise.cli.main import parse_args
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc

parent_dir = Path(__file__).parent


def test_structure_info_wo_options():
    parsed_args = parse_args(["si"])
    expected = Namespace(
        poscar="POSCAR",
        symprec=defaults.symmetry_length_tolerance,
        angle_tolerance=defaults.symmetry_angle_tolerance,
        show_conventional=False,
        show_primitive=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_structure_info_w_options():
    parsed_args = parse_args(["si",
                              "-p", "a",
                              "-s", "1",
                              "-a", "2",
                              "-c"])
    expected = Namespace(
        poscar="a",
        symprec=1.0,
        angle_tolerance=2.0,
        show_conventional=True,
        show_primitive=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_get_poscars_wo_options():
    parsed_args = parse_args(["gp", "-m", "mp-1234"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        poscar="POSCAR",
        prior_info=Path("prior_info.yaml"),
        mpid="mp-1234",
        func=parsed_args.func)
    assert parsed_args == expected


def test_get_poscars_w_options():
    parsed_args = parse_args(["get_poscar",
                              "-p", "a",
                              "-m", "123",
                              "-pi", "b.yaml"])

    expected = Namespace(
        poscar="a",
        prior_info=Path("b.yaml"),
        mpid="123",
        func=parsed_args.func)
    assert parsed_args == expected


def test_vasp_set_wo_options():
    parsed_args = parse_args(["vs"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        poscar=Path("POSCAR"),
        task=defaults.task,
        xc=defaults.xc,
        kpt_density=None,
        overridden_potcar=defaults.overridden_potcar,
        user_incar_settings=None,
        prev_dir=None,
        vasprun=defaults.vasprun,
        outcar=defaults.outcar,
        options=None,
        uniform_kpt_mode=False,
        file_transfer_type=None,
        func=parsed_args.func,
        )
    assert parsed_args == expected


def test_vasp_set_w_options():
    parsed_args = parse_args(["vs",
                              "--poscar", "POSCAR-tmp",
                              "-t", "band",
                              "-x", "pbesol",
                              "-k", "4.2",
                              "--potcar", "Mg_pv", "O_h",
                              "--user_incar_settings", "LREAD", "F",
                              "-d", "c",
                              "--vasprun", "vasprun_1",
                              "--outcar", "OUTCAR_1",
                              "--options", "encut", "800",
                              "--uniform_kpt_mode",
                              "--file_transfer", "WAVECAR", "C",
                              ])

    expected = Namespace(
        poscar=Path("POSCAR-tmp"),
        task=Task.band,
        xc=Xc.pbesol,
        kpt_density=4.2,
        overridden_potcar=["Mg_pv", "O_h"],
        user_incar_settings=["LREAD", "F"],
        prev_dir=Path("c"),
        vasprun=Path("vasprun_1"),
        outcar=Path("OUTCAR_1"),
        options=["encut", "800"],
        uniform_kpt_mode=True,
        file_transfer_type=["WAVECAR", "C"],
        func=parsed_args.func,
    )

    assert parsed_args == expected


def test_plot_band_wo_options():
    parsed_args = parse_args(["pb"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        vasprun=defaults.vasprun,
        kpoints_filename="KPOINTS",
        y_range=[-10.0, 10.0],
        filename="band.pdf",
        wavecar_filename=None,
        plotly=False,
        poscar="POSCAR",
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_plot_band_w_options():
    parsed_args = parse_args(["pb",
                              "--vasprun", "vasprun_1",
                              "--kpoints", "KPOINTS_1",
                              "--y_range", "-1.0", "1.0",
                              "--filename", "band_1.pdf",
                              "--wavecar", "WAVECAR",
                              "--plotly",
                              "--poscar", "POSCAR-test"
                              ])

    expected = Namespace(
        vasprun=Path("vasprun_1"),
        kpoints_filename="KPOINTS_1",
        y_range=[-1.0, 1.0],
        filename="band_1.pdf",
        wavecar_filename="WAVECAR",
        plotly=True,
        poscar="POSCAR-test",
        func=parsed_args.func,
    )

    assert parsed_args == expected


def test_plot_dos_wo_options():
    parsed_args = parse_args(["pd"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        vasprun=defaults.vasprun,
        outcar=defaults.outcar,
        type=AtomGroupingType.non_equiv_sites,
        legend=True,
        crop_first_value=True,
        x_range=None,
        y_max_ranges=None,
        target=None,
        filename="dos.pdf",
        base_energy=None,
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_plot_dos_w_options():
    parsed_args = parse_args(["pd",
                              "--vasprun", "vasprun_1",
                              "--outcar", "OUTCAR_1",
                              "-t", "atoms",
                              "-l", "False",
                              "-c", "False",
                              "--x_range", "-1.0", "1.0",
                              "-y", "-5.0", "5.0",
                              "--target", "1", "2",
                              "--filename", "dos_1.pdf",
                              "-b", "-1"
                              ])

    expected = Namespace(
        vasprun=Path("vasprun_1"),
        outcar=Path("OUTCAR_1"),
        type=AtomGroupingType.atoms,
        legend=False,
        crop_first_value=False,
        x_range=[-1.0, 1.0],
        y_max_ranges=[-5.0, 5.0],
        target=["1", "2"],
        filename="dos_1.pdf",
        base_energy=-1.0,
        func=parsed_args.func,
    )

    assert parsed_args == expected


def test_plot_absorption_wo_options():
    parsed_args = parse_args(["pdf"])
    expected = Namespace(
        vasprun=defaults.vasprun,
        outcar=defaults.outcar,
        input_csv_name=None,
        filename=None,
        directions=["ave"],
        y_range=None,
        calc_kk=False,
        ita=0.01,
        plot_type=DieleFuncPlotType.absorption_coeff,
        to_csv=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_plot_absorption_w_options():
    parsed_args = parse_args(["pdf",
                              "--vasprun", "vasprun_1",
                              "--outcar", "OUTCAR_1",
                              "--input_csv_name", "diele_func.csv",
                              "-f", "a",
                              "-d", "ave", "xx",
                              "-y", "-5", "5",
                              "-ckk",
                              "-i", "0.1",
                              "--plot_type", "diele_func",
                              "--to_csv"])
    expected = Namespace(
        vasprun=Path("vasprun_1"),
        outcar=Path("OUTCAR_1"),
        input_csv_name="diele_func.csv",
        filename="a",
        directions=["ave", "xx"],
        y_range=[-5.0, 5.0],
        calc_kk=True,
        ita=0.1,
        plot_type=DieleFuncPlotType.diele_func,
        to_csv=True,
        func=parsed_args.func)
    assert parsed_args == expected


def test_band_edge_wo_options():
    parsed_args = parse_args(["be"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        vasprun=defaults.vasprun,
        outcar=defaults.outcar,
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_band_edge_w_options():
    parsed_args = parse_args(["be",
                              "--vasprun", "vasprun_1",
                              "--outcar", "OUTCAR_1",
                              ])

    expected = Namespace(
        vasprun=Path("vasprun_1"),
        outcar=Path("OUTCAR_1"),
        func=parsed_args.func,
    )

    assert parsed_args == expected
