# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from pathlib import Path

from vise.analyzer.atom_grouping_type import AtomGroupingType
from vise.cli.main import parse_args
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc

parent_dir = Path(__file__).parent


def test_get_poscars_wo_options():
    parsed_args = parse_args(["gp"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        poscar="POSCAR",
        mpid=None,
        func=parsed_args.func)
    assert parsed_args == expected


def test_get_poscars_w_options():
    parsed_args = parse_args(["gp",
                              "-p", "a",
                              "-m", "123"])

    expected = Namespace(
        poscar="a",
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
        charge=0.0,
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
                              "-c", "10",
                              "--user_incar_settings", "LREAD", "F",
                              "-d", "c",
                              "--vasprun", "vasprun_1",
                              "--outcar", "OUTCAR_1",
                              "--options", "encut", "800",
                              "--uniform_kpt_mode",
                              "--file_transfer_type", "WAVECAR", "C",
                              ])

    expected = Namespace(
        poscar=Path("POSCAR-tmp"),
        task=Task.band,
        xc=Xc.pbesol,
        kpt_density=4.2,
        overridden_potcar=["Mg_pv", "O_h"],
        charge=10.0,
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
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_plot_band_w_options():
    parsed_args = parse_args(["pb",
                              "--vasprun", "vasprun_1",
                              "--kpoints", "KPOINTS_1",
                              "--y_range", "-1.0", "1.0",
                              "--filename", "band_1.pdf",
                              ])

    expected = Namespace(
        vasprun=Path("vasprun_1"),
        kpoints_filename="KPOINTS_1",
        y_range=[-1.0, 1.0],
        filename="band_1.pdf",
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
