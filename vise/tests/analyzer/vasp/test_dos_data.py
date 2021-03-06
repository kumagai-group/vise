# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import numpy as np
from numpy.testing import assert_array_equal

import pytest

from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.core import Site
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.io.vasp import Vasprun

from vise.analyzer.vasp.dos_data import DosDataFromVasp
from vise.analyzer.plot_dos import DosPlotter


efermi = 0.5
energies = [-1, 0, 1]
tdos = [3, 4, 5]
pdos = [1, 2, 3]

total_dos = Dos(efermi=efermi, energies=energies,
                densities={Spin.up: tdos})
pdoses = {Site("H", [0, 0, 0]): {Orbital.s: {Spin.up: pdos},
                                 Orbital.px: {Spin.up: pdos},
                                 Orbital.py: {Spin.up: pdos},
                                 Orbital.pz: {Spin.up: pdos},
                                 Orbital.dxy: {Spin.up: pdos},
                                 Orbital.dyz: {Spin.up: pdos},
                                 Orbital.dxz: {Spin.up: pdos},
                                 Orbital.dz2: {Spin.up: pdos},
                                 Orbital.dx2: {Spin.up: pdos}}}


@pytest.fixture
def dos_data_from_vasp(mocker):
    stub_vasprun = mocker.Mock(spec=Vasprun)
    stub_vasprun.complete_dos = CompleteDos(None, total_dos, pdoses)
    return DosDataFromVasp(stub_vasprun).make_dos_data()


@pytest.fixture
def vasp_dos_data_crop_first_value(mocker):
    stub_vasprun = mocker.Mock(spec=Vasprun)
    stub_vasprun.complete_dos = CompleteDos(None, total_dos, pdoses)
    return DosDataFromVasp(stub_vasprun, crop_first_value=True).make_dos_data()


@pytest.fixture
def vasp_dos_data_energy_window(mocker):
    stub_vasprun = mocker.Mock(spec=Vasprun)
    stub_vasprun.complete_dos = CompleteDos(None, total_dos, pdoses)
    return DosDataFromVasp(stub_vasprun,
                           crop_first_value=True,
                           energy_window=[0.1, 1.1]).make_dos_data()


def test_energies(dos_data_from_vasp):
    assert dos_data_from_vasp.energies == energies


def test_total_dos(dos_data_from_vasp):
    assert_array_equal(dos_data_from_vasp.total, np.array([tdos]))


def test_pdos(dos_data_from_vasp):
    assert_array_equal(dos_data_from_vasp.pdos[0].s, np.array([pdos]))


def test_dos_for_crop_first_value(vasp_dos_data_crop_first_value):
    assert_array_equal(vasp_dos_data_crop_first_value.energies, energies[1:])
    assert_array_equal(vasp_dos_data_crop_first_value.total,
                       np.array([tdos[1:]]))
    assert_array_equal(vasp_dos_data_crop_first_value.pdos[0].s,
                       np.array([pdos[1:]]))


def test_energy_window(vasp_dos_data_energy_window):
    assert_array_equal(vasp_dos_data_energy_window.energies, energies[2])
    assert_array_equal(vasp_dos_data_energy_window.total,
                       np.array([tdos[2:]]))
    assert_array_equal(vasp_dos_data_energy_window.pdos[0].s,
                       np.array([pdos[2:]]))


def test_actual_vasp_files(test_data_files: Path):
    vasprun_file = str(test_data_files / "MgO_dos_vasprun.xml")
    vasprun = Vasprun(vasprun_file)
    dos_data = DosDataFromVasp(vasprun,
                               base_energy=1,
                               vertical_lines=[0.0]).make_dos_data()
    plot_data = dos_data.dos_plot_data({"Mg": [0], "O": [1]})
    plotter = DosPlotter(plot_data)
    plotter.construct_plot()
    plotter.plt.show()

