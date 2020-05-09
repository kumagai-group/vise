# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import numpy as np
from numpy.testing import assert_array_equal

import pytest

from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen import Spin, Site, Orbital
from pymatgen.io.vasp import Vasprun

from vise.analyzer.vasp.dos_data import VaspDosData
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
def vasp_dos_data(mocker):
    stub_vasprun = mocker.Mock(spec=Vasprun)
    stub_vasprun.complete_dos = CompleteDos(None, total_dos, pdoses)
    return VaspDosData(stub_vasprun)


@pytest.fixture
def vasp_dos_data_crop_first_value(mocker):
    stub_vasprun = mocker.Mock(spec=Vasprun)
    stub_vasprun.complete_dos = CompleteDos(None, total_dos, pdoses)
    return VaspDosData(stub_vasprun, crop_first_value=True)


def test_energies(vasp_dos_data):
    assert vasp_dos_data.energies == energies


def test_total_dos(vasp_dos_data):
    assert_array_equal(vasp_dos_data.total, np.array([tdos]))


def test_pdos(vasp_dos_data):
    assert_array_equal(vasp_dos_data.pdos[0].s, np.array([pdos]))


def test_dos_for_crop_first_value(vasp_dos_data_crop_first_value):
    assert_array_equal(vasp_dos_data_crop_first_value.energies, energies[1:])
    assert_array_equal(vasp_dos_data_crop_first_value.total,
                       np.array([tdos[1:]]))
    assert_array_equal(vasp_dos_data_crop_first_value.pdos[0].s,
                       np.array([pdos[1:]]))


def test_actual_vasp_files(test_data_files: Path):
    vasprun_file = str(test_data_files / "MgO_dos_vasprun.xml")
    vasprun = Vasprun(vasprun_file)
    vasp_dos_data = VaspDosData(vasprun)
    plot_data = vasp_dos_data.dos_plot_data({"Mg": [0], "O": [1]},
                                            base_energy=1,
                                            vertical_lines=[0.0])
    plotter = DosPlotter(plot_data)
    plotter.construct_plot()
    plotter.plt.show()

