# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
from numpy.testing import assert_array_equal
import pytest

from pymatgen import Structure

from vise.analyzer.dos_info import PDos, DosData
from vise.analyzer.plot_dos import DosPlotter

"""
TODO:
+ Create class to generate grouped_atom_indices and their names
+ Create DosPlotInfo with PDOS grouped by elements
+ Create DosPlotInfo with PDOS grouped by inequivalent sites
+ 

+ Shift energy zero to manual value
+ Allow to set vbm, cbm, efermi
+ Shift energy zero to vbm or efermi.
+ Create from FullDosInfo instance with total and atom- and orbital-decomposed pdos
+ Create DosInfo from vasp output.

DONE:
+ Create PDos class composed of orbital-decomposed pdos


"""


# @pytest.fixture
# def structure():
#     return Structure(lattice=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#                      species=["H"],
#                      coords=[[0, 0, 0]])


@pytest.fixture
def energies():
    return [-9, 0, 9]


@pytest.fixture
def total():
    return np.array([[0, 5, 0], [0, 5, 0]])


@pytest.fixture
def orbitals_f():
    return {"s":   np.array([[0, 0, 0], [1, 0, 0]], dtype=float),
            "px":  np.array([[0, 0, 1], [1, 0, 0]], dtype=float),
            "py":  np.array([[0, 0, 2], [1, 0, 0]], dtype=float),
            "pz":  np.array([[0, 0, 3], [1, 0, 0]], dtype=float),
            "dxy": np.array([[0, 1, 0], [1, 0, 0]], dtype=float),
            "dyz": np.array([[0, 2, 0], [1, 0, 0]], dtype=float),
            "dzx": np.array([[0, 3, 0], [1, 0, 0]], dtype=float),
            "dx2": np.array([[0, 4, 0], [1, 0, 0]], dtype=float),
            "dz2": np.array([[0, 5, 0], [1, 0, 0]], dtype=float),
            "f_3": np.array([[1, 0, 0], [1, 0, 0]], dtype=float),
            "f_2": np.array([[2, 0, 0], [1, 0, 0]], dtype=float),
            "f_1": np.array([[3, 0, 0], [1, 0, 0]], dtype=float),
            "f0":  np.array([[4, 0, 0], [1, 0, 0]], dtype=float),
            "f1":  np.array([[5, 0, 0], [1, 0, 0]], dtype=float),
            "f2":  np.array([[6, 0, 0], [1, 0, 0]], dtype=float),
            "f3":  np.array([[7, 0, 0], [1, 0, 0]], dtype=float)}


@pytest.fixture
def pdos_1st(orbitals_f):
    return PDos(**orbitals_f)


@pytest.fixture
def pdos_2nd(orbitals_f):
    orbitals_2 = {k: v * 2 for k, v in orbitals_f.items()}
    return PDos(**orbitals_2)


@pytest.fixture
def pdos_3rd(orbitals_f):
    orbitals_3 = {k: v * 3 for k, v in orbitals_f.items()}
    return PDos(**orbitals_3)


def test_pdos(pdos_1st, orbitals_f):
    assert_array_equal(pdos_1st.p, np.array([[0, 0, 6], [3, 0, 0]]))
    assert_array_equal(pdos_1st.d, np.array([[0, 15, 0], [5, 0, 0]]))
    assert_array_equal(pdos_1st.f, np.array([[28, 0, 0], [7, 0, 0]]))


def test_pdos_wo_f_orbital(orbitals_f):
    orbitals_f.pop("f_3")
    pdos_h_wo_f_3 = PDos(**orbitals_f)
    assert pdos_h_wo_f_3.f is None


def test_pdos_add(pdos_1st, pdos_2nd, pdos_3rd):
    pdos_sum = pdos_1st + pdos_2nd
    assert_array_equal(pdos_sum.s, pdos_3rd.s)


@pytest.fixture
def full_dos_data(energies, total, pdos_1st, pdos_2nd, pdos_3rd):
    dos_data = DosData(energies=energies,
                       total=total,
                       pdos=[pdos_1st, pdos_2nd, pdos_3rd],
                       grouped_atom_indices={"H": [0], "He": [1, 2]})
    return dos_data.dos_plot_data


@pytest.fixture
def full_dos_data_manual_lim(energies, total, pdos_1st, pdos_2nd, pdos_3rd):
    dos_data = DosData(energies=energies,
                       total=total,
                       pdos=[pdos_1st],
                       grouped_atom_indices={"H": [0]},
                       xlim=[-100, 100],
                       ylim_set=[[-20, 20], [-30, 30]])
    return dos_data.dos_plot_data


def test_full_dos_data_energies(full_dos_data, energies):
    assert full_dos_data.relative_energies == energies


def test_full_dos_data_lim(full_dos_data):
    assert full_dos_data.xlim == [-10, 10]
    assert full_dos_data.ylim_set == [[-5.5, 5.5], [-30.8, 30.8], [-154.0, 154.0]]


def test_full_dos_data_manual_lim(full_dos_data_manual_lim):
    assert full_dos_data_manual_lim.xlim == [-100, 100]
    assert full_dos_data_manual_lim.ylim_set == [[-20, 20], [-30, 30]]


def test_full_dos_data_total_dos(full_dos_data, total):
    assert full_dos_data.doses[0][0].name == "total"
    assert_array_equal(full_dos_data.doses[0][0].dos, total)


def test_full_dos_data_pdos(full_dos_data, pdos_1st, pdos_2nd, pdos_3rd):
    pdos_sum = pdos_2nd + pdos_3rd

    assert_array_equal(full_dos_data.doses[1][0].name, "H-s")
    assert_array_equal(full_dos_data.doses[1][0].dos, pdos_1st.s)
    assert_array_equal(full_dos_data.doses[1][1].name, "H-p")
    assert_array_equal(full_dos_data.doses[1][1].dos, pdos_1st.p)

    assert_array_equal(full_dos_data.doses[2][0].name, "He-s")
    assert_array_equal(full_dos_data.doses[2][0].dos, pdos_sum.s)

    plotter = DosPlotter(full_dos_data)
    plotter.construct_plot()
    plotter.plt.show()



