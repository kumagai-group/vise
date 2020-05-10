# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from vise.analyzer.dos_data import PDos, DosData, DosBySpinEnergy
from vise.analyzer.plot_dos import DosPlotter

energies = [-9, 0, 9]
total = np.array([[0, 5, 0], [0, 5, 0]])
orbitals = {"s":   np.array([[0, 0, 0], [1, 0, 0]], dtype=float),
            "px":  np.array([[0, 0, 1], [1, 0, 0]], dtype=float),
            "py":  np.array([[0, 0, 2], [1, 0, 0]], dtype=float),
            "pz":  np.array([[0, 0, 3], [1, 0, 0]], dtype=float),
            "dxy": np.array([[0, 1, 0], [1, 0, 0]], dtype=float),
            "dyz": np.array([[0, 2, 0], [1, 0, 0]], dtype=float),
            "dxz": np.array([[0, 3, 0], [1, 0, 0]], dtype=float),
            "dx2": np.array([[0, 4, 0], [1, 0, 0]], dtype=float),
            "dz2": np.array([[0, 5, 0], [1, 0, 0]], dtype=float),
            "f_3": np.array([[1, 0, 0], [1, 0, 0]], dtype=float),
            "f_2": np.array([[2, 0, 0], [1, 0, 0]], dtype=float),
            "f_1": np.array([[3, 0, 0], [1, 0, 0]], dtype=float),
            "f0":  np.array([[4, 0, 0], [1, 0, 0]], dtype=float),
            "f1":  np.array([[5, 0, 0], [1, 0, 0]], dtype=float),
            "f2":  np.array([[6, 0, 0], [1, 0, 0]], dtype=float),
            "f3":  np.array([[7, 0, 0], [1, 0, 0]], dtype=float)}
orbitals_2 = {k: v * 2 for k, v in orbitals.items()}
orbitals_3 = {k: v * 3 for k, v in orbitals.items()}

pdos_list = [PDos(**orbitals), PDos(**orbitals_2), PDos(**orbitals_3)]
reference_energy = 0.5


@pytest.fixture
def dos_data_list():
    dos_data_base_args = {"energies": energies,
                          "total": total,
                          "pdos": pdos_list}
    dos_data = DosData(**dos_data_base_args)
    dos_plot_data_w_lims = dos_data.dos_plot_data(
        grouped_atom_indices={"H": [0], "He": [1, 2]},
        xlim=[-100, 100],
        ylim_set=[[-20, 20], [-30, 30], [-40, 40]],
        base_energy=reference_energy,
        vertical_lines=[0.0, 1.0])
    dos_plot_data_wo_lims = dos_data.dos_plot_data(
        grouped_atom_indices={"H": [0], "He": [1, 2]},
        vertical_lines=[0.0, 1.0])

    return dos_data, dos_plot_data_w_lims, dos_plot_data_wo_lims


def test_pdos_s_p_d():
    p = sum([orbitals[orb] for orb in ["px", "py", "pz"]])
    d = sum([orbitals[orb] for orb in ["dxy", "dyz", "dxz", "dx2", "dz2"]])
    f = sum([orbitals[orb]
             for orb in ["f_3", "f_2", "f_1", "f0", "f1", "f2", "f3"]])
    pdos = PDos(**orbitals)
    assert_array_equal(pdos.p, p)
    assert_array_equal(pdos.d, d)
    assert_array_equal(pdos.f, f)


def test_pdos_wo_f_orbital():
    orbitals.pop("f_3")
    assert PDos(**orbitals).f is None


def test_pdos_add():
    pdos_sum = pdos_list[0] + pdos_list[1]
    assert_array_equal(pdos_sum.s, pdos_list[0].s + pdos_list[1].s)


def test_dos_data_energies(dos_data_list):
    _, dos_plot_data, _ = dos_data_list
    expected = [e - reference_energy for e in energies]
    assert dos_plot_data.relative_energies == expected


def test_dos_data_lim(dos_data_list):
    _, _, dos_plot_data_wo_lim = dos_data_list
    assert dos_plot_data_wo_lim.xlim == [-10, 10]
    assert dos_plot_data_wo_lim.ylim_set == [[-5.5, 5.5], [-154.0, 154.0], [-154.0, 154.0]]


def test_dos_data_manual_lim(dos_data_list):
    _, dos_plot_data_w_lim, _ = dos_data_list
    assert dos_plot_data_w_lim.xlim == [-100, 100]
    assert dos_plot_data_w_lim.ylim_set == [[-20, 20], [-30, 30], [-40, 40]]


def test_dos_data_total_dos(dos_data_list):
    _, dos_plot_data_w_lim, _ = dos_data_list
    assert dos_plot_data_w_lim.doses[0][0].name == "total"
    assert_array_equal(dos_plot_data_w_lim.doses[0][0].dos, total)


def test_dos_data_vertical_lines(dos_data_list):
    _, dos_plot_data_w_lim, _ = dos_data_list
    assert dos_plot_data_w_lim.vertical_lines == [-0.5, 0.5]


def test_dos_data_pdos_single(dos_data_list):
    _, dos_plot_data_w_lim, _ = dos_data_list
    pdos_sum = pdos_list[1] + pdos_list[2]

    assert_array_equal(dos_plot_data_w_lim.doses[1][0].name, "H-s")
    assert_array_equal(dos_plot_data_w_lim.doses[1][0].dos, pdos_list[0].s)
    assert_array_equal(dos_plot_data_w_lim.doses[1][1].name, "H-p")
    assert_array_equal(dos_plot_data_w_lim.doses[1][1].dos, pdos_list[0].p)

    assert_array_equal(dos_plot_data_w_lim.doses[2][0].name, "He-s")
    assert_array_equal(dos_plot_data_w_lim.doses[2][0].dos, pdos_sum.s)

    plotter = DosPlotter(dos_plot_data_w_lim)
    plotter.construct_plot()
    plotter.plt.show()


def test_orbital_dos():
    total_up = [0, 1]
    total_down = [2, 3]
    total_array = np.array([total_up, total_down])
    orbital_dos = DosBySpinEnergy("total", total_array)
    assert orbital_dos.max_dos() == max([max(total_up), max(total_down)])


def test_max_y_range():
    pdos = [PDos(s=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 px=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 py=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 pz=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 dxy=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 dyz=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 dxz=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 dx2=np.array([[0, 10, 20, 30, 40]], dtype=float),
                 dz2=np.array([[0, 10, 20, 30, 40]], dtype=float))]

    dos_data_base_args = {"energies": list(range(-2, 3)),
                          "total": np.array([[0, 1, 2, 3, 4]]),
                          "pdos": pdos}
    dos_data = DosData(**dos_data_base_args)

    dos_plot_data = dos_data.dos_plot_data(
        grouped_atom_indices={"H": [0]},
        vertical_lines=[0.0],
        xlim=[-1.1, 1.1])

    assert dos_plot_data.ylim_set == [[0, 3.3], [0, 165.0]]


def test_dos_by_spin_energy():
    dos = DosBySpinEnergy(name="test",
                          dos=np.array([[-3, -2, -1, 0],
                                        [ 1,  2,  3, 4]], dtype=float))
    mask = [True, False, False, True]
    assert dos.max_dos(mask) == 3.0


"""
TODO:
+ Allow to set vbm, cbm, efermi
+ Shift energy zero to vbm or efermi.

- + Shift energy zero to manual value
+ Create from FullDosInfo instance with total and atom- and orbital-decomposed pdos
+ Create class to generate grouped_atom_indices and their names

DONE:

"""
