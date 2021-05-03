# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task
from vise.tests.helpers.assertion import assert_msonable


def test_task_str():
    assert "dos" == str(Task.dos)


def test_task_msonable():
    assert_msonable(Task.dos)


def test_task_from_string():
    assert Task.dos == Task.from_string("dos")


def test_task_from_string_raise_error():
    with pytest.raises(AttributeError):
        Task.from_string("not_exist_member")


def test_task_name_list():
    expected = 'structure_opt, structure_opt_rough, structure_opt_tight, ' \
               'cluster_opt, phonon_force, defect, band, dos, ' \
               'dielectric_dfpt, dielectric_finite_field, dielectric_function'
    assert Task.names_string() == expected


def test_is_atom_relaxed_lattice_fixed():
    assert Task.structure_opt_tight.is_atom_relaxed_lattice_fixed is False
    assert Task.cluster_opt.is_atom_relaxed_lattice_fixed is True
    assert Task.band.is_atom_relaxed_lattice_fixed is False


def test_is_lattice_relaxed():
    assert Task.structure_opt.is_lattice_relaxed is True
    assert Task.structure_opt_rough.is_lattice_relaxed is True
    assert Task.band.is_lattice_relaxed is False


def test_is_atom_relaxed():
    assert Task.structure_opt_tight.is_atom_relaxed is True
    assert Task.cluster_opt.is_atom_relaxed is True
    assert Task.band.is_atom_relaxed is False


def test_is_dielectric():
    assert Task.dielectric_function.is_dielectric is True
    assert Task.band.is_dielectric is False


def test_is_tight_calc():
    assert Task.structure_opt_tight.is_tight_calc is True
    assert Task.band.is_tight_calc is False


def test_is_plot_task():
    assert Task.band.is_plot_task is True
    assert Task.structure_opt.is_plot_task is False


def test_is_spectrum_task():
    assert Task.dos.is_spectrum_task is True
    assert Task.band.is_spectrum_task is False


def test_default_kpt_factor():
    assert Task.band.default_kpt_factor == 1
    assert Task.dos.default_kpt_factor == 2
    assert Task.dielectric_dfpt.default_kpt_factor == 2
    assert Task.dielectric_finite_field.default_kpt_factor == 2
    assert Task.dielectric_function.default_kpt_factor == 3


def test_default_kpt_mode():
    assert Task.band.default_kpt_mode == KpointsMode.band
    assert Task.defect.default_kpt_mode == KpointsMode.uniform
    assert Task.dos.default_kpt_mode == KpointsMode.primitive


def test_requisite_num_kpt_list():
    assert Task.cluster_opt.requisite_num_kpt_list == [1, 1, 1]
    assert Task.band.requisite_num_kpt_list is None


def test_requisite_only_even_num_kpts():
    assert Task.cluster_opt.requisite_only_even_num_kpts is False
    assert Task.band.requisite_only_even_num_kpts is None


def test_requisite_kpt_shift():
    assert Task.cluster_opt.requisite_gamma_centered is True
    assert Task.dos.requisite_gamma_centered is True
    assert Task.band.requisite_gamma_centered is None


def test_need_spin():
    assert Task.defect.need_spin is True
    assert Task.band.need_spin is False


def test_fine_to_inherit_structure_from_prev():
    assert Task.structure_opt.fine_to_inherit_structure_from_prev is True
    assert Task.phonon_force.fine_to_inherit_structure_from_prev is False

