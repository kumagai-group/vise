# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.task import Task


def test_task_str():
    assert "dos" == str(Task.dos)


def test_task_from_string():
    assert Task.dos == Task.from_string("dos")


def test_task_from_string_raise_error():
    with pytest.raises(AttributeError):
        Task.from_string("not_exist")


def test_task_kpt_factor():
    assert 1 == Task.band.kpt_factor
    assert 2 == Task.dos.kpt_factor
    assert 2 == Task.dielectric_dfpt.kpt_factor
    assert 2 == Task.dielectric_finite_field.kpt_factor
    assert 3 == Task.dielectric_function.kpt_factor


def test_task_name_list():
    expected = 'structure_opt, structure_opt_rough, structure_opt_tight, ' \
               'cluster_opt, phonon_force, defect, band, dos, ' \
               'dielectric_dfpt, dielectric_finite_field, dielectric_function'
    assert expected == Task.name_list()


def test_task_is_lattice_relax():
    assert Task.structure_opt.is_lattice_relax
    assert Task.structure_opt_rough.is_lattice_relax
    assert Task.band.is_lattice_relax is False


def test_task_is_plot_task():
    assert Task.band.is_plot_task
    assert Task.structure_opt.is_plot_task is False
