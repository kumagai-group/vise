# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task


def test_task_str():
    assert "dos" == str(Task.dos)


def test_task_from_string():
    assert Task.dos == Task.from_string("dos")


def test_task_from_string_raise_error():
    with pytest.raises(AttributeError):
        Task.from_string("not_exist_member")


def test_task_name_list():
    expected = 'structure_opt, structure_opt_rough, structure_opt_tight, ' \
               'cluster_opt, phonon_force, defect, band, dos, ' \
               'dielectric_dfpt, dielectric_finite_field, dielectric_function'
    assert Task.name_list() == expected


def test_task_is_lattice_relax():
    assert Task.structure_opt.is_lattice_relaxed_task is True
    assert Task.structure_opt_rough.is_lattice_relaxed_task is True
    assert Task.band.is_lattice_relaxed_task is False


def test_is_tight_calc():
    assert Task.structure_opt_tight.is_tight_calc is True
    assert Task.band.is_tight_calc is False


def test_is_atom_relaxed_task():
    assert Task.structure_opt_tight.is_atom_relaxed_task is True
    assert Task.cluster_opt.is_atom_relaxed_task is True
    assert Task.band.is_atom_relaxed_task is False


def test_is_dielectric_task():
    assert Task.dielectric_function.is_dielectric_task is True
    assert Task.band.is_dielectric_task is False


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


def test_incar_isif():
    assert Task.structure_opt_tight.incar_isif == 3
    assert Task.defect.incar_isif == 2
    assert Task.band.incar_isif == 0
    assert Task.dielectric_function.incar_isif == 0


def test_incar_ediff():
    assert Task.structure_opt_tight.incar_ediff == 1e-8
    assert Task.structure_opt.incar_ediff == 1e-7
    assert Task.dielectric_finite_field.incar_ediff == 1e-6
    assert Task.dielectric_function.incar_ediff == 1e-5
    assert Task.structure_opt_rough.incar_ediff == 1e-4


def test_incar_ediffg():
    assert Task.structure_opt_tight.incar_ediffg_optional == -0.001
    assert Task.structure_opt.incar_ediffg_optional == -0.005
    assert Task.defect.incar_ediffg_optional == -0.03
    assert Task.structure_opt_rough.incar_ediffg_optional == -0.2
    assert Task.dos.incar_ediffg_optional is None


def test_incar_ibrion():
    assert Task.dielectric_dfpt.incar_ibrion == 8
    assert Task.band.incar_ibrion == 2


def test_incar_lreal():
    assert Task.defect.incar_lreal == "Auto"
    assert Task.band.incar_lreal is False


def test_incar_ispin():
    assert Task.defect.incar_ispin == 2
    assert Task.band.incar_ispin == 1


def test_incar_prec():
    assert Task.structure_opt_tight.incar_prec == "Accurate"
    assert Task.structure_opt.incar_prec == "Normal"


def test_incar_nsw():
    assert Task.structure_opt_tight.incar_nsw == 50
    assert Task.defect.incar_nsw == 50
    assert Task.band.incar_nsw == 0


def test_incar_potim_optional():
    assert Task.structure_opt_rough.incar_potim_optional == 0.1
    assert Task.structure_opt.incar_potim_optional is None


def test_incar_addgrid_optional():
    assert Task.structure_opt_tight.incar_addgrid_optional == 0.1
    assert Task.structure_opt.incar_addgrid_optional is None
