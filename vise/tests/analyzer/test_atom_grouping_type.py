# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.analyzer.atom_grouping_type import AtomGroupingType


def test_group_by_atoms(complex_monoclinic_structure):
    class_type = AtomGroupingType.atoms
    target = ["0,1", "2,3,4"]
    group = class_type.grouped_atom_indices(structure=complex_monoclinic_structure,
                                            target=target)
    assert group == {"0,1": [0, 1], "2,3,4": [2, 3, 4]}


def test_group_by_elements(complex_monoclinic_structure):
    class_type = AtomGroupingType.elements
    target = ["H", "He"]
    group = class_type.grouped_atom_indices(structure=complex_monoclinic_structure,
                                            target=target)
    assert group == {"H": [0], "He": [1, 2, 3, 4]}


def test_group_by_elements_wo_target(complex_monoclinic_structure):
    class_type = AtomGroupingType.elements
    group = class_type.grouped_atom_indices(structure=complex_monoclinic_structure)
    assert group == {"H": [0], "He": [1, 2, 3, 4]}


def test_group_by_non_equiv_sites(complex_monoclinic_structure):
    class_type = AtomGroupingType.non_equiv_sites
    target = None
    group = class_type.grouped_atom_indices(structure=complex_monoclinic_structure,
                                            target=target)
    assert group == {"H1_a": [0], "He1_m": [1, 2], "He2_m": [3, 4]}


def test(complex_monoclinic_structure):
    print(complex_monoclinic_structure.composition.elements)
