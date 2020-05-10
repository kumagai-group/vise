# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.analyzer.atom_grouping_type import AtomGroupingType


def test_group_by_atoms(complex_ortho_structure):
    class_type = AtomGroupingType.atoms
    target = ["0,1", "2,3,4"]
    group = class_type.grouped_atom_indices(structure=complex_ortho_structure,
                                            target=target)
    assert group == {"0,1": [0, 1], "2,3,4": [2, 3, 4]}


def test_group_by_elements(complex_ortho_structure):
    class_type = AtomGroupingType.elements
    target = ["H", "He"]
    group = class_type.grouped_atom_indices(structure=complex_ortho_structure,
                                            target=target)
    assert group == {"H": [0], "He": [1, 2, 3, 4]}


def test_group_by_non_equiv_sites(complex_ortho_structure):
    class_type = AtomGroupingType.non_equiv_sites
    target = None
    group = class_type.grouped_atom_indices(structure=complex_ortho_structure,
                                            target=target)
    assert group == {"H_a1": [0], "He_m1": [1, 2], "He_m2": [3, 4]}


"""
TODO
- Return grouped atom indices from ["1,2,3", "4,5,6"] to {"1,2,3": [1,2,3}, "4,5,6": [4, 5, 6]}



"""