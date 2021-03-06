# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Dict, List

from pymatgen.core.structure import Structure

from vise.util.enum import ExtendedEnum
from vise.util.structure_symmetrizer import StructureSymmetrizer


class AtomGroupingType(ExtendedEnum):
    non_equiv_sites = "non_equiv_sites"
    elements = "elements"
    atoms = "atoms"

    def grouped_atom_indices(self, structure, target=None):
        if self is self.atoms:
            grouped_indices = group_by_atoms(structure, target)
        elif self is self.elements:
            grouped_indices = group_by_elements(structure, target)
        else:
            grouped_indices = group_by_non_equiv_sites(structure)

        return grouped_indices


def group_by_atoms(structure: Structure, target):
    result = {}
    for num_str in target:
        num_list = [int(x) for x in num_str.split(",")]
        result[num_str] = num_list
    max_atom_index = max([max(val) for val in result.values()])
    if max_atom_index > len(structure):
        raise ValueError("Atom index is out of range.")
    return result


def group_by_elements(structure: Structure,
                      target: List[str] = None) -> Dict[str, List[int]]:
    result = {}
    target = target or [str(e) for e in structure.composition.elements]
    for elem in target:
        result[elem] = \
            [i for i, site in enumerate(structure) if elem in site]
    return result


def group_by_non_equiv_sites(structure: Structure):
    symmetrizer = StructureSymmetrizer(structure)
    return symmetrizer.grouped_atom_indices()


