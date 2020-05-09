# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pymatgen.core.structure import Structure

from vise.util.enum import ExtendedEnum
from vise.util.structure_symmetrizer import StructureSymmetrizer


class AtomGroupingType(ExtendedEnum):
    non_equiv_sites = "non_equiv_sites"
    elements = "elements"
    atoms = "atoms"

    def grouped_atom_indices(self, structure, target):
        if self is self.atoms:
            grouped_indices = self._group_by_atoms(structure, target)
        elif self is self.elements:
            grouped_indices = self._group_by_elements(structure, target)
        else:
            grouped_indices = self._group_by_non_equiv_sites(structure, target)

        return grouped_indices

    @staticmethod
    def _group_by_atoms(structure, target):
        result = {}
        for num_str in target:
            num_list = [int(x) for x in num_str.split(",")]
            result[num_str] = num_list
        max_atom_index = max([max(val) for val in result.values()])
        if max_atom_index > len(structure):
            raise ValueError("Atom index is out of range.")
        return result

    @staticmethod
    def _group_by_elements(structure, target):
        result = {}
        for elem in target:
            result[elem] = \
                [i for i, site in enumerate(structure) if elem in site]
        return result

    @staticmethod
    def _group_by_non_equiv_sites(structure: Structure, target):
        symmetrizer = StructureSymmetrizer(structure)
        return symmetrizer.grouped_atom_indices()


