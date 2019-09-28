# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Callable, Union

from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest


class ViseTest(PymatgenTest):
    """ Extends PymatgenTest with some path modification. """
    MODULE_DIR = Path(__file__).absolute().parent
    TEST_FILES_DIR = MODULE_DIR / ".." / "test_files"
    POSCARS_DIR = TEST_FILES_DIR / "poscars"
#    CORE_DIR = TEST_FILES_DIR

    @classmethod
    def get_structure_by_name(cls, name: str):
        filename = cls.POSCARS_DIR / ("POSCAR-" + name)
        return Structure.from_file(filename)

    @classmethod
    def get_structure_by_sg(cls, sg: int):
        filename = cls.POSCARS_DIR / "poscar_by_sg" / ("POSCAR_" + str(sg))
        return Structure.from_file(filename)

    @classmethod
    def get_object_by_name(cls, method: Callable, name: Union[list, str]):
        """ return cls by passing classmethod returning cls """
        if isinstance(name, str):
            name = [name]

        filename = cls.TEST_FILES_DIR
        for n in name:
            filename = filename / n

        return method(filename)

    @classmethod
    def get_filename(cls, name):
        return cls.TEST_FILES_DIR / name
