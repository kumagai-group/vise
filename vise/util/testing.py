# -*- coding: utf-8 -*-
from pathlib import Path

from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest


class ViseTest(PymatgenTest):
    """Extends PymatgenTest with some path modification. """
    MODULE_DIR = Path(__file__).absolute().parent
    TEST_FILES_DIR = MODULE_DIR / ".." / "test_files"
    POSCARS_DIR = TEST_FILES_DIR / "poscars"

    @classmethod
    def get_structure_by_name(cls, name: str) -> Structure:
        filename = cls.POSCARS_DIR / ("POSCAR-" + name)
        return Structure.from_file(filename)

    @classmethod
    def get_structure_by_sg(cls, sg: int) -> Structure:
        filename = cls.POSCARS_DIR / "poscar_by_sg" / ("POSCAR_" + str(sg))
        return Structure.from_file(filename)

    @classmethod
    def get_filename(cls, name) -> Path:
        return cls.TEST_FILES_DIR / name
