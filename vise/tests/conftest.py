# -*- coding: utf-8 -*-

from pathlib import Path

import pytest
from pymatgen.core.structure import Structure


@pytest.fixture(scope="session")
def test_data_files():
    return Path(__file__).parent / "test_data_files"


@pytest.fixture(scope="session")
def sc_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    results = Structure(lattice=lattice, species=["H"], coords=coords)

    return results


@pytest.fixture(scope="session")
def mc_structure():
    lattice =[[ 6.0, 3.0, 0.0],
              [ 6.0,-3.0, 0.0],
              [-4.0, 0.0, 7.0]]
    coords = [[0.0, 0.0, 0.0]]
    structure = Structure(lattice=lattice, species=["H"], coords=coords)
    return structure


@pytest.fixture(scope="session")
def complex_ortho_structure():
    lattice =[[10.0,  0.0,  0.0],
              [ 0.0, 20.0,  0.0],
              [-4.0,  0.0, 30.0]]
    coords = [[0.0, 0.0, 0.0],
              [0.1, 0.0, 0.0],
              [0.9, 0.0, 0.0],
              [0.2, 0.0, 0.0],
              [0.8, 0.0, 0.0],
              ]
    structure = Structure(lattice=lattice,
                          species=["H", "He", "He", "He", "He"],
                          coords=coords)
    return structure

