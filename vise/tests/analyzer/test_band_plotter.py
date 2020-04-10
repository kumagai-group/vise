# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
import unittest
import warnings
import pytest

from vise.util.testing import ViseTest
from vise.analyzer.band_plotter import labels_to_unicode

# from vise.analyzer.band_plotter import (
#     labels_to_unicode, PrettyBSPlotter, ModBSPlotter, make_bs_sym_line)


def test_labels_to_unicode():
    d = {"GAMMA": {"SIGMA": {"DELTA": "S_0"}}}
    assert {'Γ': {'Σ': {'Δ': '{\\rm S}_0'}}} == labels_to_unicode(d)


# @pytest.mark.filterwarnings(
#     'ignore::UserWarning',  # needed when POTCAR is absent.
#     'ignore::pymatgen.io.vasp.outputs.UnconvergedVASPWarning')
# def test_nhse_compared_band_display(test_data_files: Path):
#     kpoints = test_data_files / "MgO_band_KPOINTS"
#     vasprun1 = test_data_files / "MgO_band_vasprun.xml"
#     vasprun2 = test_data_files / "MgO_band_nhse_vasprun.xml"
#     plotter = PrettyBSPlotter.from_vasp_files(str(kpoints),
#                                               str(vasprun1),
#                                               str(vasprun2),
#                                               absolute=False)
#     plotter.show()
#
#
# @pytest.mark.filterwarnings('ignore::UserWarning')
# def test_ferromagnetic_band_display(test_data_files: Path):
#     ko2_band = make_bs_sym_line(
#         kpoints=str(test_data_files / "KO2_band_KPOINTS"),
#         vasprun=str(test_data_files / "KO2_band_vasprun.xml"))
#     plotter = PrettyBSPlotter(band=ko2_band, absolute=False)
#     plotter.show()
#
