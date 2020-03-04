# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
from copy import deepcopy

from argparse import Namespace

from unittest.mock import patch

from pymatgen.core.structure import IStructure

from vise.util.testing import ViseTest
from vise.cli.main_function import (
    vasp_settings_from_args, get_poscar_from_mp, vasp_set, vasp_run,
    chempotdiag, plot_band, plot_dos, band_gap)
from vise.cli.tests.test_main import default_vasp_args, symprec_args
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.config import KPT_DENSITY


logger = get_logger(__name__)


vasp_args = default_vasp_args


class VaspSettingsFromArgsTest(ViseTest):
    def setUp(self) -> None:
        vasp = {"ldauu": ["Mn", "5", "O", "3"],
                "ldaul": ["Mn", "3", "O", "2"],
                "vise_opts": ["only_even", "True"],
                "user_incar_settings": ["POTIM", "0.4"],
                "additional_user_incar_settings": ["ALGO", "Fast"],
                "potcar_set": ["Mn_pv"],
                "potcar_set_name": "normal",
                "charge": 1}
        self.args = Namespace(**vasp)

    def test(self):
        user_incar_settings, vis_base_kwargs = \
            vasp_settings_from_args(self.args)

        expected = {"POTIM": 0.4, "ALGO": "Fast"}
        self.assertEqual(expected, user_incar_settings)

        expected = {"ldauu": {"Mn": 5, "O": 3},
                    "ldaul": {"Mn": 3, "O": 2},
                    "only_even": True,
                    "override_potcar_set": {"Mn": "Mn_pv"},
                    "potcar_set_name": "normal",
                    "charge": 1}
        self.assertEqual(expected, vis_base_kwargs)


class GetPoscarFromMpTest(ViseTest):
    def setUp(self) -> None:
        self.args_num = Namespace(number=110, poscar="POSCAR")
        self.kwargs = {"elements": ["Mg", "O"],
                       "e_above_hull": 1.5,
                       "molecules": True}
        self.args_elements = Namespace(**self.kwargs)
        self.args_none = Namespace()

    def test_number(self):
        get_poscar_from_mp(self.args_num)
        expected = IStructure.from_str("""Mg1
1.0
2.922478 0.000000 -1.033252
-1.461239 2.530940 -1.033252
0.000000 0.000000 3.099756
Mg
1
direct
0.000000 0.000000 0.000000 Mg""", fmt="poscar")
        actual = IStructure.from_file("POSCAR")
        self.assertEqual(expected, actual)

    @patch('vise.cli.main_function.make_poscars_from_mp')
    def test_mp_poscars(self, mock_make_poscars):
        get_poscar_from_mp(self.args_elements)
        mock_make_poscars.assert_called_with(**self.kwargs)

    def test_warning(self):
        expected = ['WARNING:vise.cli.main_function:Set mp number or elements']
        with self.assertLogs(level='WARNING') as cm:
            get_poscar_from_mp(self.args_none)
            self.assertEqual(expected, cm.output)

    def tearDown(self) -> None:
        if os.path.exists("POSCAR"):
            os.remove("POSCAR")


class VaspSetTest(ViseTest):
    def setUp(self) -> None:
        self.min_default_kwargs = {"poscar": "POSCAR",
                                   "kpt_density": KPT_DENSITY,
                                   "standardize_structure": True,
                                   "prior_info": True,
                                   "dirs": ["."],
                                   }
        self.args_print = Namespace(print=True,
                                    json="vasp_input_set.json",
                                    **self.min_default_kwargs,
                                    **vasp_args, **symprec_args)
        self.prev_dir = "a"
        self.args_prev = Namespace(print=False,
                                   json=None,
                                   prev_dir=self.prev_dir,
                                   **self.min_default_kwargs,
                                   **vasp_args, **symprec_args)

        self.args_normal = Namespace(print=False,
                                     json=None,
                                     **self.min_default_kwargs,
                                     prev_dir=None,
                                     **vasp_args, **symprec_args)

    @patch('vise.cli.main_function.ViseInputSet.load_json')
    def test_print(self, mock):
        vasp_set(self.args_print)
        mock.assert_called_with("vasp_input_set.json")

    @patch('vise.cli.main_function.ViseInputSet.from_prev_calc')
    def test_prev_dir(self, mock):
        vasp_set(self.args_prev)
        kwargs = {"override_potcar_set": vasp_args["potcar_set"] or {},
                  "potcar_set_name": vasp_args["potcar_set_name"],
                  "ldauu": vasp_args["ldauu"] or {},
                  "ldaul": vasp_args["ldaul"] or {},
                  "charge": vasp_args["charge"],
                  "symprec": symprec_args["symprec"],
                  "angle_tolerance": symprec_args["angle_tolerance"],
                  "kpt_density": self.min_default_kwargs["kpt_density"],
                  "standardize_structure":
                      self.min_default_kwargs["standardize_structure"],
                  "files_to_transfer":
                      {"CHGCAR": "C", "WAVECAR": "M", "WAVEDER": "M"},
                  "task": Task.from_string(vasp_args["task"]),
                  "xc": Xc.from_string(vasp_args["xc"])}

        mock.assert_called_with("a", **kwargs)

    @patch('vise.cli.main_function.Structure.from_file')
    @patch('vise.cli.main_function.ViseInputSet.make_input')
    def test_input_set(self, mock, mock_structure):
        vasp_set(self.args_normal)

        user_incar_settings = {}

        kwargs = {"structure": mock_structure(),
                  "override_potcar_set": vasp_args["potcar_set"] or {},
                  "potcar_set_name": vasp_args["potcar_set_name"],
                  "ldauu": vasp_args["ldauu"] or {},
                  "ldaul": vasp_args["ldaul"] or {},
                  "charge": vasp_args["charge"],
                  "symprec": symprec_args["symprec"],
                  "angle_tolerance": symprec_args["angle_tolerance"],
                  "kpt_density": self.min_default_kwargs["kpt_density"],
                  "standardize_structure":
                      self.min_default_kwargs["standardize_structure"],
                  "task": Task.from_string(vasp_args["task"]),
                  "xc": Xc.from_string(vasp_args["xc"]),
                  "user_incar_settings": user_incar_settings}

        mock.assert_called_with(**kwargs)


class VaspRunTest(ViseTest):
    def setUp(self) -> None:
        self.min_default_kwargs = {"poscar": "POSCAR",
                                   "kpt_density": KPT_DENSITY,
                                   "standardize_structure": True,
                                   "prior_info": True,
                                   }
        self.args_print = Namespace(json="vasp_input_set.json",
                                    **self.min_default_kwargs,
                                    **vasp_args, **symprec_args)


class TestChemPotDiag(ViseTest):
    def setUp(self) -> None:
        self.kwargs_1 = {
            "draw_phase_diagram": True,
            "vasp_dirs": None,
            "vasprun": None,
            "elements": ["Mg", "O"],
            "target_comp": "MgO",
            "filename": None,
            "parse_gas": True,
            "partial_pressures": ["O", "1e+5"],
            "temperature": 1000.0}

        self.from_mp_pd_filename = Namespace(**self.kwargs_1)

        self.kwargs_2 = {
            "draw_phase_diagram": False,
            "vasp_dirs": ["a_dir", "b_dir"],
            "vasprun": "vasprun.xml.finish",
            "elements": None,
            "target_comp": "MgO",
            "filename": "cpd.pdf",
            "parse_gas": False,
            "partial_pressures": ["O", "1e+5"],
            "temperature": 0.0}

        self.from_dir_cpd = Namespace(**self.kwargs_2)

    @patch('vise.cli.main_function.PDPlotter')
    @patch('vise.cli.main_function.PhaseDiagram')
    @patch('vise.cli.main_function.FreeEnergyEntrySet.from_mp')
    def test_from_mp_pd_filename(self, mock_entry_set, mock_pd, mock_pd_plot):
        chempotdiag(self.from_mp_pd_filename)
        mock_entry_set.assert_called_with(self.kwargs_1["elements"])
        mock_pd.assert_called_with(entries=mock_entry_set().entries)
        mock_pd_plot.assert_called_with(mock_pd())

    @patch('vise.cli.main_function.ChemPotDiag.from_phase_diagram')
    @patch('vise.cli.main_function.PhaseDiagram')
    @patch('vise.cli.main_function.FreeEnergyEntrySet.from_vasp_files')
    def test_from_dir_cpd(self, mock_entry_set, mock_pd, mock_cpd):
        chempotdiag(self.from_dir_cpd)
        mock_entry_set.assert_called_with(
            directory_paths=self.kwargs_2["vasp_dirs"],
            vasprun=self.kwargs_2["vasprun"],
            parse_gas=self.kwargs_2["parse_gas"],
            temperature=self.kwargs_2["temperature"],
            partial_pressures=self.kwargs_2["partial_pressures"])
        mock_pd.assert_called_with(entries=mock_entry_set().entries)
        mock_cpd.assert_called_with(mock_pd(),
                                    target_comp=self.kwargs_2["target_comp"])


class TestPlotBand(ViseTest):
    def setUp(self) -> None:
        self.kwargs = {"kpoints": "KPOINTS",
                       "vasprun": "1/vasprun.xml",
                       "vasprun2": "2/vasprun.xml",
                       "absolute": True,
                       "y_range": False,
                       "legend": False}
        self.args = Namespace(filename="band.pdf", **self.kwargs, **symprec_args)

    @patch('vise.cli.main_function.PrettyBSPlotter.from_vasp_files')
    def test_get_band_plot(self, mock):
        plot_band(self.args)
        kwargs = deepcopy(self.kwargs)
        kwargs["kpoints_filenames"] = kwargs.pop("kpoints")
        kwargs["vasprun_filenames"] = kwargs.pop("vasprun")
        kwargs["vasprun2_filenames"] = kwargs.pop("vasprun2")
        mock.assert_called_with(**kwargs, **symprec_args)


class TestGetDosPlotter(ViseTest):
    def setUp(self) -> None:
        self.kwargs = {"vasprun": "vasprun.xml",
                       "cbm_vbm": [3, 1],
                       "pdos_type": "element",
                       "specific": None,
                       "orbital": True,
                       "x_range": False,
                       "ymaxs": [2.0, 3.0],
                       "absolute": True,
                       "legend": False,
                       "crop_first_value": True}
        self.args = Namespace(filename="dos.pdf", **self.kwargs, **symprec_args)

    @patch('vise.cli.main_function.get_dos_plot')
    def test_get_dos_plot(self, mock):
        plot_dos(self.args)
        kwargs = deepcopy(self.kwargs)
        kwargs["zero_at_efermi"] = not kwargs.pop("absolute")
        kwargs["xlim"] = kwargs.pop("x_range")
        mock.assert_called_with(**kwargs, **symprec_args)


class TestBandGap(ViseTest):
    def setUp(self) -> None:
        self.kwargs = {"vasprun": "vasprun.xml", "outcar": "OUTCAR"}
        self.args = Namespace(**self.kwargs)

    @patch('vise.cli.main_function.band_gap_properties',
           return_value=({"a": 1}, {"b": 2}, {"c": 3}))
    def test_succeed(self, mock):
        band_gap(self.args)
        mock.assert_called_with(**self.kwargs)

