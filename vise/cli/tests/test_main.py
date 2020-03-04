# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
import os
import shutil
from pathlib import Path

from vise.util.testing import ViseTest
from vise.custodian_extension.jobs import ViseVaspJob
from vise.cli.main import (
    simple_override, parse_args)
from vise.cli.main_tools import get_default_args
from vise.util.mp_tools import make_poscars_from_mp
from vise.input_set.input_set import ViseInputSet
from vise.config import (
    SYMMETRY_TOLERANCE, ANGLE_TOL, KPT_DENSITY, KPT_INIT_DENSITY, TIMEOUT)


parent_dir = Path(__file__).parent


class SimpleOverrideTest(ViseTest):
    def setUp(self) -> None:
        vise_yaml = parent_dir / ".." / ".." / "test_files" / "vise.yaml"
        shutil.copy(str(vise_yaml), "vise.yaml")

    def test(self):
        d = {"xc": "pbesol", "symprec": 1}
        simple_override(d, ["xc"])
        expected = {'xc': 'hse', 'symprec': 1}

        self.assertEqual(expected, d)

    def tearDown(self) -> None:
        os.remove("vise.yaml")


vis_defaults = get_default_args(ViseInputSet.make_input)
vis_defaults.update(ViseInputSet.TASK_OPTIONS)
vis_defaults.update(ViseInputSet.XC_OPTIONS)

kpt_conv_defaults = get_default_args(ViseVaspJob.kpt_converge)

# func is a pointer so need to point the same address.
# The directly passed args from argpaser.
default_vasp_args = {
    "potcar_set": None,
    "potcar_set_name": vis_defaults["potcar_set_name"],
    "xc": vis_defaults["xc"],
    "task": vis_defaults["task"],
    "vise_opts": None,
    "user_incar_settings": None,
    "additional_user_incar_settings": None,
    "ldauu": vis_defaults["ldauu"],
    "ldaul": vis_defaults["ldaul"],
    "charge": 0.0}

symprec_args = {"symprec": SYMMETRY_TOLERANCE,  "angle_tolerance": ANGLE_TOL}


class MainGetPoscarsTest(ViseTest):

    def test_get_poscars_wo_options(self):
        actual = parse_args(["gp"])
        print(actual)
        d = get_default_args(make_poscars_from_mp)
        # func is a pointer so need to point the same address.
        expected = Namespace(
            poscar=None,
            number=None,
            elements=None,
            e_above_hull=d["e_above_hull"],
            molecules=d["molecules"],
            func=actual.func)
#        self.assertEqual(expected, actual)

    def test_get_poscars_w_options(self):
        parsed_args = parse_args(["gp",
                                  "-p", "POSCAR",
                                  "-n", "123",
                                  "-e", "Mg", "O",
                                  "--e_above_hull", "0.1",
                                  "--molecules", "F"])

        expected = Namespace(
            poscar="POSCAR",
            number=123,
            elements=["Mg", "O"],
            e_above_hull=0.1,
            molecules=False,
            func=parsed_args.func)
        self.assertEqual(expected, parsed_args)


class MainVaspSetTest(ViseTest):
    def test_vasp_set_wo_options(self):
        parsed_args = parse_args(["vs"])
        # func is a pointer so need to point the same address.
        expected = Namespace(
            print=False,
            json="vise.json",
            poscar="POSCAR",
            kpt_density=KPT_DENSITY,
            standardize_structure=True,
            prior_info=True,
            dirs=["."],
            prev_dir=None,
            func=parsed_args.func,
            **default_vasp_args, **symprec_args)
        self.assertEqual(expected, parsed_args)

    def test_vasp_set_yaml(self):
        vise_yaml = parent_dir / ".." / ".." / "test_files" / "vise.yaml"
        shutil.copy(str(vise_yaml), str("vise.yaml"))
        parsed_args = parse_args(["vs",
                                  "--ldauu", "Mg", "5",
                                  "-d", "c"])
        os.remove("vise.yaml")
        expected = Namespace(
            print=False,
            potcar_set=default_vasp_args["potcar_set"],
            potcar_set_name="gw",
            xc="hse",
            task=default_vasp_args["task"],
            vise_opts=default_vasp_args["vise_opts"],
            user_incar_settings=default_vasp_args["user_incar_settings"],
            additional_user_incar_settings=
            default_vasp_args["additional_user_incar_settings"],
            ldauu=["Mg", "5"],
            ldaul=["Mn", "2", "O", "1"],
            charge=default_vasp_args["charge"],
            symprec=2,
            angle_tolerance=symprec_args["angle_tolerance"],
            json="vise.json",
            poscar="POSCAR",
            kpt_density=KPT_DENSITY,
            standardize_structure=True,
            prior_info=True,
            dirs=["."],
            prev_dir="c",
            func=parsed_args.func)

        self.assertEqual(expected, parsed_args)

    # Here, test vasp_parser and symprec_parser as well.
    def test_vasp_set_w_options(self):
        parsed_args = parse_args(["vs",
                                  "--print",
                                  "--potcar", "Mg_pv", "O_h",
                                  "--potcar_set_name", "gw",
                                  "-x", "pbesol",
                                  "-t", "band",
                                  "--vise_opts", "encut", "800",
                                  "--user_incar_settings", "LREAD", "F",
                                  "-auis", "ALGO", "D",
                                  "--ldauu", "Mg", "5",
                                  "--ldaul", "Mg", "1",
                                  "-c", "10",
                                  "--symprec", "0.0001",
                                  "--angle_tolerance", "10",
                                  "--json", "tmp.json",
                                  "--poscar", "POSCAR-tmp",
                                  "-k", "4.2",
                                  "-s", "F",
                                  "-pi", "T",
                                  "--dirs", "a", "b",
                                  "-d", "c"])

        expected = Namespace(
            print=True,
            potcar_set=["Mg_pv", "O_h"],
            potcar_set_name="gw",
            xc="pbesol",
            task="band",
            vise_opts=["encut", "800"],
            user_incar_settings=["LREAD", "F"],
            additional_user_incar_settings=["ALGO", "D"],
            ldauu=["Mg", "5"],
            ldaul=["Mg", "1"],
            charge=10.0,
            symprec=0.0001,
            angle_tolerance=10.0,
            json="tmp.json",
            poscar="POSCAR-tmp",
            kpt_density=4.2,
            standardize_structure=False,
            prior_info=True,
            dirs=["a", "b"],
            prev_dir="c",
            func=parsed_args.func)

        self.assertEqual(expected, parsed_args)


class MainVaspRunTest(ViseTest):
    def test_vasp_run_wo_options(self):
        parsed_args = parse_args(["vr", "-v", "vasp", "cmd"])
        # func is a pointer so need to point the same address.
        expected = Namespace(
            print=False,
            json_file="str_opt.json",
            vasp_cmd=["vasp", "cmd"],
            handler_name="default",
            timeout=TIMEOUT,
            remove_wavecar=False,
            max_relax_num=kpt_conv_defaults["max_relax_num"],
            left_files=kpt_conv_defaults["left_files"],
            func=parsed_args.func,
            **symprec_args)
        self.assertEqual(expected, parsed_args)

    def test_vasp_set_w_options(self):
        parsed_args = parse_args(["vr",
                                  "--print",
                                  "--json_file", "test.json",
                                  "-v", "vasp", "cmd",
                                  "--handler_name", "dielectric",
                                  "--timeout", "1000",
                                  "--remove_wavecar",
                                  "--max_relax_num", "10",
                                  "--left_files", "POSCAR", "PCDAT",
                                  ])
        expected = Namespace(
            print=True,
            json_file="test.json",
            vasp_cmd=["vasp", "cmd"],
            handler_name="dielectric",
            timeout=1000,
            remove_wavecar=True,
            max_relax_num=10,
            left_files=["POSCAR", "PCDAT"],
            func=parsed_args.func,
            **symprec_args)
        self.assertEqual(expected, parsed_args)


class MainKptConvTest(ViseTest):
    def test_kpt_conv_wo_options(self):
        parsed_args = parse_args(["kc", "-v", "vasp", "cmd"])
        # func is a pointer so need to point the same address.
        expected = Namespace(
            print=False,
            json_file="kpt_conv.json",
            vasp_cmd=["vasp", "cmd"],
            initial_kpt_density=KPT_INIT_DENSITY,
            handler_name="default",
            timeout=TIMEOUT,
            remove_wavecar=False,
            max_relax_num=kpt_conv_defaults["max_relax_num"],
            convergence_criterion=kpt_conv_defaults["convergence_criterion"],
            left_files=kpt_conv_defaults["left_files"],
            func=parsed_args.func,
            **default_vasp_args, **symprec_args)
        self.assertEqual(expected, parsed_args)

    def test_kpt_conv_w_options(self):
        parsed_args = parse_args(["kc",
                                  "--json_file", "test.json",
                                  "-v", "vasp", "cmd",
                                  "-ikd", "3.5",
                                  "--handler_name", "dielectric",
                                  "--timeout", "1000",
                                  "--remove_wavecar",
                                  "--max_relax_num", "10",
                                  "--criterion", "0.01",
                                  "--left_files", "POSCAR", "PCDAT"])
        expected = Namespace(
            print=False,
            json_file="test.json",
            vasp_cmd=["vasp", "cmd"],
            initial_kpt_density=3.5,
            handler_name="dielectric",
            timeout=1000,
            remove_wavecar=True,
            max_relax_num=10,
            convergence_criterion=0.01,
            left_files=["POSCAR", "PCDAT"],
            func=parsed_args.func,
            **default_vasp_args, **symprec_args)
        self.assertEqual(expected, parsed_args)


class MainChemPotDiagTest(ViseTest):
    def test_vasp_run_wo_options(self):
        parsed_args = parse_args(["cpd"])
        expected = Namespace(
            draw_phase_diagram=False,
            vasp_dirs=None,
            vasprun="vasprun.xml",
            elements=None,
            target_comp=None,
            filename=None,
            parse_gas=True,
            partial_pressures=None,
            temperature=0.0,
            func=parsed_args.func)
        self.assertEqual(expected, parsed_args)

    def test_vasp_run_w_options(self):
        parsed_args = parse_args(["cpd",
                                  "-dpd",
                                  "-d", "a_dir", "b_dir",
                                  "-v", "vasprun.xml.finish",
                                  "-e", "Mg", "O",
                                  "-c", "MgO",
                                  "-f", "cpd.pdf",
                                  "-pg", "F",
                                  "-pp", "O", "1e+5",
                                  "-t", "1000"])

        expected = Namespace(
            draw_phase_diagram=True,
            vasp_dirs=["a_dir", "b_dir"],
            vasprun="vasprun.xml.finish",
            elements=["Mg", "O"],
            target_comp="MgO",
            filename="cpd.pdf",
            parse_gas=False,
            partial_pressures=["O", "1e+5"],
            temperature=1000.0,
            func=parsed_args.func)
        self.assertEqual(expected, parsed_args)


# class MainPlotBandTest(ViseTest):
#     def
