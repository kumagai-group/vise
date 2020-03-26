# -*- coding: utf-8 -*-
import shutil
import tempfile
from pathlib import Path
import os

from pymatgen.core.structure import Structure

from vise.util.testing import ViseTest
from vise.custodian_extension.jobs import (
    rm_wavecar, StructureOptResult, KptConvResult, ViseVaspJob)


parent_dir = Path(__file__).parent


class RmWavecarTest(ViseTest):

    def test_remove_wavecar_at_current_dir(self):
        # Whether WAVECAR at the current directory exists or not.
        with tempfile.TemporaryDirectory() as tmp_dirname:
            os.chdir(tmp_dirname)
            Path('WAVECAR').touch()
            rm_wavecar(remove_current=True)
            self.assertFalse(Path("WAVECAR").is_file())

    def test_remove_wavecar_at_subdir(self):
        # Create a file in the temporary directory
        with tempfile.TemporaryDirectory() as tmp_dirname:
            os.chdir(tmp_dirname)
            tmp = Path("tmp")
            tmp.mkdir()
            f = tmp / 'WAVECAR'
            f.touch()
            rm_wavecar(remove_current=False, remove_subdirectories=True)
            self.assertFalse(f.is_file())


class StructureOptResultTest(ViseTest):

    def setUp(self) -> None:
        p = parent_dir / "MgO" / "kpt7x7x7_pre-sg225_pos-sg225" / "files"
        poscar = p / "POSCAR.orig"
        contcar = p / "CONTCAR"
        initial_structure = Structure.from_file(poscar)
        final_structure = Structure.from_file(contcar)
        self.result = StructureOptResult(uuid=1234,
                                         energy_per_atom=-5.95584738,
                                         num_kpt=[7, 7, 7],
                                         final_structure=final_structure,
                                         final_sg=225,
                                         kpt_density=2.5,
                                         initial_structure=initial_structure,
                                         initial_sg=225)

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.result.to_json_file(tmp_file.name)
        prior_info_from_json = StructureOptResult.load_json(tmp_file.name)
        self.assertEqual(prior_info_from_json.as_dict(), self.result.as_dict())

    def test_msonable(self):
        self.assertMSONable(self.result)

    def test_from_dir(self):
        from_dir = parent_dir / "MgO" / "kpt7x7x7_pre-sg225_pos-sg225"
        with tempfile.TemporaryDirectory() as tmp_dirname:

            for f in ["KPOINTS.orig", "POSCAR.orig", "CONTCAR", "vasprun.xml"]:
                shutil.copy(str(from_dir / "files" / f),
                            "/".join([tmp_dirname, f]))

            with self.assertRaises(FileNotFoundError):
                StructureOptResult.from_dir(dir_name=parent_dir / tmp_dirname,
                                            kpoints="KPOINTS.orig",
                                            poscar="POSCAR.orig",
                                            contcar="CONTCAR",
                                            vasprun="vasprun.xml")

            shutil.copy(str(from_dir / "vise.json"), tmp_dirname)

            r = StructureOptResult.from_dir(dir_name=parent_dir / tmp_dirname,
                                            kpoints="KPOINTS.orig",
                                            poscar="POSCAR.orig",
                                            contcar="CONTCAR",
                                            vasprun="vasprun.xml")
            r.uuid = 1234
            self.assertEqual(r.as_dict(), self.result.as_dict())

    def test_is_sg_changed(self):
        self.assertFalse(self.result.is_sg_changed)


class KptConvResultTest(ViseTest):
    def setUp(self) -> None:
        p = parent_dir / "MgO" / "kpt7x7x7_pre-sg225_pos-sg225" / "files"
        poscar = p / "POSCAR.orig"
        contcar = p / "CONTCAR"
        initial_structure = Structure.from_file(poscar)
        final_structure = Structure.from_file(contcar)
        sor1 = StructureOptResult(uuid=307128884155949620271696925608470838559,
                                  energy_per_atom=-5.95584738,
                                  num_kpt=[7, 7, 7],
                                  final_structure=final_structure,
                                  final_sg=225,
                                  kpt_density=2.5,
                                  initial_structure=initial_structure,
                                  initial_sg=225)
        sor2 = StructureOptResult(uuid=64782881856798798965584516606915682920,
                                  energy_per_atom=-5.95587987,
                                  num_kpt=[8, 8, 8],
                                  final_structure=final_structure,
                                  final_sg=225,
                                  kpt_density=3.0,
                                  initial_structure=initial_structure,
                                  initial_sg=225,
                                  prev_structure_opt_uuid=307128884155949620271696925608470838559)
        sor3 = StructureOptResult(uuid=213375410158051128414604616357653492786,
                                  energy_per_atom=-5.95588858,
                                  num_kpt=[10, 10, 10],
                                  final_structure=final_structure,
                                  final_sg=225,
                                  kpt_density=3.5999999999999996,
                                  initial_structure=initial_structure,
                                  initial_sg=225,
                                  prev_structure_opt_uuid=64782881856798798965584516606915682920)

        self.kpt_conv = KptConvResult(str_opts=[sor1, sor2, sor3],
                                      convergence_criterion=0.01,
                                      num_kpt_check=2,
                                      symprec=0.01,
                                      angle_tolerance=5)

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.kpt_conv.to_json_file(tmp_file.name)
        prior_info_from_json = StructureOptResult.load_json(tmp_file.name)
        self.assertEqual(
            prior_info_from_json.as_dict(), self.kpt_conv.as_dict())

    def test_from_dirs(self):
        kpt_conv = KptConvResult.from_dirs(
            convergence_criterion=0.01,
            num_kpt_check=2,
            symprec=0.01,
            angle_tolerance=5,
            dirs=[parent_dir / "MgO/kpt8x8x8_pre-sg225_pos-sg225",
                  parent_dir / "MgO/kpt7x7x7_pre-sg225_pos-sg225",
                  parent_dir / "MgO/kpt10x10x10_pre-sg225_pos-sg225"])

        self.assertEqual(kpt_conv.as_dict(), self.kpt_conv.as_dict())

    def test_print(self):
        print(self.kpt_conv)


class ViseVaspJobTest(ViseTest):

    def setUp(self) -> None:
        self.vise_vasp_job = ViseVaspJob(vasp_cmd=["vasp"])

    def test_msonable(self):
        self.assertMSONable(self.vise_vasp_job)

