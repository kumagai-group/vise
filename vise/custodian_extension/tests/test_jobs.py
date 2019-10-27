# -*- coding: utf-8 -*-
import shutil
import tempfile
from pathlib import Path
from glob import glob

from pymatgen.core.structure import Structure

from vise.util.testing import ViseTest
from vise.custodian_extension.jobs import rm_wavecar, StructureOptResult, KptConvResult

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class RmWavecarTest(ViseTest):

    def setUp(self) -> None:
        # Create a temporary directory
        Path("tmp").mkdir()
        self.test_dir = Path("tmp")

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def test_remove_wavecar_at_current_dir(self):
        # Whether WAVECAR at the current directory exists or not.
        Path('WAVECAR').touch()
        rm_wavecar(remove_current=True)
        self.assertFalse(Path("WAVECAR").is_file())

    def test_remove_wavecar_at_subdir(self):
        # Create a file in the temporary directory
        f = self.test_dir / 'WAVECAR'
        f.touch()

        rm_wavecar(remove_current=False, remove_subdirectories=True)
        self.assertFalse(f.is_file())


class StructureOptResultTest(ViseTest):

    def setUp(self) -> None:
        # Create a temporary directory
        Path("tmp").mkdir()
        self.test_dir = Path("tmp")

        p = Path("MgO") / "kpt7x7x7_pre-sg225_pos-sg225" / "files"
        poscar = p / "POSCAR.orig"
        contcar = p / "CONTCAR"
        initial_structure = Structure.from_file(poscar)
        final_structure = Structure.from_file(contcar)
        self.result = StructureOptResult(uuid=1234,
                                         energy_per_atom=-5.9558474,
                                         num_kpt=[7, 7, 7],
                                         final_structure=final_structure,
                                         final_sg=225,
                                         kpt_density=2.5,
                                         initial_structure=initial_structure,
                                         initial_sg=225)

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.result.to_json_file(tmp_file.name)
        prior_info_from_json = StructureOptResult.load_json(tmp_file.name)
        self.assertEqual(prior_info_from_json.as_dict(), self.result.as_dict())

    def test_msonable(self):
        self.assertMSONable(self.result)

    def test_from_dir(self):
        from_dir = Path("MgO") / "kpt7x7x7_pre-sg225_pos-sg225"
        for f in ["KPOINTS.orig", "POSCAR.orig", "CONTCAR", "vasprun.xml"]:
            shutil.copy(str(from_dir / "files" / f), str(self.test_dir / f))

        with self.assertRaises(FileNotFoundError):
            StructureOptResult.from_dir(dir_name="tmp",
                                        kpoints="KPOINTS.orig",
                                        poscar="POSCAR.orig",
                                        contcar="CONTCAR",
                                        vasprun="vasprun.xml")

        shutil.copy(str(from_dir / "vise.json"), str(self.test_dir))

        r = StructureOptResult.from_dir(dir_name="tmp",
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
        # Create a temporary directory
#        Path("tmp").mkdir()
#        self.test_dir = Path("tmp")

        p = Path("MgO") / "kpt7x7x7_pre-sg225_pos-sg225" / "files"
        poscar = p / "POSCAR.orig"
        contcar = p / "CONTCAR"
        initial_structure = Structure.from_file(poscar)
        final_structure = Structure.from_file(contcar)
        sor1 = StructureOptResult(uuid=60100624528902743369273616597818673779,
                                  energy_per_atom=-5.9558474,
                                  num_kpt=[7, 7, 7],
                                  final_structure=final_structure,
                                  final_sg=225,
                                  kpt_density=2.5,
                                  initial_structure=initial_structure,
                                  initial_sg=225)
        sor2 = StructureOptResult(uuid=97337091404756212407260371561032442781,
                                  energy_per_atom=-5.955879895,
                                  num_kpt=[8, 8, 8],
                                  final_structure=final_structure,
                                  final_sg=225,
                                  kpt_density=3.0,
                                  initial_structure=initial_structure,
                                  initial_sg=225,
                                  prev_structure_opt_uuid=60100624528902743369273616597818673779)
        sor3 = StructureOptResult(uuid=34316653945151013128859614474515815310,
                                  energy_per_atom=-5.955888605,
                                  num_kpt=[10, 10, 10],
                                  final_structure=final_structure,
                                  final_sg=225,
                                  kpt_density=3.6,
                                  initial_structure=initial_structure,
                                  initial_sg=225,
                                  prev_structure_opt_uuid=97337091404756212407260371561032442781)

        self.kpt_conv = KptConvResult(str_opts=[sor1, sor2, sor3],
                                      convergence_energy_criterion=0.01,
                                      num_kpt_check=2,
                                      symprec=0.01,
                                      angle_tolerance=5)

    def test_json(self):
        """ round trip test of to_json and from_json """
        print(self.kpt_conv.as_dict())

        # tmp_file = tempfile.NamedTemporaryFile()
        # self.result.to_json_file(tmp_file.name)
        # prior_info_from_json = StructureOptResult.load_json(tmp_file.name)
        # self.assertEqual(prior_info_from_json.as_dict(), self.result.as_dict())

