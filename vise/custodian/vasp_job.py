# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
from functools import reduce
from glob import glob
from operator import mul
from typing import Union, Optional
from pathlib import Path
import numpy as np
from uuid import uuid4

from custodian.ansible.actions import FileActions, DictActions
from custodian.ansible.interpreter import Modder
from custodian.custodian import Job
from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VASP_BACKUP_FILES

from monty.json import  MSONable
from monty.serialization import dumpfn, loadfn
from monty.shutil import decompress_dir

from pymatgen.io.vasp import VaspInput, Poscar, Kpoints, Vasprun, Potcar
from pymatgen.core.structure import IStructure, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vise.util.error_classes import EnergyNotConvergedError, \
    KptNotConvergedError
from vise.input_set.incar import ViseIncar
from vise.input_set.input_set import ViseInputSet
from vise.util.logger import get_logger
from vise.config import SYMMETRY_TOLERANCE, ANGLE_TOL, KPT_INIT_DENSITY, KPT_FACTOR
"""
This module implements basic kinds of jobs for VASP runs.
Copied from custodian"s VaspJob and modified.
"""

logger = get_logger(__name__)


__author__ = "Yu Kumagai, Akira Takahashi"
__maintainer__ = "Yu Kumagai"


VASP_INPUT_FILES = {"INCAR", "POSCAR", "POTCAR", "KPOINTS"}
VASP_SAVED_FILES = ["INCAR",  "vasprun.xml", "CONTCAR", "OUTCAR"]


class GeomOptResult(MSONable):
    def __init__(self,
                 uuid: int,
                 num_kpt: list,
                 kpt_density: float,
                 initial_structure: Structure,
                 final_structure: Structure,
                 initial_sg: int,
                 final_sg: int,
                 energy_atom: float,
                 prev_geom_opt_uuid: Optional[int] = None):
        self.uuid = uuid
        self.num_kpt = list(num_kpt)
        self.kpt_density = kpt_density
        self.initial_structure = initial_structure
        self.final_structure = final_structure
        self.initial_sg = initial_sg
        self.final_sg = final_sg
        self.energy_atom = energy_atom
        self.prev_geom_opt_uuid = prev_geom_opt_uuid

    @classmethod
    def from_directory(cls,
                       dir_name: str,
                       symprec: float = SYMMETRY_TOLERANCE,
                       angle_tol: float = ANGLE_TOL,
                       prev_geom_opt: Optional["GeomOptResult"] = None):
        k = Kpoints.from_file(Path(dir_name, "KPOINTS"))
        num_kpt = k.kpts[0]

        vise = ViseInputSet.load_json(Path(dir_name, "vise.json"))
        kpt_denisty = vise.kwargs("kpt_density")

        final_structure = IStructure.from_file(Path(dir_name, "CONTCAR-finish"))
        sga = SpacegroupAnalyzer(final_structure,
                                 symprec=symprec, angle_tolerance=angle_tol)
        symmetry_dataset = sga.get_symmetry_dataset()
        final_sg = symmetry_dataset["number"]

        v = Vasprun(Path(dir_name, "vasprun.xml-finish"))
        energy_atom = v.final_energy / len(final_structure)

        if prev_geom_opt:
            prev_geom_opt_uuid = prev_geom_opt.uuid
            initial_structure = prev_geom_opt.final_structure
            initial_sg = prev_geom_opt.final_sg
        else:
            prev_geom_opt_uuid = None
            initial_structure = IStructure.from_file(Path(dir_name, "POSCAR"))
            sga = SpacegroupAnalyzer(initial_structure,
                                     symprec=symprec, angle_tolerance=angle_tol)
            initial_symmetry_dataset = sga.get_symmetry_dataset()
            initial_sg = initial_symmetry_dataset["number"]

        uuid = int(uuid4())

        return cls(uuid, num_kpt, kpt_denisty, initial_structure,
                   final_structure, initial_sg, final_sg, energy_atom,
                   prev_geom_opt_uuid)

    def is_sg_changed(self):
        return self.initial_sg == self.final_sg

    def dirname(self):
        return "kpt{}x{}x{}_sg{}_sg{}".format(*self.num_kpt, self.initial_sg, self.final_sg)

        # kpt_s, spg_str = tuple(p.split("-")[1].split("_"))
        # kpt = \
        #     tuple(int(k) for k in kpt_s.replace("kpt", "").split("x"))
        # space_group = int(spg_str.replace("spg", ""))
        # # read vasprun
        # v = Vasprun("{}/vasprun.xml".format(dir_name))
        # structure = v.final_structure
        # energy_per_atom = v.final_energy / structure.num_sites
        # return cls(kpt, space_group, dir_name, structure,
        #            energy_per_atom)

    # def __str__(self):
    #     return "{0[0]:2d} x {0[1]:2d} x {0[2]:2d} total-kpt= {1:4d} " \
    #            "sg= {2:3d} " \
    #            "a= {3[0]:>7.4f} " \
    #            "b= {3[1]:>7.4f} " \
    #            "c= {3[2]:>7.4f} " \
    #            "alpha= {4[0]:>7.3f} " \
    #            "beta= {4[1]:>7.3f} " \
    #            "gamma= {4[2]:>7.3f} " \
    #            "energy={5:>10.4f} eV/atom\n". \
    #         format(self.kpt, reduce(mul, self.kpt, 1), self.sg,
    #                self.structure.lattice.abc,
    #                self.structure.lattice.angles,
    #                self.energy)


class KptConvResult(list):

    def __init__(self, energy_tol, num):
        self.energy_tol = energy_tol
        self.num = num
        self.sg_changed = []
        super().__init__()

    def append(self, p_object: GeomOptResult):
        self.sg_changed.append(p_object.is_sg_changed)
        super().append(p_object)

    @classmethod
    def from_dirs(cls, energ_tol=10, num=10):
        xx = []
        lll = cls(energ_tol, num)
        lll += [x for x in xx if not x.prev_geom_opt_uuid]

        for i in range(len(xx)):
            lll += [x for x in xx if x.prev_geom_opt_uuid == lll[-1].uuid]

        return lll

    @property
    def is_converged(self):
        for i in self.num:
            x = -(i + 1)
            if abs(self[x].energy_aotm - self[x - 1].energy_aotm) > self.energy_tol:
                return False
        return True

    @property
    def is_sg_changed(self):
        return any(self.sg_changed)


def make_prefix(kpt: list,
                sg: int,
                max_prefix_num: int = 10,
                other_prefix: str = "") -> str:

    kpt_str = "kpt{}x{}x{}".format(*kpt)
    prefix = f"{other_prefix}{kpt_str}_sg{sg}/"

    if not os.path.exists(prefix):
        return prefix
    else:
        n = 1
        prefix_w_num = prefix
        while os.path.exists(prefix):
            if n <= max_prefix_num:
                prefix_w_num = "-".join([prefix, str(n)])
                n += 1
            else:
                raise ValueError(f"{n}th space group occur (>MAX_NUM {max_prefix_num}.)")

        return prefix_w_num


def remove_duplicated_files(result_dir):
    # post process (remove duplicated files)
    to_remove = ["CHG",
                 "CHGCAR",
                 "CONTCAR",
                 "DOSCAR",
                 "EIGENVAL",
                 "IBZKPT",
                 "INCAR",
                 "KPOINTS",
                 "OSZICAR",
                 "OUTCAR",
                 "POSCAR",
                 "WAVECAR",
                 "vasp.out",
                 "vasprun.xml"]
    for file_ in to_remove:
        if os.path.exists(f"{result_dir}/{file_}"):
            os.remove(file_)


class ViseVaspJob(VaspJob):
    def __init__(self,
                 vasp_cmd: list,
                 gamma_vasp_cmd: Optional[list] = None,
                 output_file: str = "vasp.out",
                 stderr_file: str = "std_err.txt",
                 suffix: str = "",
                 final: bool = True,
                 backup: bool = True,
                 auto_continue: bool = False):
        """
        Override constructor to close some options as some options use Incar
        instead of VisaIncar.

        Args: See docstrings of VaspJob.
        """

        if gamma_vasp_cmd is None:
            gamma_vasp_cmd = vasp_cmd[:-1]
            gamma_vasp_cmd.append(vasp_cmd[-1].replace("std", "gam"))

        # Note that only list instance is accepted for vasp_cmd in VaspJob at
        # ver.2019.8.24.
        super().__init__(vasp_cmd=vasp_cmd,
                         output_file=output_file,
                         stderr_file=stderr_file,
                         suffix=suffix,
                         final=final,
                         backup=backup,
                         auto_npar=False,
                         auto_gamma=True,
                         settings_override=None,
                         gamma_vasp_cmd=gamma_vasp_cmd,
                         copy_magmom=False,
                         auto_continue=auto_continue)

    def postprocess(self):
        """
        """
        for f in VASP_SAVED_FILES + [self.output_file]:
            if os.path.exists(f):
                shutil.copy(f, "{}{}".format(f, self.suffix))

        # Remove continuation so if a subsequent job is run in
        # the same directory, will not restart this job.
        if os.path.exists("continue.json"):
            os.remove("continue.json")

    @classmethod
    def geom_opt_run(cls,
                     vasp_cmd: list,
                     gamma_vasp_cmd: Optional[list] = None,
                     max_relax_num: int = 10,
                     removes_wavecar: bool = False,
                     vaspout: str = "vasp.out"):
        """Geometry optimization

        Args:
            vasp_cmd (str/list):
                Vasp command.
            gamma_vasp_cmd: (str/list/None)
                Gamma version of vasp command.
            max_relax_num (int):
                Maximum number of vasp calculations.
            removes_wavecar: bool = False,
                Whether to remov the WAVECAR files after geometry optimization.
            vaspout: str = "vasp.out"):
                Name of the file showing the standard output.

        Yield:
            ViseVaspJob
        """
        backup = True
        for job_number in range(1, max_relax_num + 1):
            yield cls(vasp_cmd=vasp_cmd,
                      gamma_vasp_cmd=gamma_vasp_cmd,
                      final=True,
                      backup=backup,
                      suffix=f"-{job_number}")
            shutil.copy("-".join(["CONTCAR", str(job_number)]), "POSCAR")
            backup = False

            if len(Vasprun("vasprun.xml").ionic_steps) == 1:
                break
        else:
            raise EnergyNotConvergedError("get_runs_geom not converged")

        for f in [vaspout, "OUTCAR", "vasprun.xml", "CONTCAR"]:
            shutil.move(f"{f}-{job_number}", f"{f}-finish")

        for f in glob("*"):
            if os.stat(f).st_size == 0:
                os.remove(f)

        if removes_wavecar:
            os.remove("WAVECAR")

    # TODO: move files except for inupt output files.

    @classmethod
    def kpt_converge(cls,
                     vasp_cmd: list,
                     gamma_vasp_cmd: Optional[list] = None,
                     max_relax_num: int = 10,
                     max_iter_kpt: int = 10,
                     energy_tol: float = 0.003,
                     removes_wavecar: bool = False,
                     vaspout: str = "vasp.out",
                     copy_wavecar_chgcar_to_next_kpt=False):
        """ kpt convergence

        Args:
            vasp_cmd:
            gamma_vasp_cmd:
            max_relax_num:
            removes_wavecar:
            vaspout:
                See docstrings of geom_opt_run.
            max_iter_kpt (str):
                Max kpoint iteraction number
        """
        history = []
        # read previous calc
        for dir_ in glob("*-kpt*x*x*_spg*/"):
            logger.info(f"calculated_directory {dir_} is found")
            if force_override:
                logger.info(f"force_override, then remove {dir_}")
                shutil.rmtree(dir_)
            else:
                if os.path.exists(f"{dir_}/finish"):
                    logger.info(f"append history {dir_}, and keep dir")
                    history.append(
                        [CalcResult.from_directory(
                            "{}/finish/".format(dir_))])
                else:
                    logger.info("not found finish, then remove {}".format(dir_))
                    shutil.rmtree(dir_)
        history.sort(key=lambda e: e[0].job_number)
        if os.path.exists("energy.txt"):
            shutil.move("energy.txt", "old_energy.txt")
        for calc_res_ in history:
            calc_res = calc_res_[0]
            with open("energy.txt", "a") as fw:
                fw.write(str(calc_res))
        logger.info("history = {}".format(history))
        kpt_converged = False

        kpt_density = getattr(history[-1][0], "kpt_density", INIT_KPT_DENS)

        while (not kpt_converged) and (len(history) < max_iter_kpt):

            calc_count = len(history) + 1
            oba_vis = \
                ViseInputSet.from_prev_calc(".",
                                            parse_calc_results=False,
                                            standardize_structure=True,
                                            parse_incar=True,
                                            parse_kpoints=False,
                                            kpt_density=kpt_density)
            sg = oba_vis.sg
            k = oba_vis.kpoints.kpts[0]

            # change kpoints
            if len(history) == 0:  # initial calculation
                pass
            elif sg != history[-1][-1].sg:  # space group changed
                logger.info("space group changed from {} to {}".
                            format(history[-1][-1].sg, sg))

                # When sg is changed, WAVECAR must be removed.
                if os.path.exists("WAVECAR"):
                    os.remove("WAVECAR")

                kpt_density = INIT_KPT_DENS
            else:  # increase kpt
                pre_k = history[-1][-1].kpt
                while any(n < pre_n for n, pre_n in zip(k, pre_k)) or \
                        k == pre_k:
                    kpt_density *= KPT_INCREASE_RATE
                    oba_vis = \
                        ViseInputSet.from_prev_calc(".",
                                                    parse_calc_results=False,
                                                    standardize_structure=True,
                                                    parse_incar=True,
                                                    parse_kpoints=False,
                                                    kpt_density=kpt_density)
                    k = oba_vis.kpoints.kpts[0]

            # do vasp calculation
            oba_vis.write_input(".")
            geom_history = []
            prefix = make_prefix(k, sg, other_prefix=str(calc_count))
            for repeat, job, in enumerate(
                    cls.get_runs_geom(vasp_cmd, max_relax_num,
                                      gamma_vasp_cmd=gamma_vasp_cmd,
                                      suffix=prefix,
                                      removes_outputs_in_current_dir=False,
                                      removes_wavecar=removes_wavecar)):
                logger.info("prefix(in kpt) = {}".format(prefix))
                logger.info("repeat = {}".format(repeat))
                yield job
                dir_ = prefix + "/repeat-{}/".format(repeat+1)
                v = Vasprun("{}/vasprun.xml".format(dir_))

                st = v.final_structure
                energy = v.final_energy / st.num_sites
                result = CalcResult(k, sg, dir_, st, energy)
                geom_history.append(result)

            # Move WAVECAR & CHGCAR to ?-kpt?x?x?-sg???
            # Copy such huge files might not be good ideas.
            if copy_wavecar_chgcar_to_next_kpt:
                for f in ["WAVECAR", "CHGCAR"]:
                    if os.stat(f).st_size > 0:
                        shutil.copy(f, prefix)

            with open("energy.txt", "a") as fw:
                fw.write(str(geom_history[-1]))
            history.append(geom_history)

            # check convergence
            def is_converged(calc_res1, calc_res2):
                if calc_res1.sg != calc_res2.sg:
                    return False
                if abs(calc_res1.energy - calc_res2.energy) > energy_tol:
                    return False
                lat1, lat2 = \
                    calc_res1.structure.lattice, calc_res2.structure.lattice
                if not np.allclose(lat1.matrix, lat2.matrix, atol=SYMMETRY_TOLERANCE):
                    return False
                return True

            # 1st or 2nd calc
            if len(history) <= 2:
                continue
            # converged
            elif is_converged(history[-1][-1], history[-3][-1]) and \
                    is_converged(history[-2][-1], history[-3][-1]):
                logger.info("converged, energy, pre_energy, pre_pre_energy "
                            "= {}, {}, {}".
                            format(history[-1][-1].energy,
                                   history[-2][-1].energy,
                                   history[-3][-1].energy))
                os.symlink(history[-3][-1].dir_name, "finish")
                logger.info("Geometry optimization complete")
                converged_k = history[-3][-1].kpt
                with open("energy.txt", "a") as fw:
                    fw.write("Converged k-points (Criterion: {0} eV/atom): "
                             "{1[0]:2d} x "
                             "{1[1]:2d} x "
                             "{1[2]:2d} ".format(energy_tol, converged_k))
                # if removes_outputs_in_current_dir:
                remove_duplicated_files("finish")
                for f in ["PCDAT", "REPORT", "XDATCAR", "CHG", "CHGCAR",
                          "WAVECAR"]:
                    if not os.path.exists(f):
                        continue
                    if os.stat(f).st_size == 0:
                        os.remove(f)
                    else:
                        # leave CHGCAR and WAVECAR at home
                        if f in ["PCDAT", "REPORT", "XDATCAR", "CHG"]:
                            shutil.move(f, "finish/")

                # Keep WAVECAR obtained at the converged kpoints.
                for dir_ in glob("*-kpt*x*x*_spg*"):
                    kpt_conv_dir = history[-3][-1].dir_name.split("/")[0]
                    for f in ["CHGCAR", "WAVECAR"]:
                        ff = os.path.join(dir_, f)
                        if os.path.exists(ff) and os.stat(f).st_size > 0:
                            if dir_ == kpt_conv_dir:
                                print("{} is moved from {}".format(f, dir_))
                                if os.path.exists(f):
                                    os.remove(f)
                                shutil.move(ff, ".")
                            else:
                                os.remove(ff)

                for f in ["OUTCAR", "vasprun.xml", "KPOINTS", "INCAR"]:
                    os.symlink("finish/{}".format(f), "./{}".format(f))

                shutil.copy("finish/CONTCAR", "POSCAR-finish")
                shutil.copy("org/POSCAR", "POSCAR")
                break
            # not converged
            else:
                logger.info("not converged, energy, pre_energy, "
                            "pre_pre_energy = {}, {}, {}".
                            format(history[-1][-1].energy,
                                   history[-2][-1].energy,
                                   history[-3][-1].energy))
        if len(history) == max_iter_kpt:
            raise KptNotConvergedError()

