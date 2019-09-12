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

from custodian.ansible.actions import FileActions, DictActions
from custodian.ansible.interpreter import Modder
from custodian.custodian import Job
from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VASP_BACKUP_FILES

from monty.os.path import zpath
from monty.serialization import dumpfn, loadfn
from monty.shutil import decompress_dir

from pymatgen.io.vasp import VaspInput, Poscar, Kpoints, Vasprun, Potcar
from vise.util.error_classes import EnergyNotConvergedError, \
    KptNotConvergedError

from vise.input_set.vise_incar import ViseIncar
from vise.input_set.vise_input_set import ViseInputSet
from vise.util.logger import get_logger
from vise.config import SYMMETRY_TOLERANCE

"""
This module implements basic kinds of jobs for VASP runs.
Copied from custodian's VaspJob and modified.
"""

logger = get_logger(__name__)


__author__ = "Akira Takahashi, Yu Kumagai"
__maintainer__ = "Yu Kumagai"


VASP_INPUT_FILES = {"INCAR", "POSCAR", "POTCAR", "KPOINTS"}

VASP_OUTPUT_FILES = ['DOSCAR', 'INCAR', 'KPOINTS',  'POSCAR', 'PROCAR',
                     'vasprun.xml', 'CHG', 'EIGENVAL', 'OSZICAR', 'CONTCAR',
                     'IBZKPT', 'OUTCAR']


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
    def __init__(self, vasp_cmd, gamma_vasp_cmd=None, output_file="vasp.out",
                 stderr_file="std_err.txt", suffix="", final=True,
                 backup=True, auto_continue=False):
        """
        Override constructor to close some options as some options use Incar
        instead of VisaIncar.

        Args: See docstrings of VaspJob.
        """
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

    @classmethod
    def get_runs_geom(cls,
                      vasp_cmd: Union[str, list],
                      max_relax: int,
                      gamma_vasp_cmd: Optional[str, list] = None,
                      suffix: str = "",
                      removes_outputs_in_current_dir: bool = False,
                      removes_wavecar: bool = False,
                      vaspout: str = "vasp.out"):

        if isinstance(vasp_cmd, str):
            vasp_cmd = vasp_cmd.split()

        backup = True
        for job_number in range(max_relax + 1):
            suffix = Path(suffix) if suffix else Path(".")
            last_dir = suffix / f"repeat-{job_number + 1}"
            yield cls(vasp_cmd=vasp_cmd,
                      gamma_vasp_cmd=gamma_vasp_cmd,
                      backup=backup,
                      suffix=last_dir)
            backup = False
            if len(Vasprun("vasprun.xml").ionic_steps) == 1:
                break
        else:
            raise EnergyNotConvergedError("get_runs_geom not converged")

        for n in range(1, job_number + 1):

            rdir = suffix / f"repeat-{n}"
            for f in [vaspout, "OUTCAR", "vasprun.xml", "CONTCAR"]:
                dst_f = f if f != "CONTCAR" else "POSCAR"
                if n != job_number:
                    n_str = str(n)
                else:
                    n_str = "finish"
                if "." in dst_f:
                    dst = dst_f.replace(".", "-{}.".format(n_str),
                                        1)
                else:
                    dst = "{}-{}".format(dst_f, n_str)
                shutil.copyfile("{}/{}".format(rdir, f),
                                "{}{}".format(suffix, dst))
            if n != job_number:
                shutil.rmtree(rdir)
            else:
                for f in glob("{}/*".format(rdir)):
                    if os.stat(f).st_size == 0:
                        os.remove(f)

        finish_link = "{}finish".format(temp_p)
        os.symlink("repeat-{}/".format(job_number), finish_link)
        # link finish
        if removes_outputs_in_current_dir:
            remove_duplicated_files("finish")
            # remove the following files if it's empty; otherwise move
            # to finish
            for f in ["PCDAT", "REPORT", "XDATCAR",
                      "CHG", "CHGCAR", "WAVECAR"]:
                if not os.path.exists(f):
                    continue
                if os.stat(f).st_size == 0:
                    os.remove(f)
                else:
                    # leave CHGCAR and WAVECAR at home
                    print("aaa", f)
                    if f in ["PCDAT", "REPORT", "XDATCAR", "CHG"]:
                        shutil.move(f, "finish/")
            for f in ["OUTCAR", "vasprun.xml"]:
                os.symlink("finish/{}".format(f), "./{}".format(f))

        if removes_wavecar:
            os.remove(f"{last_dir}/WAVECAR")

    @classmethod
    def a_get_runs_geom(cls,
                      vasp_cmd: Union[str, list],
                      max_relax: int,
                      gamma_vasp_cmd: Optional[str, list] = None,
                      prefix: str = "",
                      removes_outputs_in_current_dir: bool = False,
                      removes_wavecar: bool = False,
                      vaspout: str = "vasp.out"):

        vasp_cmd = vasp_cmd.split() if isinstance(vasp_cmd, str) else vasp_cmd

        converged = False
        # job_number = max_relax_number: if converged , if not, simply break
        for job_number in range(max_relax + 1):
            # set calculation condition
            if job_number == 0:
                backup = True
#                settings = []
            else:
                v = Vasprun("vasprun.xml")
                backup = False
                if len(v.ionic_steps) == 1:
                    converged = True
#                settings = [{"file": "CONTCAR",
#                             "action": {"_file_copy": {"dest": "POSCAR"}}}]

            if converged:
                # TODO: can be more simple?
                # move and remove unnecessary file
                temp_p = "{}/".format(prefix) if prefix else ""
                for n in range(1, job_number + 1):
                    rdir = "{}repeat-{}".format(temp_p, n)
                    for f in [vaspout, "OUTCAR", "vasprun.xml", "CONTCAR"]:
                        dst_f = f if f != "CONTCAR" else "POSCAR"
                        if n != job_number:
                            n_str = str(n)
                        else:
                            n_str = "finish"
                        if "." in dst_f:
                            dst = dst_f.replace(".", "-{}.".format(n_str), 1)
                        else:
                            dst = "{}-{}".format(dst_f, n_str)
                        shutil.copyfile("{}/{}".format(rdir, f),
                                        "{}{}".format(temp_p, dst))
                    if n != job_number:
                        shutil.rmtree(rdir)
                    else:
                        for f in glob("{}/*".format(rdir)):
                            if os.stat(f).st_size == 0:
                                os.remove(f)
                finish_link = "{}finish".format(temp_p)
                os.symlink("repeat-{}/".format(job_number), finish_link)
                # link finish
                if removes_outputs_in_current_dir:
                    remove_duplicated_files("finish")
                    # remove the following files if it's empty; otherwise move
                    # to finish
                    for f in ["PCDAT", "REPORT", "XDATCAR",
                              "CHG", "CHGCAR", "WAVECAR"]:
                        if not os.path.exists(f):
                            continue
                        if os.stat(f).st_size == 0:
                            os.remove(f)
                        else:
                            # leave CHGCAR and WAVECAR at home
                            print("aaa", f)
                            if f in ["PCDAT", "REPORT", "XDATCAR", "CHG"]:
                                shutil.move(f, "finish/")
                    for f in ["OUTCAR", "vasprun.xml"]:
                        os.symlink("finish/{}".format(f), "./{}".format(f))
                if removes_wavecar:
                    os.remove("{}/WAVECAR".format(last_dir))
                break
            elif job_number != max_relax:
                if prefix:
                    last_dir = "{}/repeat-{}/".format(prefix, job_number + 1)
                else:
                    last_dir = "repeat-{}/".format(job_number + 1)
                yield cls(vasp_cmd=vasp_cmd,
                          gamma_vasp_cmd=gamma_vasp_cmd,
                          backup=backup,
                          prefix=last_dir,
                          settings_override=settings)

            else:  # not converged and reach max_relax_number
                raise EnergyNotConvergedError("get_runs_geom not converged")

    @classmethod
    def kpt_converge(cls, command, max_relax, gamma_vasp_cmd=None,
                     energy_tol=0.003,
                     max_iter_kpt=10,
                     removes_wavecar=False,
                     force_override=False,
                     copy_wavecar_chgcar_to_kpt=False):
        """
        kpt converge

        """
        class CalcResult:
            def __init__(self, kpt, space_group, dir_name, structure,
                         energy_per_atom):

                self.kpt = kpt
                self.sg = space_group
                self.dir_name = dir_name
                self.structure = structure
                self.energy = energy_per_atom

            @classmethod
            def from_directory(cls, dir_name):
                # read prefix like "1-kpt8x8x8_spg225/"
                p = dir_name.split("/")[0]
                # prefix = str(prefix)
                kpt_s, spg_str = tuple(p.split("-")[1].split("_"))
                kpt = \
                    tuple(int(k) for k in kpt_s.replace("kpt", "").split("x"))
                space_group = int(spg_str.replace("spg", ""))
                # read vasprun
                v = Vasprun("{}/vasprun.xml".format(dir_name))
                structure = v.final_structure
                energy_per_atom = v.final_energy / structure.num_sites
                return cls(kpt, space_group, dir_name, structure,
                           energy_per_atom)

            @property
            def kpt_density(self):
                kpt_file = "{}/KPOINTS".format(self.dir_name)
                return 2.5

            @property
            def job_number(self):
                return int(self.dir_name.split("-")[0])

            def __str__(self):
                return '{0[0]:2d} x {0[1]:2d} x {0[2]:2d} total-kpt= {1:4d} ' \
                       'spg= {2:3d} ' \
                       'a= {3[0]:>7.4f} ' \
                       'b= {3[1]:>7.4f} ' \
                       'c= {3[2]:>7.4f} ' \
                       'alpha= {4[0]:>7.3f} ' \
                       'beta= {4[1]:>7.3f} ' \
                       'gamma= {4[2]:>7.3f} ' \
                       'energy={5:>10.4f} eV/atom\n'. \
                    format(self.kpt, reduce(mul, self.kpt, 1), self.sg,
                           self.structure.lattice.abc,
                           self.structure.lattice.angles,
                           self.energy)

        # read ViseSet from directory
        INIT_KPT_DENS = 2.5
        KPT_INCREASE_RATE = 1.2

        history = []
        # read previous calc
        for dir_ in glob("*-kpt*x*x*_spg*/"):
            logger.info("calculated_directory {} is found".format(dir_))
            if force_override:
                logger.info("force_override, then remove {}".format(dir_))
                shutil.rmtree(dir_)
            else:
                if os.path.exists("{}/finish".format(dir_)):
                    logger.info("append history {}, and keep dir".format(dir_))
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
        try:
            kpt_density = history[-1][0].kpt_density
        except:
            kpt_density = INIT_KPT_DENS

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
                    cls.get_runs_geom(command, max_relax,
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
            if copy_wavecar_chgcar_to_kpt:
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
                logger.info('converged, energy, pre_energy, pre_pre_energy '
                            '= {}, {}, {}'.
                            format(history[-1][-1].energy,
                                   history[-2][-1].energy,
                                   history[-3][-1].energy))
                os.symlink(history[-3][-1].dir_name, "finish")
                logger.info("Geometry optimization complete")
                converged_k = history[-3][-1].kpt
                with open("energy.txt", "a") as fw:
                    fw.write('Converged k-points (Criterion: {0} eV/atom): '
                             '{1[0]:2d} x '
                             '{1[1]:2d} x '
                             '{1[2]:2d} '.format(energy_tol, converged_k))
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

