# coding: utf-8

from __future__ import unicode_literals, division
from functools import reduce
import subprocess
from glob import glob
import numpy as np
from operator import mul
import os
import shutil
import six
import logging

from pymatgen.io.vasp import VaspInput, Poscar, Kpoints, Vasprun, Potcar
from monty.os.path import zpath
from monty.shutil import decompress_dir
from monty.serialization import dumpfn, loadfn

from custodian.ansible.actions import FileActions, DictActions
from custodian.ansible.interpreter import Modder
from custodian.custodian import Job
from custodian.utils import backup
from custodian.vasp.handlers import VASP_BACKUP_FILES

from obadb.vasp.input_set import ObaSet, ObaIncar
from obadb.vasp.kpoints import get_kpt_dens_from_file
from obadb.util.structure_handler import SYMPREC

"""
This module implements basic kinds of jobs for VASP runs.
Copied from custodian's VaspJob and modified.
"""


logger = logging.getLogger(__name__)


__author__ = "Akira Takahashi, Yu Kumagai"
__maintainer__ = "Akira Takahashi"


VASP_INPUT_FILES = {"INCAR", "POSCAR", "POTCAR", "KPOINTS"}

VASP_OUTPUT_FILES = ['DOSCAR', 'INCAR', 'KPOINTS',  'POSCAR', 'PROCAR',
                     'vasprun.xml', 'CHG', 'EIGENVAL', 'OSZICAR', 'CONTCAR',
                     'IBZKPT', 'OUTCAR']
# VASP_OUTPUT_FILES = ['DOSCAR', 'POSCAR', 'PROCAR',
#                      'vasprun.xml', 'CHGCAR', 'CHG', 'EIGENVAL', 'OSZICAR',
#                      'WAVECAR', 'CONTCAR', 'IBZKPT', 'OUTCAR']


def is_empty_file(filename):
    return os.stat(filename).st_size == 0


def kpt_str(kpt):
    return "kpt{}x{}x{}".format(kpt[0], kpt[1], kpt[2])


def make_prefix(kpt, spg, other_prefix=""):
    if other_prefix:
        other_prefix += "-"
    prefix = "{}{}_spg{}/".\
        format(other_prefix, kpt_str(kpt), spg)
    if not os.path.exists(prefix):
        return prefix
    else:
        n = 1
        MAX_NUM = 10
        while os.path.exists(prefix) and n <= MAX_NUM:
            prefix = "{}{}_spg{}-{}/". \
                format(other_prefix, kpt_str(kpt), spg, n)
            n += 1
        if n > MAX_NUM:
            raise ValueError("{} th space group occur (>MAX_NUM {})".
                             format(n, MAX_NUM))
        return prefix


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
        if os.path.exists("{}/{}".format(result_dir, file_)):
            os.remove(file_)


class ObaVaspModder(Modder):
    def __init__(self, actions=None, strict=True, vi=None):
        """
        Initializes a Modder for VaspInput sets

        Args:
            actions ([Action]): A sequence of supported actions. See
                :mod:`custodian.ansible.actions`. Default is None,
                which means DictActions and FileActions are supported.
            strict (bool): Indicating whether to use strict mode. In non-strict
                mode, unsupported actions are simply ignored without any
                errors raised. In strict mode, if an unsupported action is
                supplied, a ValueError is raised. Defaults to True.
            vi (VaspInput): A VaspInput object from the current directory.
                Initialized automatically if not passed (but passing it will
                avoid having to reparse the directory).
        """
        self.vi = vi or ObaVaspInput.from_directory('.')
        actions = actions or [FileActions, DictActions]
        super(ObaVaspModder, self).__init__(actions, strict)

    def apply_actions(self, actions):
        """
        Applies a list of actions to the Vasp Input Set and rewrites modified
        files.
        Args:
            actions [dict]: A list of actions of the form {'file': filename,
                'action': moddermodification} or {'dict': vaspinput_key,
                'action': moddermodification}
        """
        modified = []
        for a in actions:
            if "dict" in a:
                k = a["dict"]
                modified.append(k)
                self.vi[k] = self.modify_object(a["action"], self.vi[k])
            elif "file" in a:
                self.modify(a["action"], a["file"])
            else:
                raise ValueError("Unrecognized format: {}".format(a))
        for f in modified:
            if f == "INCAR":
                oba_incar = ObaIncar.from_dict(self.vi["INCAR"].as_dict())
                oba_incar.write_file("INCAR")
            else:
                self.vi[f].write_file(f)


class ObaVaspInput(VaspInput):

    @staticmethod
    def from_directory(input_dir, optional_files=None):
        """
        Read in a set of VASP input from a directory. Note that only the
        standard INCAR, POSCAR, POTCAR and KPOINTS files are read unless
        optional_filenames is specified.

        Args:
            input_dir (str): Directory to read VASP input from.
            optional_files (dict): Optional files to read in as well as a
                dict of {filename: Object type}. Object type must have a
                static method from_file.
        """
        sub_d = {}
        for fname, ftype in [("INCAR", ObaIncar), ("KPOINTS", Kpoints),
                             ("POSCAR", Poscar), ("POTCAR", Potcar)]:
            fullzpath = zpath(os.path.join(input_dir, fname))
            sub_d[fname.lower()] = ftype.from_file(fullzpath)
        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = \
                    ftype.from_file(os.path.join(input_dir, fname))
        return ObaVaspInput(**sub_d)


class ObaVaspJob(Job):
    """
    A basic vasp job. Just runs whatever is in the directory. But conceivably
    can be a complex processing of inputs etc. with initialization.
    """

    def __init__(self, vasp_cmd, gamma_vasp_cmd=None, output_file="vasp.out",
                 stderr_file="std_err.txt", prefix="", suffix="",
                 backup=True, settings_override=None, auto_continue=False):
        """
        This constructor is necessarily complex due to the need for
        flexibility. For standard kinds of runs, it's often better to use one
        of the static constructors. The defaults are usually fine too.

        Args:
            vasp_cmd (str): Command to run vasp as a list of args. For example,
                if you are using mpirun, it can be something like
                ["mpirun", "pvasp.5.2.11"]
            output_file (str): Name of file to direct standard out to.
                Defaults to "vasp.out".
            stderr_file (str): Name of file to direct standard error to.
                Defaults to "std_err.txt".
            prefix (str): A prefix to be appended to the beginning output.
                E.g., to rename all VASP output from say vasp.output to
                kpt1x1x1-relax1/vasp.out, provide "kpt1x1x1-relax1/" as the
                prefix. If prefix includes "/", the directory will be made.
            suffix (str): A suffix to be appended to the final output. E.g.,
                to rename all VASP output from say vasp.out to
                vasp.out.relax1, provide ".relax1" as the suffix.
            backup (bool): Whether to backup the initial input files. If True,
                the INCAR, KPOINTS, POSCAR and POTCAR will be copied with a
                ".org" appended. Defaults to True.
            settings_override ([dict]): An ansible style list of dict to
                override changes. For example, to set ISTART=1 for subsequent
                runs and to copy the CONTCAR to the POSCAR, you will provide::

                    [{"dict": "INCAR", "action": {"_set": {"ISTART": 1}}},
                     {"file": "CONTCAR",
                      "action": {"_file_copy": {"dest": "POSCAR"}}}]
            auto_continue (bool): Whether to automatically continue a run
                if a STOPCAR is present. This is very usefull if using the
                wall-time handler which will write a read-only STOPCAR to
                prevent VASP from deleting it once it finishes
        """
        self.vasp_cmd = vasp_cmd
        if isinstance(self.vasp_cmd, six.string_types):
            self.vasp_cmd = self.vasp_cmd.split()
        self.gamma_vasp_cmd = gamma_vasp_cmd
        if isinstance(self.gamma_vasp_cmd, six.string_types):
            self.gamma_vasp_cmd = self.gamma_vasp_cmd.split()
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.backup = backup
        self.prefix = prefix
        self.suffix = suffix
        self.settings_override = settings_override
        self.auto_continue = auto_continue

    def setup(self):
        """
        Performs initial setup for VaspJob, including overriding any settings
        and backing up.
        """
        decompress_dir('.')

        if self.backup:
            for f in VASP_INPUT_FILES:
                dst_name = "org/{}".format(f)
                if os.path.dirname(dst_name):
                    os.makedirs(os.path.dirname(dst_name), exist_ok=True)
                shutil.copy(f, dst_name)

        if self.auto_continue:
            if os.path.exists("continue.json"):
                actions = loadfn("continue.json").get("actions")
                logger.info("Continuing previous VaspJob. Actions: {}".format(actions))
                backup(VASP_BACKUP_FILES, prefix="prev_run")
                ObaVaspModder().apply_actions(actions)

            else:
                # Default functionality is to copy CONTCAR to POSCAR and set
                # ISTART to 1 in the INCAR, but other actions can be specified
                if self.auto_continue is True:
                    actions = [{"file": "CONTCAR",
                                "action": {"_file_copy": {"dest": "POSCAR"}}},
                               {"dict": "INCAR",
                                "action": {"_set": {"ISTART": 1}}}]
                else:
                    actions = self.auto_continue
                dumpfn({"actions": actions}, "continue.json")

        if self.settings_override is not None:
            ObaVaspModder().apply_actions(self.settings_override)

    def run(self):
        """
        Perform the actual VASP run.

        Returns:
            (subprocess.Popen) Used for monitoring.
        """
        # check existence and validation of input files
        for file_name, check_cls in [("POSCAR", Poscar),
                                     ("INCAR", ObaIncar),
                                     ("KPOINTS", Kpoints),
                                     ("POTCAR", Potcar)]:
            if not os.path.exists(file_name):
                raise FileNotFoundError("{} was not found".format(file_name))
            try:
                check_cls.from_file(file_name)
            except:
                raise ValueError("Invalid {}".format(file_name))

        # determine whether vasp_std or vasp_gam is used
        cmd = list(self.vasp_cmd)
        if self.gamma_vasp_cmd:
            kpoints = Kpoints.from_file("KPOINTS")
            if kpoints.style == Kpoints.supported_modes.Gamma \
                    and tuple(kpoints.kpts[0]) == (1, 1, 1):
                cmd = list(self.gamma_vasp_cmd)

        # do vasp calculation
        logger.info("Running {}".format(" ".join(cmd)))
        with open(self.output_file, 'w') as f_std, \
                open(self.stderr_file, "w", buffering=1) as f_err:
            # use line buffering for stderr
            p = subprocess.Popen(cmd, stdout=f_std, stderr=f_err)
        return p

    def postprocess(self):
        """
        Postprocessing includes renaming and gzipping where necessary.
        Also copies the magmom to the incar if necessary
        """
        if self.prefix or self.suffix:
            for f in VASP_OUTPUT_FILES + [self.output_file]:
                if os.path.exists(f):
                    dst_name = "{}{}{}".format(self.prefix, f, self.suffix)
                    if os.path.dirname(dst_name):
                        os.makedirs(os.path.dirname(dst_name), exist_ok=True)
                    shutil.copy(f, dst_name)

        # Remove continuation so if a subsequent job is run in
        # the same directory, will not restart this job.
        if os.path.exists("continue.json"):
            os.remove("continue.json")

    def terminate(self):
        for k in self.vasp_cmd:
            if "vasp" in k:
                try:
                    os.system("killall %s" % k)
                except:
                    pass

    @classmethod
    def get_runs_geom(cls, command, max_relax, gamma_vasp_cmd=None, prefix="",
                      removes_outputs_in_current_dir=False,
                      removes_wavecar=False,
                      vaspout="vasp.out"):

        if isinstance(command, six.string_types):
            vasp_command = command.split()
        else:
            vasp_command = command
        converged = False
        # job_number = max_relax_number: if converged , if not, simply break
        for job_number in range(max_relax + 1):
            # set calculation condition
            if job_number == 0:
                backup = True
                settings = []
            else:
                v = Vasprun("vasprun.xml")
                backup = False
                if len(v.ionic_steps) == 1:
                    converged = True
                settings = [{"file": "CONTCAR",
                             "action": {"_file_copy": {"dest": "POSCAR"}}}]

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
                            if is_empty_file(f):
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
                        if is_empty_file(f):
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
                yield cls(vasp_command, gamma_vasp_cmd=gamma_vasp_cmd,
                          backup=backup, prefix=last_dir,
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
                return get_kpt_dens_from_file(kpt_file)

            @property
            def job_number(self):
                return int(self.dir_name.split("-")[0])

            def __str__(self):
                return '{0[0]:2d} x {0[1]:2d} x {0[2]:2d} total-kpt= {1:4d} '\
                        'spg= {2:3d} '\
                        'a= {3[0]:>7.4f} '\
                        'b= {3[1]:>7.4f} '\
                        'c= {3[2]:>7.4f} '\
                        'alpha= {4[0]:>7.3f} '\
                        'beta= {4[1]:>7.3f} '\
                        'gamma= {4[2]:>7.3f} '\
                        'energy={5:>10.4f} eV/atom\n'.\
                        format(self.kpt, reduce(mul, self.kpt, 1), self.sg,
                               self.structure.lattice.abc,
                               self.structure.lattice.angles,
                               self.energy)

        # read ObaSet from directory
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
                ObaSet.from_prev_calc(".",
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
                while any(n < pre_n for n, pre_n in zip(k, pre_k)) or\
                        k == pre_k:
                    kpt_density *= KPT_INCREASE_RATE
                    oba_vis = \
                        ObaSet.from_prev_calc(".",
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
                                      prefix=prefix,
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
                    if not is_empty_file(f):
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
                if not np.allclose(lat1.matrix, lat2.matrix, atol=SYMPREC):
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
                logging.info("Geometry optimization complete")
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
                    if is_empty_file(f):
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
                        if os.path.exists(ff) and not is_empty_file(ff):
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


class EnergyNotConvergedError(Exception):
    pass


class KptNotConvergedError(Exception):
    pass
