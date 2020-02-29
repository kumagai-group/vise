# -*- coding: utf-8 -*-
import subprocess
from collections import Counter
from copy import deepcopy

from custodian.custodian import ErrorHandler
from custodian.utils import backup
from custodian.vasp import handlers as orig_handlers
from pymatgen.io.vasp import Vasprun
from vise.input_set.incar import ViseIncar
from vise.input_set.vasp_input import ViseVaspInput
from vise.custodian_extension.modder import ViseVaspModder
from pymatgen.io.vasp.inputs import VaspInput, Incar, Kpoints
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.transformations.standard_transformations import \
    SupercellTransformation

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

# Monkey patch
orig_handlers.Incar = ViseIncar
orig_handlers.VaspModder = ViseVaspModder
orig_handlers.VaspInput = ViseVaspInput


# When the original error handlers are modified, add *Vise* to the class name.
class ViseVaspErrorHandler(orig_handlers.VaspErrorHandler):

    error_msgs = deepcopy(orig_handlers.VaspErrorHandler.error_msgs)
    error_msgs["plane_wave_coeff"] = \
        ["ERROR: while reading WAVECAR, plane wave coefficients changed"]
    error_msgs.pop("brmix", None)

    def __init__(self,  output_filename="vasp.out", natoms_large_cell=100):
        super().__init__(output_filename=output_filename,
                         natoms_large_cell=natoms_large_cell,
                         errors_subset_to_catch=ViseVaspErrorHandler.error_msgs)

    def check(self):
        incar = Incar.from_file("INCAR")
        self.errors = set()
        errors = set()
        with open(self.output_filename, "r") as f:
            for line in f:
                l = line.strip()
                for err, msgs in self.error_msgs.items():
                    if err in self.errors_subset_to_catch:
                        for msg in msgs:
                            if l.find(msg) != -1:
                                # this checks if we want to run a charged
                                # computation (e.g., defects) if yes we don't
                                # want to kill it because there is a change in
                                # e-density (brmix error)
                                if err == "brmix" and 'NELECT' in incar:
                                    continue
                                self.errors.add(err)
                                errors.add(msg)
        for msg in errors:
            self.logger.error(msg, extra={"incar": incar.as_dict()})
        return len(self.errors) > 0

    def correct(self):

        backup(orig_handlers.VASP_BACKUP_FILES | {self.output_filename})
        actions = []
        vi = VaspInput.from_directory(".")

        if self.errors.intersection(["tet", "dentet"]):
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"ISMEAR": 0}}})

        if "inv_rot_mat" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"SYMPREC": 1e-8}}})

        # ----- added ---------------------------------------------------
        if "plane_wave_coeff" in self.errors:
            actions.append({"file": "WAVECAR",
                            "action": {
                                "_file_delete": {'mode': "actual"}}})
            actions.append({"file": "CHGCAR",
                            "action": {
                                "_file_delete": {'mode': "actual"}}})
        # ---------------------------------------------------------------

        if "zpotrf" in self.errors:
            # Usually caused by short bond distances. If on the first step,
            # volume needs to be increased. Otherwise, it was due to a step
            # being too big and POTIM should be decreased.  If a static run
            # try turning off symmetry.
            try:
                oszicar = Oszicar("OSZICAR")
                nsteps = len(oszicar.ionic_steps)
            except:
                nsteps = 0

            if nsteps >= 1:
                potim = float(vi["INCAR"].get("POTIM", 0.5)) / 2.0
                actions.append(
                    {"dict": "INCAR",
                     "action": {"_set": {"ISYM": 0, "POTIM": potim}}})
            elif vi["INCAR"].get("NSW", 0) == 0 \
                    or vi["INCAR"].get("ISIF", 0) in range(3):
                actions.append(
                    {"dict": "INCAR", "action": {"_set": {"ISYM": 0}}})
            else:
                s = vi["POSCAR"].structure
                s.apply_strain(0.2)
                actions.append({"dict": "POSCAR",
                                "action": {
                                    "_set": {"structure": s.as_dict()}}})

            # Based on VASP forum's recommendation, you should delete the
            # CHGCAR and WAVECAR when dealing with this error.
            if vi["INCAR"].get("ICHARG", 0) < 10:
                actions.append({"file": "CHGCAR",
                                "action": {
                                    "_file_delete": {'mode': "actual"}}})
                actions.append({"file": "WAVECAR",
                                "action": {
                                    "_file_delete": {'mode': "actual"}}})

        if self.errors.intersection(["subspacematrix"]):
            if self.error_count["subspacematrix"] == 0:
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"LREAL": False}}})
            else:
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"PREC": "Accurate"}}})
            self.error_count["subspacematrix"] += 1

        if self.errors.intersection(
                ["rspher", "real_optlay", "nicht_konv"]):
            s = vi["POSCAR"].structure
            if len(s) < self.natoms_large_cell:
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"LREAL": False}}})
            else:
                # for large supercell, try an in-between option LREAL = True
                # prior to LREAL = False
                if self.error_count['real_optlay'] == 0:
                    # use real space projectors generated by pot
                    actions.append({"dict": "INCAR",
                                    "action": {"_set": {"LREAL": True}}})
                elif self.error_count['real_optlay'] == 1:
                    actions.append({"dict": "INCAR",
                                    "action": {"_set": {"LREAL": False}}})
                self.error_count['real_optlay'] += 1

        if self.errors.intersection(["tetirr", "incorrect_shift"]):

            # --Modified------------------------------------------------------
            if vi["KPOINTS"].style == Kpoints.supported_modes.Monkhorst or \
                    vi["KPOINTS"].kpts_shift != [0.0, 0.0, 0.0]:
                actions.append({"dict": "KPOINTS",
                                "action": {
                                    "_set": {"generation_style": "Gamma",
                                             "usershift": [0.0, 0.0,
                                                           0.0]}}})
            # -----------------------------------------------------------
        if "rot_matrix" in self.errors:
            # --Modified------------------------------------------------------
            if vi["KPOINTS"].style == Kpoints.supported_modes.Monkhorst or \
                    vi["KPOINTS"].kpts_shift != [0.0, 0.0, 0.0]:
                actions.append({"dict": "KPOINTS",
                                "action": {
                                    "_set": {"generation_style": "Gamma",
                                             "usershift": [0.0, 0.0,
                                                           0.0]}}})
            # -----------------------------------------------------------
            else:
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"ISYM": 0}}})

        if "amin" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"AMIN": "0.01"}}})

        if "triple_product" in self.errors:
            s = vi["POSCAR"].structure
            trans = SupercellTransformation(
                ((1, 0, 0), (0, 0, 1), (0, 1, 0)))
            new_s = trans.apply_transformation(s)
            actions.append({"dict": "POSCAR",
                            "action": {
                                "_set": {"structure": new_s.as_dict()}},
                            "transformation": trans.as_dict()})

        if "pricel" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {
                                "_set": {"SYMPREC": 1e-8, "ISYM": 0}}})

        if "brions" in self.errors:
            potim = float(vi["INCAR"].get("POTIM", 0.5)) + 0.1
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"POTIM": potim}}})

        if "zbrent" in self.errors:
            # Modified so as not to use IBRION=1 as it does not show the
            # eigenvalues in vasprun.xml >>>>>>>>>>>>
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"ADDGRID": True}}})
            actions.append({"file": "CONTCAR",
                            "action": {"_file_copy": {"dest": "POSCAR"}}})
        #            actions.append({"dict": "INCAR",
        #                            "action": {"_set": {"IBRION": 1}}})
        #            actions.append({"file": "CONTCAR",
        #                            "action": {"_file_copy": {"dest": "POSCAR"}}})
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if "too_few_bands" in self.errors:
            if "NBANDS" in vi["INCAR"]:
                nbands = int(vi["INCAR"]["NBANDS"])
            else:
                with open("OUTCAR") as f:
                    for line in f:
                        if "NBANDS" in line:
                            try:
                                d = line.split("=")
                                nbands = int(d[-1].strip())
                                break
                            except (IndexError, ValueError):
                                pass

        # >>>>>>>>>>>>
        # Modified since when the nbands is less than 9, it is not incremented
            actions.append({"dict": "INCAR",
                            "action": {
                                "_set": {"NBANDS": int(1.1 * nbands) + 2}}})
                   #            "_set": {"NBANDS": int(1.1 * nbands)}}})
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if "pssyevx" in self.errors:
            actions.append({"dict": "INCAR", "action":
                {"_set": {"ALGO": "Normal"}}})
        if "eddrmm" in self.errors:
            # RMM algorithm is not stable for this calculation
            if vi["INCAR"].get("ALGO", "Normal") in ["Fast", "VeryFast"]:
                actions.append({"dict": "INCAR", "action":
                    {"_set": {"ALGO": "Normal"}}})
            else:
                potim = float(vi["INCAR"].get("POTIM", 0.5)) / 2.0
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"POTIM": potim}}})
            if vi["INCAR"].get("ICHARG", 0) < 10:
                actions.append({"file": "CHGCAR",
                                "action": {
                                    "_file_delete": {'mode': "actual"}}})
                actions.append({"file": "WAVECAR",
                                "action": {
                                    "_file_delete": {'mode': "actual"}}})

        if "edddav" in self.errors:
            if vi["INCAR"].get("ICHARG", 0) < 10:
                actions.append({"file": "CHGCAR",
                                "action": {
                                    "_file_delete": {'mode': "actual"}}})
            actions.append({"dict": "INCAR", "action":
                {"_set": {"ALGO": "All"}}})

        if "grad_not_orth" in self.errors:
            if vi["INCAR"].get("ISMEAR", 1) < 0:
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"ISMEAR": "0"}}})

        if "zheev" in self.errors:
            if vi["INCAR"].get("ALGO", "Fast").lower() != "exact":
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"ALGO": "Exact"}}})
        if "elf_kpar" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"KPAR": 1}}})

        if "rhosyg" in self.errors:
            if vi["INCAR"].get("SYMPREC", 1e-4) == 1e-4:
                actions.append({"dict": "INCAR",
                                "action": {"_set": {"ISYM": 0}}})
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"SYMPREC": 1e-4}}})

        if "posmap" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"SYMPREC": 1e-6}}})

        if "point_group" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"ISYM": 0}}})

        ViseVaspModder(vi=vi).apply_actions(actions)
        return {"errors": list(self.errors), "actions": actions}


class ViseUnconvergedErrorHandler(orig_handlers.UnconvergedErrorHandler):

    def correct(self):
        backup(orig_handlers.VASP_BACKUP_FILES)
        v = Vasprun(self.output_filename)
        actions = [{"file":   "CONTCAR",
                    "action": {"_file_copy": {"dest": "POSCAR"}}}]
        if not v.converged_electronic:
            # For SCAN try switching to CG for the electronic minimization
            if "SCAN" in v.incar.get("METAGGA", "").upper():
                new_settings = {"ALGO": "All"}
            else:
                new_settings = {"ISTART":   1,
                                "ALGO":     "Normal",
                                "NELMDL":   -6,
                                "BMIX":     0.001,
                                "AMIX_MAG": 0.8,
                                "BMIX_MAG": 0.001}

            if all([v.incar.get(k, "") == val for k, val in
                    new_settings.items()]):
                return {"errors": ["Unconverged"], "actions": None}

            actions.append({"dict":   "INCAR",
                            "action": {"_set": new_settings}})

        # Instead of changing the IBRION, we reduce the EDIFF since we usually
        # use the constant EDIFF value as default.
        elif not v.converged_ionic:
            ediff = v.incar["EDIFF"]
            actions.append({"dict":   "INCAR",
                            "action": {"_set": {"EDIFF": ediff * 0.5}}})
            actions.append({"dict":   "INCAR",
                            "action": {"_set": {"ADDGRID": True}}})
        ViseVaspModder().apply_actions(actions)
        return {"errors": ["Unconverged"], "actions": actions}


class DielectricMaxIterationErrorHandler(ErrorHandler):
    """Detects if the SCF is not attained. """

    is_monitor = True

    def __init__(self, oszicar="OSZICAR", incar="INCAR"):
        self.error_count = Counter()
        self.oszicar = oszicar
        self.incar = incar

    def check(self):

        incar = Incar.from_file(self.incar)
        nelm = incar.get("NELM", 60)
        unconverged_line = f":{nelm:4d}"
        with open(self.oszicar, "r") as fr:
            return any(unconverged_line in line for line in fr)

    def correct(self):
        # Uncorrectable error. Just return None for actions.
        return {"errors": ["No_DFPT_convergence"], "actions": None}


class MemoryOverflowHandler(ErrorHandler):
    """Detects if the memory is overflowed. """

    is_monitor = True

    def __init__(self, memory_usage_limit=0.85):
        """
        Initializes the handler with the output file to check.

        Args:
            memory_usage_limit (float):
        """
        self.memory_usage_limit = memory_usage_limit

    def check(self):
        res = subprocess.check_output('free').split()
        memory_usage = int(res[8]) / int(res[7])
        if memory_usage > self.memory_usage_limit:
            return True

    def correct(self):
        # Uncorrectable error. Just return None for actions.
        return {"errors": ["Too_much_memory_usage"], "actions": None}


class TooLongTimeCalcErrorHandler(ErrorHandler):
    """Detects if calculation runtime is longer than the given timeout."""

    is_monitor = True

    def __init__(self, timeout=129600):
        """

        60 * 60 * 36 = 129600 (36 hours)

        Args:
            timeout (int):
        """
        self.timeout = timeout

    def check(self):
        now_time = int(subprocess.check_output(['date', '+%s']))
        incar_time = \
            int(subprocess.check_output(['date', '+%s', '-r', 'INCAR']))
        incar_age = now_time - incar_time
        if incar_age > self.timeout:
            return True

    def correct(self):
        # Uncorrectable error. Just return None for actions.
        return {"errors": ["Too_long_calc"], "actions": None}


class DivergingEnergyErrorHandler(ErrorHandler):

    def __init__(self, output_filename="OSZICAR", energy_criterion=10**6):
        """Initializes the handler with the output file to check.

        Args:
            output_filename (str): This is the OSZICAR file. Change
                this only if it is different from the default (unlikely).
        """
        self.output_filename = output_filename
        self.energy_criterion = energy_criterion

    def check(self):
        oszicar = Oszicar(self.output_filename)
        esteps = oszicar.electronic_steps
        # OSZICAR file can be empty, thus we need try-except here.
        try:
            max_energy = max([es["E"] for es in esteps[-1]])
        except (IndexError, KeyError):
            return False

        if max_energy > self.energy_criterion:
            return True

    def correct(self):
        # Uncorrectable error. Just return None for actions.
        return {"errors": ["Energy_diverging"], "actions": None}


class ReturnErrorHandler(ErrorHandler):
    """Return Error for test."""
    is_monitor = True

    def __init__(self):
        pass

    def check(self):
        return True

    def correct(self):
        # Uncorrectable error. Just return None for actions.
        return {"errors": ["Always return Error with this."], "actions": None}



