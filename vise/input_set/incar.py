# -*- coding: utf-8 -*-
import itertools
import os
import re
from copy import deepcopy
from math import ceil

from monty.io import zopen
from monty.serialization import loadfn
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp import Incar
from pymatgen.util.io_utils import clean_lines
from tabulate import tabulate
from vise.input_set.xc import Xc, LDA_OR_GGA, HYBRID_FUNCTIONAL
from vise.input_set.task import Task
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Potcar
from vise.input_set.sets.element_parameters import unoccupied_bands
from math import ceil

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
incar_flags_yaml = os.path.join(MODULE_DIR, "incar_flags.yaml")

# This incar_flags should be OrderedDict, but from python 3.6, dict uses
# order-preserving semantics. Furthermore, it does not affect vasp result.
incar_flags = loadfn(incar_flags_yaml)


def calc_nbands(structure: Structure, potcar: Potcar) -> int:
    """
    Calculate the total number of bands required for the unoccupied related
    properties such as optical absorption, band structure, and DOS.

    Args:
        structure (Structure/IStructure):
        potcar (Potcar):
    """

    if not all([s == p.element for s, p in zip(structure.symbol_set, potcar)]):
        raise ValueError("The structure and POTCAR file are not compatible.")

    comp = structure.composition
    nbands = sum([comp[c] * (p.nelectrons / 2 + unoccupied_bands[str(c)])
                  for c, p in zip(comp, potcar)])
    return ceil(nbands)


class ObaIncar(Incar):
    """
    Incar class modified for pretty writing of INCAR file.
    Since from_file and from_string methods in Incar class use Incar class
    constructor, we need to override them.
    """

    @staticmethod
    def from_file(filename):
        """
        Reads an Incar object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            ObaIncar object
        """
        with zopen(filename, "rt") as f:
            return ObaIncar.from_string(f.read())

    @classmethod
    def from_dict(cls, d):
        if d.get("MAGMOM") and isinstance(d["MAGMOM"][0], dict):
            d["MAGMOM"] = [Magmom.from_dict(m) for m in d["MAGMOM"]]

        return cls({k: v for k, v in d.items() if k not in ("@module",
                                                            "@class")})

    @staticmethod
    def from_string(string: str):
        """ Reads an Incar object from a string.

        This method is different from that of Incar superclass at it does not
        support ";" semantic wich split the incar flaggs.

        Args:
            string (str): Incar string

        Returns:
            ObaIncar object
        """
        lines = list(clean_lines(string.splitlines()))
        params = {}
        for line in lines:
            # YK: Support the split of ";" semantic in INCAR file
            for sub_line in line.split(";"):
                m = re.match(r'(\w+)\s*=\s*(.*)', sub_line)
                if m:
                    key = m.group(1).strip()
                    val = m.group(2).strip()
                    val = ObaIncar.proc_val(key, val)
                    params[key] = val

        return ObaIncar(params)

    def __add__(self, other):
        """
        Add all the values of another INCAR object to this object.
        Facilitates the use of "standard" INCARs.
        """
        params = {k: v for k, v in self.items()}
        for k, v in other.items():
            if k in self and v != self[k]:
                raise ValueError("Incars have conflicting values!")
            else:
                params[k] = v
        return ObaIncar(params)

    def get_string(self, sort_keys=False, pretty=False):
        """ This method is overridden for the pretty printing. """
        lines = []
        incar_keys = deepcopy(self)
        for key, val in incar_flags.items():
            comment = False
            blank_line = False
            ll = []

            for v in val:
                if v in incar_keys:
                    if comment is False:
                        lines.append(f"# {key} \n")
                        comment = True
                    if v == "MAGMOM" and isinstance(self[v], list):
                        value = []

                        if (isinstance(self[v][0], list) or
                            isinstance(self[v][0], Magmom)) and \
                                (self.get("LSORBIT") or
                                 self.get("LNONCOLLINEAR")):
                            value.append(
                                " ".join(str(i) for j in self[v] for i in j))
                        elif self.get("LSORBIT") or self.get("LNONCOLLINEAR"):
                            for m, g in itertools.groupby(self[v]):
                                value.append("3*{}*{}".format(len(tuple(g)), m))
                        # YK: Add this
                        elif len(self[v]) == 1:
                            value.append(self[v][0])
                        else:
                            # float() to ensure backwards compatibility between
                            # float magmoms and Magmom objects
                            try:
                                for m, g in itertools.groupby(self[v],
                                                              lambda x: float(
                                                                  x)):
                                    value.append(
                                        "{}*{}".format(len(tuple(g)), m))
                            except ValueError:
                                raise ValueError("MAGMOM could be improper.")

                        ll.append([v, " ".join(value)])
                    elif isinstance(self[v], list):
                        ll.append([v, " ".join([str(i) for i in self[v]])])
                    else:
                        ll.append([v, self[v]])
                    blank_line = True
                    incar_keys.pop(v)
            if blank_line:
                lines.append(str(tabulate([[l[0], "=", l[1]] for l in ll],
                                          tablefmt="plain")))
                lines.append("\n\n")

        for mson_key in ["@module", "@class"]:
            if mson_key in incar_keys:
                incar_keys.pop(mson_key)

        if len(incar_keys) > 0:
            raise ValueError(
                "{} are not valid in INCAR.".format(incar_keys.keys()))

        return "".join(lines)


def make_incar_setting(structure: Structure,
                       xc: Xc,
                       task: Task,
                       rough: bool,
                       encut: float,
                       vbm_cbm: list,
                       potcar: Potcar,
                       only_even_kpt: bool,
                       ionic_contribution: bool,
                       factor: int,
                       encut_factor_str_opt: float):
    """Make incar settings from input parameters.

    Args:
        structure (Structure):
            Input structure.
        xc (Xc):
            Exchange-correlation treatment.
        task (Task):
            Task
        rough (bool):
            Whether to use rough criterion for convergence.
            Set EDIFF = 1e-5, EDIFFG = -0.04, NKRED = 2 in INCAR.
        encut (float):
            Cutoff energy in eV.
        vbm_cbm (list):
            The valence band maximum (vbm) and conduction band minimum (cbm)
            with [vbm, cbm] format. Used for determining the dos and
            absorption spectra region.
        potcar (Potcar):
            Potcar object
        only_even_kpt (bool):
            
        ionic_contribution: (bool):
        factor (int):
        encut_factor_str_opt (float):


    """

    # Incar
    incar_settings = {}
    # - xc
    if xc is Xc.pbesol:
        incar_settings.update({"GGA": "PS"})
    elif xc is Xc.scan:
        incar_settings.update({"METAGGA": "SCAN"})
    elif xc in (Xc.pbe0, Xc.hse):
        incar_settings.update({"LHFCALC": True,
                               "PRECFOCK": "Fast",
                               "ALGO": "D",
                               "AEXX": 0.25,
                               "TIME": 0.5,
                               "PREC": "A",
                               "LWAVE": True})

    if xc is Xc.hse:
        incar_settings.update({"HFSCREEN": 0.208})

    # - rough
    if rough:
        # POTIM is 1/2.5 of the default, which allows us to relax the atoms
        # slowly, and useful for prototype calculations.
        incar_settings.update({"EDIFF": 1e-4,
                               "EDIFFG": -0.2,
                               "POTIM": 0.2})
        # This is fine as long as the num of kpoints are even numbers.
        if xc in HYBRID_FUNCTIONAL and only_even_kpt:
            incar_settings.update({"NKRED": 2})

    # - algo
    # To attain the convergence of the unoccupied bands, ALGO=N is needed.
    # This change is not switched on if ALGO is already set above.
    if task in (Task.band, Task.dos, Task.dielectric,
                Task.dielectric_function) and "ALGO" not in incar_settings:
        incar_settings.update({"ALGO": "N"})
    elif task in (Task.gw_pre_calc2,):
        incar_settings.update({"ALGO": "Exact"})

    # - lreal & prec
    # When the number of atoms is less than 40, LREAL = False
    if len(structure) <= 40:
        incar_settings.update({"LREAL": False})

    # - encut
    if task in (Task.structure_opt,):
        enmax = round(max([p.enmax for p in potcar]) * encut_factor_str_opt, 4)
        incar_settings.update({"ENCUT": enmax})

    # Always overwrite ENCUT when explicitly given.
    if encut:
        incar_settings.update({"ENCUT": encut})

    # - lorbit
    if task in (Task.band, Task.dos, Task.defect):
        incar_settings.update({"LORBIT": 12})

    # - nbands
    if task in (Task.band, Task.dos, Task.dielectric_function):
        nbands = calc_nbands(structure, potcar)
        incar_settings.update({"NBANDS": nbands})

    # - emin, emax, nedos
    if task in (Task.dos, Task.dielectric_function):
        step_size = 0.01
        if vbm_cbm:
            emin = ceil(vbm_cbm[0]) - 15 - step_size
            emax = ceil(vbm_cbm[1]) + 15
        else:
            emin = -20 - step_size
            emax = 20
        nedos = int(round((emax - emin) / step_size, 0)) + 1
        incar_settings.update({"EMIN": emin, "EMAX": emax, "NEDOS": nedos})

    # - ediff
    if task in (Task.defect,):
        incar_settings.update({"EDIFFG": -0.04})

    # - ispin
    if task in (Task.defect,):
        incar_settings.update({"ISPIN": 2})

    # - ismear
    # Note1: Task.dielectric assumes insulators.
    # Note2: band_gap belongs to <class 'numpy.bool_'>, not build-in bool.
    # Note3: Tetrahedron method is valid when num_kpoints >= 4.
    #        Thus, num_kpoints will be checked later.
    # ISMEAR = -5 is now switched off following the vasp warning below
    # | ALGO = A and IALGO = 5X tend to fail with the tetrahedron method |
    # | (e.g.Bloechls method ISMEAR=-5 is not variational) |
    # | please switch to IMSEAR = 0 - n, except for DOS calculations |
    # if (band_gap and band_gap > BAND_GAP_EXISTENCE_CRITERION) or \
    #         (vbm_cbm and
    #          vbm_cbm[1] - vbm_cbm[0] > BAND_GAP_EXISTENCE_CRITERION) or \
    #         task in (Task.dielectric, ):
    #     incar_settings.update({"ISMEAR": -5})
    # - nsw
    # NSW = 1 (not NSW = 0) is a must for dielectric constant.
    if task in (Task.single_point, Task.band, Task.dos, Task.dielectric,
                Task.dielectric_function, Task.gw_pre_calc2) or \
            xc in (Xc.gw0,):
        incar_settings.update({"NSW": 1})

    # Basically INCAR settings for dielectric constants are the same for
    # SCAN and HSE due to the VASP implementation.

    # - lepsilon
    if task in (Task.dielectric,) and xc in LDA_OR_GGA:
        incar_settings.update({"LEPSILON": True})

    # - lrpa
    if task in (Task.dielectric,) and xc in BEYOND_GGA:
        incar_settings.update({"LRPA": False})

    # -lcalceps
    if task in (Task.dielectric,) and xc in BEYOND_GGA:
        incar_settings.update({"LCALCEPS": True})

    # - loptics
    # Note: CSHIFT doesn't affect the results when using the modified vasp.
    if task in (Task.dielectric_function,):
        incar_settings.update({"LOPTICS": True})

    # - isif
    if task in (Task.single_point, Task.band, Task.dos, Task.dielectric,
                Task.dielectric_function):
        incar_settings.update({"ISIF": 0})
    elif task in (Task.defect,):
        incar_settings.update({"ISIF": 2})

    # - ibrion
    # The hybrid functional is not allowed for ionic contribution calc.
    if xc in BEYOND_GGA and task is Task.dielectric and ionic_contribution:
        raise ValueError("The {} functional/approximation is not allowed "
                         "for calculating ionic contribution of "
                         "dielectric constant.".format(str(xc)))

    if task is Task.dielectric and ionic_contribution:
        incar_settings.update({"IBRION": 8})

    # - nkred
    if not factor:
        factor = 1
        if task in (Task.dos, Task.dielectric):
            factor = 2
        elif task in (Task.dielectric_function,):
            factor = 3
        if xc in HYBRID_FUNCTIONAL and factor > 1:
            incar_settings.update({"NKRED": factor})

    return incar_settings, factor


def make_gw_incar_setting(structure: Structure,
                          xc: Xc,
                          task: Task,
                          rough: bool,
                          encut: float,
                          vbm_cbm: list,
                          potcar: Potcar,
                          only_even_kpt: bool,
                          ionic_contribution: bool,
                          factor: int,
                          encut_factor_str_opt: float):
    """Make incar settings from input parameters.

    Args:
        structure (Structure):
            Input structure.
        xc (Xc):
            Exchange-correlation treatment.
        task (Task):
            Task
        rough (bool):
            Whether to use rough criterion for convergence.
            Set EDIFF = 1e-5, EDIFFG = -0.04, NKRED = 2 in INCAR.
        encut (float):
            Cutoff energy in eV.
        vbm_cbm (list):
            The valence band maximum (vbm) and conduction band minimum (cbm)
            with [vbm, cbm] format. Used for determining the dos and
            absorption spectra region.
        potcar (Potcar):
            Potcar object
        only_even_kpt (bool):

        ionic_contribution: (bool):
        factor (int):
        encut_factor_str_opt (float):


    """

    # Incar
    incar_settings = {}
    # - xc
    if xc is Xc.pbesol:
        incar_settings.update({"GGA": "PS"})
    elif xc is Xc.scan:
        incar_settings.update({"METAGGA": "SCAN"})
    elif xc in (Xc.pbe0, Xc.hse):
        incar_settings.update({"LHFCALC": True,
                               "PRECFOCK": "Fast",
                               "ALGO": "D",
                               "AEXX": 0.25,
                               "TIME": 0.5,
                               "PREC": "A",
                               "LWAVE": True})
    elif xc in (Xc.gw0,):
        incar_settings.update({"NELM": 4,
                               "ALGO": "GW0",
                               "PRECFOCK": "Fast",
                               "NOMEGA": 100,
                               "OMEGAMAX": 30,
                               "ANTIRES": 1,
                               "NMAXFOCKAE": 2,
                               "LPEAD": False,
                               "LOPTICS": False})

    if xc is Xc.hse:
        incar_settings.update({"HFSCREEN": 0.208})

    # - rough
    if rough:
        # POTIM is 1/2.5 of the default, which allows us to relax the atoms
        # slowly, and useful for prototype calculations.
        incar_settings.update({"EDIFF": 1e-4,
                               "EDIFFG": -0.2,
                               "POTIM": 0.2})
        # This is fine as long as the num of kpoints are even numbers.
        if xc in HYBRID_FUNCTIONAL and only_even_kpt:
            incar_settings.update({"NKRED": 2})

    # - algo
    # To attain the convergence of the unoccupied bands, ALGO=N is needed.
    # This change is not switched on if ALGO is already set above.
    if task in (Task.band, Task.dos, Task.dielectric,
                Task.dielectric_function) and "ALGO" not in incar_settings:
        incar_settings.update({"ALGO": "N"})
    elif task in (Task.gw_pre_calc2,):
        incar_settings.update({"ALGO": "Exact"})

    # - lreal & prec
    # When the number of atoms is less than 40, LREAL = False
    if len(structure) <= 40:
        incar_settings.update({"LREAL": False})
    if task in (Task.gw_pre_calc1, Task.gw_pre_calc2):
        incar_settings.update({"LREAL": False})
        incar_settings.update({"PREC": "A"})

    # - encut
    if task in (Task.structure_opt, Task.gw_pre_calc1):
        enmax = round(max([p.enmax for p in potcar]) * encut_factor_str_opt, 4)
        incar_settings.update({"ENCUT": enmax})

    # Always overwrite ENCUT when explicitly given.
    if encut:
        incar_settings.update({"ENCUT": encut})

    # - lorbit
    if task in (Task.band, Task.dos, Task.defect):
        incar_settings.update({"LORBIT": 12})

    # - nbands
    # For Task.gw_pre_calc2 and xc.gw,
    # NBANDS = (the number of unoccupied bands approximately).
    # TODO: Check if this is fine for gw_pre_calc2
    if task in (Task.band, Task.dos, Task.dielectric_function) and \
            xc not in (Xc.gw0, ):
        nbands = calc_nbands(structure, potcar)
        incar_settings.update({"NBANDS": nbands})
    if xc in (Xc.gw0,):
        nbandsgw = calc_nbands(structure, potcar)
        incar_settings.update({"NBANDSGW": nbandsgw})

    # - emin, emax, nedos
    if task in (Task.dos, Task.dielectric_function):
        step_size = 0.01
        if vbm_cbm:
            emin = ceil(vbm_cbm[0]) - 15 - step_size
            emax = ceil(vbm_cbm[1]) + 15
        else:
            emin = -20 - step_size
            emax = 20
        nedos = int(round((emax - emin) / step_size, 0)) + 1
        incar_settings.update({"EMIN": emin, "EMAX": emax, "NEDOS": nedos})

    # - nelem & nelemin
    if task in (Task.gw_pre_calc2,):
        incar_settings.update({"NELM": 1})
        incar_settings.update({"NELMIN": 1})

    # - ediff
    if task in (Task.gw_pre_calc1, Task.gw_pre_calc2):
        incar_settings.update({"EDIFF": 1e-8})

    # - ediff
    if task in (Task.defect,):
        incar_settings.update({"EDIFFG": -0.04})

    # - ispin
    if task in (Task.defect,):
        incar_settings.update({"ISPIN": 2})

    # - lcharg
    if task in (Task.gw_pre_calc1,):
        if xc not in HYBRID_FUNCTIONAL:
            incar_settings.update({"LCHARG": True})
    elif task in (Task.gw_pre_calc2,):
        incar_settings.update({"LCHARG": False})

    # - icharg
    if task in (Task.gw_pre_calc2,) and xc in DFT_FUNCTIONAL:
        incar_settings.update({"ICHARG": 11})

    # - lwave
    if task in (Task.gw_pre_calc2,):
        incar_settings.update({"LWAVE": True})

    # - ismear
    # Note1: Task.dielectric assumes insulators.
    # Note2: band_gap belongs to <class 'numpy.bool_'>, not build-in bool.
    # Note3: Tetrahedron method is valid when num_kpoints >= 4.
    #        Thus, num_kpoints will be checked later.
    # ISMEAR = -5 is now switched off following the vasp warning below
    # | ALGO = A and IALGO = 5X tend to fail with the tetrahedron method |
    # | (e.g.Bloechls method ISMEAR=-5 is not variational) |
    # | please switch to IMSEAR = 0 - n, except for DOS calculations |
    # if (band_gap and band_gap > BAND_GAP_EXISTENCE_CRITERION) or \
    #         (vbm_cbm and
    #          vbm_cbm[1] - vbm_cbm[0] > BAND_GAP_EXISTENCE_CRITERION) or \
    #         task in (Task.dielectric, ):
    #     incar_settings.update({"ISMEAR": -5})
    # - nsw
    # NSW = 1 (not NSW = 0) is a must for dielectric constant.
    if task in (Task.single_point, Task.band, Task.dos, Task.dielectric,
                Task.dielectric_function, Task.gw_pre_calc2) or \
            xc in (Xc.gw0,):
        incar_settings.update({"NSW": 1})

    # Basically INCAR settings for dielectric constants are the same for
    # SCAN and HSE due to the VASP implementation.

    # - lepsilon
    if task in (Task.dielectric,) and xc in LDA_OR_GGA:
        incar_settings.update({"LEPSILON": True})

    # - lrpa
    if task in (Task.dielectric,) and xc in BEYOND_GGA:
        incar_settings.update({"LRPA": False})

    # -lcalceps
    if task in (Task.dielectric,) and xc in BEYOND_GGA:
        incar_settings.update({"LCALCEPS": True})

    # - loptics
    # Note: CSHIFT doesn't affect the results when using the modified vasp.
    if task in (Task.dielectric_function, Task.gw_pre_calc2):
        incar_settings.update({"LOPTICS": True})

    # - isif
    if task in (Task.single_point, Task.band, Task.dos, Task.dielectric,
                Task.dielectric_function, Task.gw_pre_calc2):
        incar_settings.update({"ISIF": 0})
    elif task in (Task.defect,):
        incar_settings.update({"ISIF": 2})

    # - ibrion
    # The hybrid functional is not allowed for ionic contribution calc.
    if xc in BEYOND_GGA and task is Task.dielectric and ionic_contribution:
        raise ValueError("The {} functional/approximation is not allowed "
                         "for calculating ionic contribution of "
                         "dielectric constant.".format(str(xc)))

    if task is Task.dielectric and ionic_contribution:
        incar_settings.update({"IBRION": 8})

    # - nkred
    if not factor:
        factor = 1
        if task in (Task.dos, Task.dielectric):
            factor = 2
        elif task in (Task.dielectric_function,):
            factor = 3
        if xc in HYBRID_FUNCTIONAL and factor > 1:
            incar_settings.update({"NKRED": factor})

    return incar_settings, factor
