# -*- coding: utf-8 -*-
import os

from pymatgen.io.vasp import Poscar, Outcar, Vasprun

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def check_vasp_output(directory_path, contcar_name="CONTCAR",
                      outcar_name="OUTCAR", vasprun_name="vasprun.xml"):
    """
    """
    p = {"contcar": True, "outcar": True, "vasprun": True, "all": None}

    try:
        Poscar.from_file(os.path.join(directory_path, contcar_name))
    except IOError:
        p["contcar"] = None
    except:
        p["contcar"] = False

    try:
        Outcar(os.path.join(directory_path, outcar_name))
    except IOError:
        p["outcar"] = None
    except:
        p["outcar"] = False

    try:
        Vasprun(os.path.join(directory_path, vasprun_name))
    except IOError:
        p["vasprun"] = None
    except:
        p["vasprun"] = False

    if p["contcar"] and p["outcar"] and p["vasprun"]:
        p["all"] = True

    return p


def vasp_convergence_ionic(directory_path, vasprun_name="vasprun.xml"):
    """ """
    vasprun = Vasprun(os.path.join(directory_path, vasprun_name))
    return vasprun.converged_ionic


def vasp_convergence_electronic(directory_path, vasprun_name="vasprun.xml"):
    """ """
    vasprun = Vasprun(os.path.join(directory_path, vasprun_name))
    return vasprun.converged_electronic


class VaspNotConvergedError(Exception):
    pass
