# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os
import warnings
from pathlib import Path
from xml.etree.ElementTree import ParseError

from pymatgen.core import Element
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from vise.atom_energies.make_atom_vasp_set import is_target_element
from vise.util.logger import get_logger


logger = get_logger(__name__)

warnings.simplefilter('ignore', UnknownPotcarWarning)


def make_energy_yaml():
    dirs = [f.name for f in os.scandir(".") if f.is_dir()]

    for e in list(Element):
        e = str(e)
        if e not in dirs or is_target_element(e) is False:
            continue

        try:
            v = Vasprun(Path(e) / "vasprun.xml")

            if v.converged_electronic is False:
                logger.warning(f"Calculation for {e} is not converged.")
                continue
        except ParseError:
            logger.warning(f"Parsing vasprun.xml for {e} failed.")
            continue

        outcar = Outcar(Path(e) / "OUTCAR")
        print(f"{e + ':':<3} {outcar.final_energy:11.8f}")


if __name__ == '__main__':
    make_energy_yaml()
