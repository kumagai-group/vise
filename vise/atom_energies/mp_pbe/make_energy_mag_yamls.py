# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os
import warnings
from pathlib import Path

from pymatgen.core import Element
from pymatgen.io.vasp import Oszicar

warnings.simplefilter('ignore')


def make_energy_mag_yamls():

    filename_e = "../mp_atom_energy.yaml"
    filename_m = "../mp_atom_mag.yaml"
    try:
        os.remove(filename_e)
    except FileNotFoundError:
        pass
    try:
        os.remove(filename_m)
    except FileNotFoundError:
        pass

    energies = []
    mags = []
    for e in Element:
        # Z(Pr) = 59, Z(Yb) = 70, Z(Pa) = 91, Z(No) = 102
        if 59 <= e.Z <= 70 or 84 <= e.Z:
            continue

        if e == Element.Ti:
            d = "algo_n/Ti_mag2"
        elif e == Element.La:
            d = "algo_n/La_mag1"
        elif str(e) in ["B", "P", "Ce", "Co", "Fe", "Ni", "Sb", "Si", "Tc", "Th"]:
            d = f"algo_n/{e}"
        else:
            d = f"algo_d/{e}"

        o = Oszicar(Path(d) / "OSZICAR")

        s = str(e)+':'
        energies.append(f"{s:<3} {o.final_energy:11.8f}")
        mags.append(f"{s:<3} {abs(o.ionic_steps[-1]['mag']):5.2f}")

    Path(filename_e).write_text("\n".join(energies))
    Path(filename_m).write_text("\n".join(mags))


if __name__ == '__main__':
    make_energy_mag_yamls()
