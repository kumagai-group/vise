# -*- coding: utf-8 -*-

from pymatgen.io.vasp import inputs as orig_inputs
from vise.input_set.incar import ViseIncar

# Monkey patch
orig_inputs.Incar = ViseIncar
ViseVaspInput = orig_inputs.VaspInput
