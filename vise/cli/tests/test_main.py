# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import unittest
import numpy as np
from unittest.mock import patch

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element
from pymatgen.analysis.phase_diagram import (
    PDEntry, PhaseDiagram, CompoundPhaseDiagram, PDPlotter)

from vise.chempotdiag.chem_pot_diag import ChemPotDiag, sort_coords
from vise.util.testing import ViseTest
from vise.cli.tests.mymodule import parse_args


class TestMain(ViseTest):

    def test(self):
        with self.assertRaises(SystemExit):
            pa = parse_args(["bg", "-u", "vasprun"])
            print(pa)