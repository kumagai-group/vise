# -*- coding: utf-8 -*-
from copy import deepcopy

from pymatgen.io.vasp import Potcar

from vise.input_set.settings_potcar import XcTaskPotcar
from vise.input_set.xc import Xc

from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class XcTaskPotcarTest(ViseTest):
    def test_gga(self) -> None:
        xc_task_potcar = XcTaskPotcar.from_options(xc=Xc.pbe,
                                                   symbol_list=["Mg", "O"])

        expected = Potcar(["Mg", "O"], functional="PBE_54").as_dict()
        self.assertEqual(expected, xc_task_potcar.potcar.as_dict())
        self.assertEqual(400.0, xc_task_potcar.max_enmax)

    # def test_lda(self) -> None:
    #     xc_task_potcar = XcTaskPotcar.from_options(xc=Xc.lda,
    #                                                symbol_list=["Mg", "O"])

        # expected = Potcar(["Mg", "O"], functional="LDA").as_dict()
        # self.assertEqual(expected, xc_task_potcar.potcar.as_dict())

    def test_mp_relax_set(self) -> None:
        xc_task_potcar = \
            XcTaskPotcar.from_options(xc=Xc.pbe,
                                      symbol_list=["Mg", "O"],
                                      potcar_set_name="mp_relax_set")

        expected = Potcar(["Mg_pv", "O"], functional="PBE_54").as_dict()
        self.assertEqual(expected, xc_task_potcar.potcar.as_dict())

    def test_override(self) -> None:
        xc_task_potcar = \
            XcTaskPotcar.from_options(xc=Xc.pbe,
                                      symbol_list=["Mg", "O"],
                                      override_potcar_set={"Mg": "Mg_pv"})

        expected = Potcar(["Mg_pv", "O"], functional="PBE_54").as_dict()
        self.assertEqual(expected, xc_task_potcar.potcar.as_dict())
