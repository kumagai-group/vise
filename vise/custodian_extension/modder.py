# -*- coding: utf-8 -*-

from custodian.vasp.interpreter import VaspModder

from vise.input_set.vasp_input import ViseVaspInput


class ViseVaspModder(VaspModder):
    def __init__(self, actions=None, strict=True, vi=None):
        """ Initializes a Modder for VaspInput sets

        Args:
            actions ([Action]):
            strict (bool)
            vi (vise.input_set.vasp_input.ViseVaspInput)
                See docstrings of constructor of VaspModder
        """
        if not vi:
            vi = ViseVaspInput.from_directory(".")
        super().__init__(actions, strict, vi)

