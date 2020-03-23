# -*- coding: utf-8 -*-

from custodian.custodian import Validator
from pymatgen.io.vasp import Vasprun
from collections import deque

import os
import logging


class VasprunXMLValidator(Validator):
    """
    Checks that a valid vasprun.xml was generated
    """

    def __init__(self,
                 vasprun_xml_file="vasprun.xml.finish",
                 output_file="vasp.out.finish",
                 stderr_file="std_err.txt"):
        """
        Args:
            output_file (str): Name of file VASP standard output is directed to.
                Defaults to "vasp.out".
            stderr_file (str): Name of file VASP standard error is direct to.
                Defaults to "std_err.txt".
        """
        self.vasprun_xml_file = vasprun_xml_file
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.logger = logging.getLogger(self.__class__.__name__)

    def check(self):
        try:
            Vasprun(self.vasprun_xml_file)
        except Exception:
            exception_context = {}

            if os.path.exists(self.output_file):
                with open(self.output_file, "r") as output_file:
                    output_file_tail = deque(output_file, maxlen=10)
                exception_context["output_file_tail"] = "".join(output_file_tail)

            if os.path.exists(self.stderr_file):
                with open(self.stderr_file, "r") as stderr_file:
                    stderr_file_tail = deque(stderr_file, maxlen=10)
                exception_context["stderr_file_tail"] = "".join(stderr_file_tail)

            if os.path.exists(self.vasprun_xml_file):
                stat = os.stat("vasprun.xml")
                exception_context["vasprun_st_size"] = stat.st_size
                exception_context["vasprun_st_atime"] = stat.st_atime
                exception_context["vasprun_st_mtime"] = stat.st_mtime
                exception_context["vasprun_st_ctime"] = stat.st_ctime

                with open("vasprun.xml", "r") as vasprun:
                    vasprun_tail = deque(vasprun, maxlen=10)
                exception_context["vasprun_tail"] = "".join(vasprun_tail)

            self.logger.error("Failed to load vasprun.xml",
                              exc_info=True, extra=exception_context)
            return True

        return False


class VaspFilesValidator(Validator):
    """
    Check for existence of some of the files that VASP
        normally create upon running.
    """

    def __init__(self):
        pass

    def check(self):
        for vfile in ["CONTCAR.finish", "OUTCAR.finish"]:
            if not os.path.exists(vfile):
                return True
        return False
