# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.


import pytest

from vise import __version__
from vise.input_set.task import Task
from vise.input_set.vise_log import ViseLog
from vise.input_set.xc import Xc
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def vise_log():
    return ViseLog(version=__version__, task=Task.structure_opt, xc=Xc.pbe,
                   input_options={"charge": 1}, user_incar_settings={"NSW": 2})


def test_vise_log(vise_log, tmpdir):
    expected_text = """input_options:
  charge: 1
task: structure_opt
user_incar_settings:
  NSW: 2
version: 0.6.3
xc: pbe
"""
    assert_yaml_roundtrip(vise_log, tmpdir, expected_text)


def test_vise_log_no_incar_settings(tmpdir):
    vise_log_no_incar_settings = ViseLog(version=__version__,
                                         task=Task.structure_opt, xc=Xc.pbe,
                                         input_options={"charge": 1})
    expected_text = """input_options:
  charge: 1
task: structure_opt
version: 0.6.3
xc: pbe
"""
    assert_yaml_roundtrip(vise_log_no_incar_settings, tmpdir, expected_text)


