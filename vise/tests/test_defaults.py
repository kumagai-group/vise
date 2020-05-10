# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os

import pytest


@pytest.fixture
def modified_defaults(tmpdir):
    os.chdir(tmpdir)
    tmpdir.join("vise.yaml").write("""
    ldauu: 
        Mn: 5""")
    # Need to import here as Defaults is a singleton class.
    from vise.defaults import Defaults

    return Defaults()


@pytest.mark.skip(reason='default_override')
def test_user_settings_in_defaults(modified_defaults):
    print(modified_defaults.user_settings)
    assert modified_defaults.user_settings["ldauu"] == {"Mn": 5}
    assert modified_defaults.ldauu == {"Mn": 5}


