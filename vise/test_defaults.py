# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
import tempfile



"""
TODO
* Parse vise.yaml file and create related attributes.
DONE

"""


def test_user_settings_in_defaults():
    with tempfile.TemporaryDirectory() as dirname:

        os.chdir(dirname)
        with open('vise.yaml', mode="w") as f:
            f.write("""
        ldauu: 
            Mn: 5""")
        # Need to import here as Defaults is a singleton class.
        from vise.defaults import Defaults
        defaults = Defaults()

    assert defaults.user_settings["ldauu"] == {"Mn": 5}
    assert defaults.ldauu == {"Mn": 5}


