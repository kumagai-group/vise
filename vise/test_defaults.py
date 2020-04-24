# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
import tempfile
from pathlib import Path

from vise.defaults import UserSettings

# @patch("vise.config.get_user_settings",
#        ({"SYMMETRY_TOLERANCE": 10, "xc": "a"}, "b"))
# class ConfigTest(ViseTest):
#     from vise.config import \
#         SYMMETRY_TOLERANCE, MAIN_SETTINGS, VISE_YAML_FILES

    # def test_override(self):
    #     self.assertEqual(10, self.SYMMETRY_TOLERANCE)
    #     self.assertEqual("a", self.MAIN_SETTINGS["xc"])
    #     self.assertEqual("b", self.VISE_YAML_FILES[0])


"""Todo
* Method to return the file names and their respective setting keys.
* Parse hidden directory at home directory when receive the name.
* When the key includes "/", the absolute path is added as a prefix.
  E.g., unitcell/unitcell.json -> /something/../unitcell/unitcell.json
* The value of "potcar_dict: Mg_pv O_h" is a single string of
  "Mg_pv O_h", which is suited for main default value.
* Allow hidden files, e.g., .test.yaml

DONE
- Get current working directory Path.
- Find all the yaml files with the given name to the root directory.
"""


def test_get_current_working_directory():
    with tempfile.TemporaryDirectory() as dirname:
        os.chdir(dirname)
        path_dirname = Path.cwd()
        user_settings = UserSettings(yaml_filename="test.yaml")
    assert user_settings.cwd == path_dirname


def test_get_yaml_filenames():
    with tempfile.TemporaryDirectory() as dirname:

        os.chdir(dirname)
        path_second_last_dirname = Path.cwd()
        Path("test.yaml").touch()

        Path("test_dir").mkdir()
        os.chdir("test_dir")
        Path("test.yaml").touch()

        user_settings = UserSettings(yaml_filename="test.yaml")
        actual = user_settings.yaml_file_list

    path_last_dirname = path_second_last_dirname / "test_dir"

    expected = [path_second_last_dirname / "test.yaml",
                path_last_dirname / "test.yaml"]

    assert actual == expected


# def test_user_settings():
#     with tempfile.TemporaryDirectory() as dirname:

        # os.chdir(dirname)
        # with open('test.yaml', mode="w") as f:
        #     f.write("""# This file is needed for the unittest.
        # float: 0.01
        # dict:
        #   d: 5
        #   e: 3
        # f_str: g h
        # i_dir: ~/j""")

        # user_settings = UserSettings(yaml_filename="test.yaml")
        # actual = user_settings.yaml_filenames

