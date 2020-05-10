# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os

import pytest

from vise.defaults import UserSettings


@pytest.fixture
def user_settings(tmpdir):
    second_last_dir = tmpdir
    os.chdir(second_last_dir)
    second_last_dir.join("test.yaml").write("""
        float: 0.01
        dict:
          d1: 5
          d2: 3""")

    terminal_dir = tmpdir.join("terminal_dir")
    terminal_dir.mkdir()
    os.chdir("terminal_dir")
    terminal_dir.join("test.yaml").write("""
        str: s1 s2
        dict:
          d1: 10 
          d2: 3""")

    return UserSettings("test.yaml"), second_last_dir, terminal_dir


def test_get_yaml_filenames(user_settings):
    user_settings, second_last_dir, terminal_dir = user_settings
    actual = user_settings.yaml_files_from_root_dir
    expected = [second_last_dir / "test.yaml", terminal_dir / "test.yaml"]

    assert actual == expected


def test_overridden_user_settings(user_settings):
    user_settings, _, _ = user_settings
    actual = user_settings.user_settings
    expected = {"float": 0.01,
                "dict": {"d1": 10, "d2": 3},
                "str": "s1 s2"}
    assert actual == expected


def test_parse_hidden_file(tmpdir):
    os.chdir(tmpdir)
    tmpdir.join("test.yaml").write("""
    float: 0.01""")
    user_settings = UserSettings(yaml_filename="test.yaml")
    actual = user_settings.user_settings
    expected = {"float": 0.01}
    assert actual == expected


def test_add_absolute_path(tmpdir):
    os.chdir(tmpdir)
    tmpdir.join("test.yaml").write("""
    path: a/b""")
    user_settings = UserSettings(yaml_filename="test.yaml")
    actual = user_settings.user_settings["path"]
    expected = tmpdir / "a/b"
    assert actual == expected
