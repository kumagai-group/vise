# -*- coding: utf-8 -*-

import re
from pathlib import Path
from typing import Optional

import yaml

from pydefect.util.tools import is_str_int, is_str_digit


def potcar_str2dict(potcar_list: Optional[str]) -> dict:
    """Sanitize the string type potcar setting to dict.

    An example is "Mg_pv O_h" -> {"Mg": "Mg_pv", "O": "O_h"}
    If potcar_list is None, {} is returned.

    Args:
         potcar_list (str/None):

    Returns:
         Dictionary of potcar types.
    """
    if potcar_list is None:
        return {}
    elif isinstance(potcar_list, str):
        potcar_list = [potcar_list]

    d = {}
    for p in potcar_list:
        element = p.split("_")[0]
        if element in d:
            raise ValueError("Multiple POTCAR files for an element are not "
                             "supported yet.")
        d[element] = p
    return d


def list2dict(flattened_list: Optional[list], key_candidates: list) -> dict:
    """Sanitize the list to dict with keys in the flags

    If a string in l does not exist in key_candidates, raise ValueError.

    key_candidates = ["ENCUT", "MAGMOM", "LWAVE", ...]
    list2dict(["ENCUT", "500", "MAGMOM", "4", "4", "LWAVE", "F"]) =
                        {"ENCUT": 500, "MAGMOM": [4, 4], "LWAVE": False}

    arg_list = ["ENCUT", "500", "MAGMAM", "4", "4"]
    raise ValueError

    Args:
        flattened_list (list): Input list
        key_candidates (list): List of key candidates, e.g., INCAR flags.
    Return:
        Sanitized dict
    """
    flattened_list = flattened_list or []

    d = {}
    key = None
    value_list = []

    def insert():
        if not value_list:
            raise ValueError(f"Invalid input: {flattened_list}.")
        if len(value_list) == 1:
            d[key] = value_list[0]
        else:
            d[key] = value_list

    for string in flattened_list:
        if key is None and string not in key_candidates:
            raise ValueError(f"Keys are invalid: {flattened_list}.")
        elif string in key_candidates:
            if key:
                insert()
                key = None
                value_list = []
            key = string
        else:
            if string.lower() == "true" or string == "T":
                value = True
            elif string.lower() == "false" or string == "F":
                value = False
            elif is_str_int(string):
                value = int(string)
            elif is_str_digit(string):
                value = float(string)
            else:
                value = string

            value_list.append(value)
    else:
        if key:
            insert()

    return d


def get_user_settings(yaml_filename: str,
                      setting_keys: list) -> dict:
    """Get the user specifying settings written in yaml_filename

    Note1: The yaml_filename is explored in the parent folders up to home
           or root directory until it's found. If it does not exist, empty
           dictionary is returned.
    Note2: When the key includes "/", the absolute path is added as a prefix.
           E.g., unitcell/unitcell.json -> /something/../unitcell/unitcell.json
    Note3: The value of "potcar_set: Mg_pv O_h" is "Mg_pv O_h" string, which
           is suited when used for main default value.

    Args:
        yaml_filename (str):
            User setting yaml filename.
        setting_keys (list):
            Only setting_keys are valid as input keys, otherwise raise
            ValueError.

    Returns:
        Dictionary of configs.
    """

    config_path = Path.cwd()
    home = Path.home()

    while True:
        if config_path == home or config_path == Path("/"):
            return {}

        f = config_path / yaml_filename
        if f.exists():
            with open(f, "r") as f:
                user_settings = yaml.load(f, Loader=yaml.FullLoader)
            break

        else:
            config_path = config_path.parent

    # Add full path
    for k, v in user_settings.items():
        if k not in setting_keys:
            raise ValueError(f"Key {k} in {yaml_filename} is invalid."
                             f"The candidate keys are {setting_keys}")
        if isinstance(v, str) and re.match(r'\S*/\S*', v):
            user_settings[k] = str(config_path / v)

    return user_settings


def dict2list(d: dict) -> list:
    """Sanitize the string type potcar setting to dict.

    The string is also separated by space. An example is
    dict2list({"a": 1, "b": "2 3 4", "c": True}) =
                                 ["a", "1", "b", "2", "3", "4", "c", "True"]

    Args:
         d (dict)

    Return:
         list of flattened dict
    """

    d = d if d else {}
    flattened_list = []
    for k, v in d.items():
        flattened_list.append(k)
        if isinstance(v, str):
            flattened_list.extend(v.split())
        else:
            flattened_list.append(str(v))

    return flattened_list
