# -*- coding: utf-8 -*-

from inspect import signature, _empty
from pathlib import Path
import re
from typing import Optional, Callable, List, Union, Tuple

import yaml

from pymatgen.core.periodic_table import Element

from vise.util.tools import is_str_int, is_str_digit, str2bool
from distutils.util import strtobool


def potcar_str2dict(potcar_list: Union[str, List[str], None]) -> dict:
    """Sanitize potcar names to dict.

    If potcar_list is None, {} is returned.
    Examples
        ["Mg_pv", O_h"] -> {"Mg": "Mg_pv", "O": "O_h"}
        "Mg_pv" -> {"Mg": "Mg_pv"}
        None -> {}

    Args:
         potcar_list (str/list/None):
            List of Potcar names or a name, including element names before the
            hyphen.

    Returns:
         Dictionary with element names as keys and potcar names as values.
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
        try:
            Element(element)
        except ValueError:
            raise ValueError(f"First part of POTCAR file name {p} is not an "
                             f"element supported yet.")

        d[element] = p
    return d


def list2dict(flattened_list: Optional[list], key_candidates: list) -> dict:
    """Sanitize a list to a dict with keys in the flags

    If a string in list does not exist in key_candidates, raise ValueError.

    Examples:
        flattened_list = ["ENCUT", "500", "MAGMOM", "4", "4", "LWAVE", "F"]
        key_candidates = ["ENCUT", "MAGMOM", "LWAVE", ...]
        list2dict(flattened_list, key_candidates) =
                            {"ENCUT": 500, "MAGMOM": [4, 4], "LWAVE": False}

        flattened_list = ["ENCUT", "500", "MAGMAM", "4", "4"]
        list2dict(flattened_list, key_candidates) raises ValueError

    Args:
        flattened_list (list):
            Input list.
        key_candidates (list):
            List of key candidates, e.g., INCAR flags.

    Returns:
        Sanitized dict.
    """
    flattened_list = flattened_list or []

    d = {}
    key = None  # key is None at the beginning or after inserting a value.
    value_list = []

    def insert():
        if not value_list:
            raise ValueError(f"Invalid input: {flattened_list}.")
        d[key] = value_list[0] if len(value_list) == 1 else value_list

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
            # We first need to check numbers as 0 and 1 are parse as False and
            # True with str2bool.
            if is_str_int(string):
                value = int(string)
            else:
                try:
                    value = str2bool(string)
                except ValueError:
                    if is_str_digit(string):
                        value = float(string)
                    else:
                        value = string

            value_list.append(value)
    else:
        if key:
            insert()

    return d


def dict2list(d: dict) -> list:
    """Flatten dict to list w/ string keys. String is also separated by space.

    Example:
        dict2list({"a": 1, "b": "2 3 4", "c": True}) =
                                 ["a", "1", "b", "2", "3", "4", "c", "True"]

    Args:
         d (dict):
            Input dict.

    Returns:
         List of flattened dict.
    """

    d = d or {}
    flattened_list = []
    for k, v in d.items():
        flattened_list.append(k)
        if isinstance(v, str):
            flattened_list.extend(v.split())
        else:
            flattened_list.append(str(v))

    return flattened_list


def get_user_settings(yaml_filename: str,
                      setting_keys: list) -> Tuple[dict, Optional[Path]]:
    """Get the user specifying settings written in yaml_filename

    Note1: The yaml_filename is explored in the parent folders up to home
           or root directory until it's found. If it does not exist, empty
           dictionary is returned.
    Note2: When the key includes "/", the absolute path is added as a prefix.
           E.g., unitcell/unitcell.json -> /something/../unitcell/unitcell.json
    Note3: The value of "potcar_set: Mg_pv O_h" is string of  "Mg_pv O_h", which
           is suited for main default value.

    Args:
        yaml_filename (str):
            User setting yaml filename.
        setting_keys (list):
            Only setting_keys are valid as input keys, otherwise raise
            ValueError.

    Returns:
        Tuple of dictionary of configs and Path to the config file.
    """

    config_path = Path.cwd()
    home = Path.home()

    while True:
        if config_path == home or config_path == Path("/"):
            return {}, None

        f = config_path / yaml_filename
        if f.exists():
            with open(f, "r") as fin:
                user_settings = yaml.load(fin, Loader=yaml.FullLoader)
            break

        else:
            config_path = config_path.parent

    user_settings = user_settings or {}

    # Add full path
    for k, v in user_settings.items():
        if k not in setting_keys:
            raise ValueError(f"Key {k} in {yaml_filename} is invalid."
                             f"The candidate keys are {setting_keys}")
        if isinstance(v, str) and re.match(r'\S*/\S*', v):
            user_settings[k] = str(config_path / v)

    return user_settings, f


def get_default_args(function: Callable) -> dict:
    """Get the default values of the arguments in the method/function.

    inspect._empty means no default.

    Args:
        function (Callable):
            Method or function. when class is inserted, cls.__init__ is called.

    Returns:
        Dict of default values.
    """
    defaults = {}
    signature_obj = signature(function)
    for name, param in signature_obj.parameters.items():
        if param.default != _empty:
            defaults[name] = param.default

    return defaults
