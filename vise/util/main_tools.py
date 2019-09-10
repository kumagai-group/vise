from typing import Optional

from pydefect.util.tools import is_str_int, is_str_digit


def potcar_str2dict(potcar_list: Optional[str]) -> dict:
    """Sanitize the string type potcar setting to dict.

    An exmaple is "Mg_pv O_h" -> {"Mg": "Mg_pv", "O": "O_h"}
    If potcar_list is None, {} is returned.

    Args:
         potcar_list (str/None)

    Return:
         Dictionary of potcar types.
    """
    if potcar_list is None:
        return {}
    elif isinstance(potcar_list, str):
        potcar_list = potcar_list.split()\

    d = {}
    for p in potcar_list:
        element = p.split("_")[0]
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
    flattened_list = flattened_list if flattened_list else []

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