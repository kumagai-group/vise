# contents of our original code file e.g. code.py
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os


def get_os_user_lower():
    """Simple retrieval function.
    Returns lowercase USER or raises OSError."""
    username = os.getenv("USER")

    if username is None:
        raise OSError("USER environment is not set.")

    return username.lower()