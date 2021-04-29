# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
import re
from copy import copy


def latexify(formula):
    """
    Copied from pymatgen.2021.3.9 as it will be removed from 2022.

    Generates a LaTeX formatted formula. E.g., Fe2O3 is transformed to
    Fe$_{2}$O$_{3}$.

    Args:
        formula (str): Input formula.

    Returns:
        Formula suitable for display as in LaTeX with proper subscripts.
    """
    return re.sub(r"([A-Za-z\(\)])([\d\.]+)", r"\1$_{\2}$", formula)


def numbers_to_lowercases(string: str):
    result = copy(string)
    for i, j in zip(range(10),
                    ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"]):
        result = result.replace(str(i), j)
    return result
