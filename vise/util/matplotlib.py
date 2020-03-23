# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from matplotlib.ticker import FuncFormatter


# need pos as the FuncFormatter arg requires 2 args.
# https://stackoverflow.com/questions/29139552/pyplot-remove-the-digits-of-zero-start-from-0-not-0-00

def my_formatter(var, pos):
    """convert  0.0 to 0 in the plot.

    Examples:
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
    """
    if isinstance(var, float):
        return int(var) if var.is_integer() else var
    return var


formatter = FuncFormatter(my_formatter)
