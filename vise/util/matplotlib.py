# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from matplotlib.ticker import FuncFormatter


# need pos as the FuncFormatter arg requires 2 args.
# https://stackoverflow.com/questions/29139552/pyplot-remove-the-digits-of-zero-start-from-0-not-0-00

def my_formatter(tick_value, pos):
    """convert  0.0 to 0 in the plot.

    Examples:
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
    """
    if isinstance(tick_value, float):
        rounded_value = round(tick_value, ndigits=10)
        if rounded_value.is_integer():
            return int(rounded_value)
        else:
            return rounded_value
    else:
        return str(tick_value)


float_to_int_formatter = FuncFormatter(my_formatter)
