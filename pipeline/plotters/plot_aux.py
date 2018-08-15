#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1 Initial release

Todo:
    *

"""

import numpy as np

from decimal import Decimal
from misc import significant_l


def round_number(number):
    """
    TODO Improve description


    :param number:
    :return:
    """
    number_t = significant_l(number)

    # Rounds the float value of difference
    first_digit = '%.2E' % Decimal(number - number_t)

    if first_digit[0] == '-':
        first_digit = int(first_digit[1])
    else:
        first_digit = int(first_digit[0])

    # Redondea al alza el ultimo digio
    if first_digit > 5:
        last_digit = str(number_t)[-1]
        last_digit = int(last_digit)
        last_digit += 1
        number_t = list(str(number_t))
        number_t[-1] = last_digit
        number_t = [str(i) for i in number_t]
        number_t = ''.join(number_t)
        number_t = float(number_t)
    else:
        pass  # do nothing

    return number_t


def create_ticks(seconds):
    """
    step = stp

    :param seconds:
    :return: y_ticks
    """
    divisions = 4  #

    # Gets the major step between ticks thought the difference
    # between the maximum and the minimum value of alpha
    difference = float(max(seconds)) - float(min(seconds))
    major_stp = (difference / divisions)
    #
    major_stp = float(round_number(major_stp))
    minor_stp = (major_stp / 4)

    # Gets maximum decimal position of major step
    decimals = int(str(major_stp)[::-1].find('.'))

    # Major step list starts two times before and end two times after
    # known values
    major_stps = np.arange(round(min(seconds), decimals) - major_stp * 2,
                           round(max(seconds), decimals) + major_stp * 2,
                           major_stp)
    # Minor step list starts eight times before and fend eight times
    # after know values
    minor_stps = np.arange(round(min(seconds), decimals) - minor_stp * 8,
                           round(max(seconds), decimals) + minor_stp * 8,
                           minor_stp)

    y_ticks = {'major_t': major_stps, 'minor_t': minor_stps,
               'major_s': major_stp, 'minor_s': minor_stp}

    return y_ticks
