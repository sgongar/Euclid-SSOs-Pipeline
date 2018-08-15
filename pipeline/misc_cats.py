#!/usr/bin/python
# -*- coding: utf-8 -*-

"""


Todo:
    * Improve log messages
    * Improve usability
"""

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_catalog(ccd):
    """ returns catalog from ccd name

    :param ccd:
    :return:
    """
    cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
            ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
            ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
            ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
            ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
            ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
            ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
            ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
            ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
            ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
            ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
            ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

    ccd_position = ccd[4:9]
    dither_n = ccd[11:12]

    for cat_ in cats:
        if ccd_position == cat_[0] and dither_n == cat_[1]:
            print(ccd)
