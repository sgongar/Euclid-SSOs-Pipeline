#!/usr/bin/python
# -*- coding: utf-8 -*-

"""


Todo:
    * Improve log messages
    * Improve usability
"""

from astropy.io import fits
import numpy as np

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_position(order):
    """

    :param order:
    :return:
    """
    order_d = {0: [1, 6], 1: [1, 5], 2: [1, 4], 3: [1, 3], 4: [1, 2], 5: [1, 1],
               6: [2, 6], 7: [2, 5], 8: [2, 4], 9: [2, 3], 10: [2, 2],
               11: [2, 1], 12: [3, 6], 13: [3, 5], 14: [3, 4], 15: [3, 3],
               16: [3, 2], 17: [3, 1], 18: [4, 6], 19: [4, 5], 20: [4, 4],
               21: [4, 3], 22: [4, 2], 23: [4, 1], 24: [5, 6], 25: [5, 5],
               26: [5, 4], 27: [5, 3], 28: [5, 2], 29: [5, 1], 30: [6, 6],
               31: [6, 5], 32: [6, 4], 33: [6, 3], 34: [6, 2], 35: [6, 1]}

    coords = 'x{}_y{}'.format(order_d[order][0], order_d[order][1])

    return coords


def create_ccds(logger, proc, fits_dir, fpa_dir, fpa_file):
    """

    :param logger:
    :param proc:
    :param fits_dir:
    :param fpa_dir:
    :param fpa_file:
    :return:
    """
    quadrants_d = {}
    dither = fpa_file[-6:-5]

    images_idxs = np.arange(1, 144, 4)
    hdu_list = fits.open('{}/{}'.format(fpa_dir, fpa_file))

    # todo clarify everything!
    for image_idx, idx in enumerate(images_idxs):
        coords = get_position(image_idx)
        quadrants_l = []

        for quadrant in range(0, 4, 1):
            """
            print('order {} - quadrant {} - dither {}'.format(order, quadrant,
                                                              dither))
            name = 'CCD_{}_q{}_d{}.fits'.format(fits_dir, coords, quadrant,
                                                fits_file[-6:-5])
            """
            quadrants_l.append(hdu_list[idx + quadrant])  # Base image (0) +
                                                          # quadrant

        quadrant_name = 'CCD_{}_d{}'.format(coords, dither)
        quadrants_d[quadrant_name] = quadrants_l

    for key_ in quadrants_d.keys():
        create_ccd(logger, fits_dir, quadrants_d[key_], key_)


def create_ccd(logger, fits_dir, quadrants, key_):
    """

    :param logger:
    :param quadrants:
    :param key_:
    :return:
    """
    logger.debug('Creates file {}'.format(key_))

    prex = quadrants[0].header['PRESCANX']
    ovrx = quadrants[0].header['OVRSCANX']
    try:
        ovry = quadrants[0].header['OVRSCANY']
    except:
        ovry = 0

    # Code for injected lines, when they will be implemented in ELViS
    img = np.zeros([4132, 4096], np.int16)
    for i in range(0, len(quadrants), 1):
        if i == 0:
            img[:2066, :2048] = quadrants[i].data[:, prex:-ovrx]
            #
            # Correct reference pixel coordinate
            hdr = quadrants[i].header
            hdr.set('CRPIX1', hdr['CRPIX1'] - prex)
        if i == 1:
            img[:2066, 2048:] = quadrants[i].data[:, ovrx:-prex]
        if i == 2:
            img[2066:, :2048] = quadrants[i].data[:, prex:-ovrx]
        if i == 3:
            img[2066:, 2048:] = quadrants[i].data[:, ovrx:-prex]
            #
            # Save to FITS file
            outputfits = fits.HDUList()
            outputfits.append(fits.PrimaryHDU())
            outputfits.append(fits.ImageHDU(data=img, header=hdr))
            outputfits.writeto('{}/{}.fits'.format(fits_dir, key_),
                               overwrite=True)

    return True
