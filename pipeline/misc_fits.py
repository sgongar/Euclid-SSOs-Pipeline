#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements


Todo:
    * Improve log messages
    * Improve usability
"""
import os

import misc

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_fits(unique, mag):
    """

    :param unique:
    :param mag:
    :return:
    """
    prfs_d = misc.extract_settings()
    fits_list = []

    fits_dir = '{}/{}/CCDs/'.format(prfs_d['fits_dir'], mag)

    files = os.listdir('{}'.format(fits_dir))
    for file_ in files:
        if file_[:1] == 'm' and file_[-5:] == '.fits':
            fits_list.append(file_)

    if unique:
        fits_unique = []
        for file_ in fits_list:
            fits_unique.append(file_[:21])
        fits_unique = list(set(fits_unique))

        for file_ in fits_unique:
            fits_unique[fits_unique.index(file_)] = file_ + '1.fits'

        return fits_unique
    else:
        return fits_list


def get_fpa_elvis():
    """

    :return:
    """
    prfs_d = misc.extract_settings_elvis()
    fits_list = []

    files = os.listdir('{}'.format(prfs_d['fpas_dir']))
    for file_ in files:
        if file_[-5:] == '.fits':
            fits_list.append(file_)

    return fits_list


def get_fits_d(mag_, dither):
    """ Gets a list of fits files by dither.

    :param mag_:
    :param dither:
    :return:
    """
    prfs_d = misc.extract_settings()
    fits_list = []

    files = os.listdir('{}/{}/CCDs/'.format(prfs_d['fits_dir'], mag_))
    for file_ in files:
        if file_[:1] == 'm' and file_[-5:] == '.fits':
            fits_list.append(file_)

    list_out = []
    for file_ in fits_list:
        if file_[-6:-5] == str(dither):
            list_out.append(file_)

    return list_out
