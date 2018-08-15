#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script

Versions:

Problems:
    * setup has been ignored

Todo:
    * Improve log messages
    *

"""
import os
import sys

from unittest import TestCase, main
from mock import MagicMock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import misc
import misc_fits

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestGetFitsD(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_there_are_files_of_dither_1(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        right_files_list = ['m_20-21_CCD1_d1.fits', 'm_20-21_CCD2_d1.fits',
                            'm_20-21_CCD3_d1.fits', 'm_20-21_CCD4_d1.fits']
        os.listdir = MagicMock(return_value=right_files_list)

        return_get_fits = misc_fits.get_fits_d('20-21', 1)
        statement_1 = self.assertEqual(type(return_get_fits), list)
        statement_2 = self.assertEqual(len(return_get_fits), 4)

        return statement_1 and statement_2

    def test_there_are_no_files_of_dither_1(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        right_files_list = ['m_20-21_CCD1_d2.fits', 'm_20-21_CCD2_d2.fits',
                            'm_20-21_CCD3_d2.fits', 'm_20-21_CCD4_d2.fits']
        os.listdir = MagicMock(return_value=right_files_list)

        return_get_fits = misc_fits.get_fits_d('20-21', 1)
        statement_1 = self.assertEqual(type(return_get_fits), list)
        statement_2 = self.assertEqual(len(return_get_fits), 0)

        return statement_1 and statement_2

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
