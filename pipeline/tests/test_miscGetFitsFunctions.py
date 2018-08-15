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


class TestGetFits(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_all_files_there_are_files(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        right_files_list = ['m20-21_1.fits', 'm20-21_2.fits',
                            'm20-21_3.fits', 'm20-21_4.fits']
        os.listdir = MagicMock(return_value=right_files_list)

        return_get_fits = misc_fits.get_fits(False, '20-21')
        statement_1 = self.assertEqual(type(return_get_fits), list)
        statement_2 = self.assertEqual(len(return_get_fits), 4)

        return statement_1 and statement_2

    def test_all_files_there_are_no_files(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        wrong_files_list = ['w20-21_1.fits', 'w20-21_2.fits',
                            'w20-21_3.fits', 'w20-21_4.fits']
        os.listdir = MagicMock(return_value=wrong_files_list)

        return_get_fits = misc_fits.get_fits(False, '20-21')
        statement_1 = self.assertEqual(type(return_get_fits), list)
        statement_2 = self.assertEqual(len(return_get_fits), 0)

        return statement_1 and statement_2

    def test_unique_files_there_are_files(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        wrong_files_list = ['m12345678901234567890.fits',
                            'm23456789012345678901.fits',
                            'm34567890123456789012.fits',
                            'm45678901234567890123.fits']
        os.listdir = MagicMock(return_value=wrong_files_list)

        return_get_fits = misc_fits.get_fits(True, '20-21')
        statement_1 = self.assertEqual(type(return_get_fits), list)
        statement_2 = self.assertEqual(len(return_get_fits), 4)

        return statement_1 and statement_2

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
