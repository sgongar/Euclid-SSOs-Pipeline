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


class TestGetFPAs(TestCase):
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
        prfs_d = {'fpas_dir': 'fpas_dir'}
        misc.extract_settings_elvis = MagicMock(return_value=prfs_d)
        right_files_list = ['FPA_D1.fits', 'FPA_D2.fits',
                            'FPA_D3.fits', 'FPA_D4.fits']
        os.listdir = MagicMock(return_value=right_files_list)

        return_get_fpa_elvis = misc_fits.get_fpa_elvis()

        statement_1 = self.assertEqual(type(return_get_fpa_elvis), list)
        statement_2 = self.assertEqual(len(return_get_fpa_elvis), 4)

        return statement_1 and statement_2

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
