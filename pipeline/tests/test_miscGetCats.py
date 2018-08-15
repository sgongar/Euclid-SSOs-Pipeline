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


class TestGetCats(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_there_are_catalogs(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        cat_list = ['m_20-21_cat1.cat', 'm_20-21_cat2.cat',
                    'm_20-21_cat3.cat', 'm_20-21_cat4.cat']
        os.listdir = MagicMock(return_value=cat_list)

        output = misc.get_cats()

        statement_1 = self.assertIs(type(output), dict)
        statement_2 = self.assertIs(len(output.keys()), 1)
        statement_3 = self.assertIs(len(output[output.keys()[0]]), 4)

        return statement_1 and statement_2 and statement_3

    def test_there_are_no_catalogs(self):
        """

        :return:
        """
        prfs_d = {'fits_dir': 'fits_dir'}
        misc.extract_settings = MagicMock(return_value=prfs_d)
        cat_list = ['w_20-21_cat1.cat', 'w_20-21_cat2.cat',
                    'w_20-21_cat3.cat', 'w_20-21_cat4.cat']
        os.listdir = MagicMock(return_value=cat_list)

        output = misc.get_cats()

        statement_1 = self.assertIs(type(output), dict)
        statement_2 = self.assertIs(len(output.keys()), 1)
        statement_3 = self.assertIs(len(output[output.keys()[0]]), 0)

        return statement_1 and statement_2 and statement_3

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
