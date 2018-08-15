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

import check_elvis
import cosmic_elvis
from errors import CleanFailed
import misc


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class MockedLogger:
    def __init__(self, text):
        """

        :param text:
        """
        pass

    def info(self, text):
        """

        :param text:
        :return:
        """
        pass


class TestCleanMethodFromCheckElvis(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        sys.argv = ['test_cleanMethodCheckElvis.py', '']

    def test_clean_works(self):
        """

        :return:
        """
        settings_d = {'fits_dir': 'fits_dir_mock', 'fpas_dir': 'fpas_dir_mock'}
        misc.extract_settings_elvis = MagicMock(return_value=settings_d)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        cosmic_elvis.CosmicELViS = MagicMock(return_value=True)

        sys.argv[1] = '-clean'

        self.assertTrue(check_elvis.Check)

    def test_clean_fails(self):
        """

        :return:
        """
        settings_d = {'fits_dir': 'fits_dir_mock', 'fpas_dir': 'fpas_dir_mock'}
        misc.extract_settings_elvis = MagicMock(return_value=settings_d)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        cosmic_elvis.CosmicELViS = MagicMock(return_value=False)

        sys.argv[1] = '-clean'

        self.assertRaises(CleanFailed, check_elvis.Check)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
