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
from mock import MagicMock, patch

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from check_elvis import Check
from errors import ScampFailed


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

    def debug(self, text):
        """

        :param text:
        :return:
        """
        pass


class TestScampMethodFromCheckElvis(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('scamp_aux_elvis.ScampELViS')
    def test_scamp_works(self, ScampELViS, setting_logger,
                         extract_settings_elvis):
        """

        :param ScampELViS:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        ScampELViS.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-scamp'

        self.assertTrue(Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('scamp_aux_elvis.ScampELViS')
    def test_scamp_fails(self, ScampELViS, setting_logger,
                         extract_settings_elvis):
        """

        :param ScampELViS:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        ScampELViS.return_value = False
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-scamp'

        self.assertRaises(ScampFailed, Check)

    def tearDown(self):
        """

        :return:
        """
        pass

if __name__ == '__main__':
    main()
