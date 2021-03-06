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
import misc
import misc_fits
import multiprocessing
import sextractor_aux_elvis


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
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


class MockedProcess:

    def __init__(self, target, args):
        """

        :param target:
        :param args:
        """
        pass

    def start(self):
        """

        :return:
        """
        return 'start'

    def is_alive(self):
        """

        :return:
        """
        return 'is_alive'


class TestSplitMethodFromCheckElvis(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        sys.argv = ['test_splitMethodCheckElvis.py', '']

    def test_processes_created(self):
        """

        :return:
        """
        settings_d = {'fits_dir': 'fits_dir_mock', 'fpas_dir': 'fpas_dir_mock'}
        misc.extract_settings_elvis = MagicMock(return_value=settings_d)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        fpa_l = ['fpa_1', 'fpa_2', 'fpa_3', 'fpa_4']
        misc_fits.get_fpa_elvis = MagicMock(return_value=fpa_l)
        multiprocessing.Process.side_effect = MagicMock(return_value=MockedProcess)

        sys.argv[1] = '-split'

        self.assertTrue(check_elvis.Check)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
