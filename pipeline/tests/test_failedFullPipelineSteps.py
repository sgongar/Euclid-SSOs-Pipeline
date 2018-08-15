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

import check_elvis
from errors import FullPipelineFailed, CleanFailed, SplitFailed
from errors import SextractorFailed, ScampFailed, FiltFailed, RestartFailed
from errors import ChangeTimeFailed
import times_elvis
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


class TestFullPipelineUnsuccessfulSteps(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        sys.argv = ['test_checkOptions.py', '']

    def test_restart_fails(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        check_elvis.Check.restart = MagicMock(return_value=False)

        sys.argv[1] = '-full'

        return self.assertRaises(RestartFailed, check_elvis.Check)

    def test_split_fails(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        check_elvis.Check.restart = MagicMock(return_value=True)
        check_elvis.Check.split = MagicMock(return_value=False)

        sys.argv[1] = '-full'

        return self.assertRaises(SplitFailed, check_elvis.Check)

    def test_change_times_fails(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        check_elvis.Check.restart = MagicMock(return_value=True)
        check_elvis.Check.split = MagicMock(return_value=True)
        times_elvis.change_times = MagicMock(return_value=False)

        sys.argv[1] = '-full'

        return self.assertRaises(ChangeTimeFailed, check_elvis.Check)

    def test_sextractor_fails(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        check_elvis.Check.restart = MagicMock(return_value=True)
        check_elvis.Check.split = MagicMock(return_value=True)
        times_elvis.change_times = MagicMock(return_value=True)
        check_elvis.Check.sextractor = MagicMock(return_value=False)

        sys.argv[1] = '-full'

        return self.assertRaises(SextractorFailed, check_elvis.Check)

    def test_scamp_fails(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        check_elvis.Check.restart = MagicMock(return_value=True)
        check_elvis.Check.split = MagicMock(return_value=True)
        times_elvis.change_times = MagicMock(return_value=True)
        check_elvis.Check.sextractor = MagicMock(return_value=True)
        check_elvis.Check.scamp = MagicMock(return_value=False)

        sys.argv[1] = '-full'

        return self.assertRaises(ScampFailed, check_elvis.Check)

    def test_filter_fails(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)
        check_elvis.Check.restart = MagicMock(return_value=True)
        check_elvis.Check.split = MagicMock(return_value=True)
        times_elvis.change_times = MagicMock(return_value=True)
        check_elvis.Check.sextractor = MagicMock(return_value=True)
        check_elvis.Check.scamp = MagicMock(return_value=True)
        check_elvis.Check.filt = MagicMock(return_value=False)

        sys.argv[1] = '-full'

        return self.assertRaises(FiltFailed, check_elvis.Check)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
