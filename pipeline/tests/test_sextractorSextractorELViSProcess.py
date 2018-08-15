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

import sextractor_aux_elvis
import misc


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


class TestSextractorELViSProcessCalled(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_sextractor_process_called(self):
        """

        :return:
        """
        logger = MagicMock(side_effect=MockedLogger)
        analysis_dict = {}
        analysis_d = MagicMock(return_value=analysis_dict)
        misc.extract_settings_elvis = MagicMock(return_value=True)
        sextractor_aux_elvis.SextractorELViS.sextractor_process = MagicMock(return_value=True)

        return self.assertTrue(sextractor_aux_elvis.SextractorELViS(logger, analysis_d))

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
