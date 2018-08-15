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


class TestSextractorELViSThreadsCreated(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_sextractor_threads_list_created(self):
        """

        :return:
        """
        logger = MagicMock(side_effect=MockedLogger)
        analysis_dict = {}
        analysis_d = MagicMock(return_value=analysis_dict)
        os.listdir = MagicMock(return_value=['fits_1', 'fits_2',
                                             'fits_3', 'fits_4'])
        misc.extract_settings_elvis = MagicMock(return_value=True)
        test_SextractorELViS = sextractor_aux_elvis.SextractorELViS(logger,
                                                                    analysis_d)
        print(test_SextractorELViS.sextractor_process())

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
