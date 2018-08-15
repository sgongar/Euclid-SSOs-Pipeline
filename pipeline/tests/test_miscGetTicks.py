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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from errors import BadSettings, WrongTicksList

import misc


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestExtractSettingsElvis(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_get_ticks(self):
        """

        :return:
        """
        ticks = misc.get_ticks(1, 4, 1)

        assert_1 = self.assertIs(type(ticks), list)
        assert_2 = self.assertIs(ticks[0], 1)
        assert_3 = self.assertIs(ticks[1], 2)
        assert_4 = self.assertIs(ticks[2], 3)
        assert_5 = self.assertIs(ticks[3], 4)

        return assert_1 and assert_2 and assert_3 and assert_4 and assert_5

    def test_get_ticks_wrong(self):
        """

        :return:
        """
        return self.assertRaises(WrongTicksList, misc.get_ticks, 0.3, 4, 5)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
