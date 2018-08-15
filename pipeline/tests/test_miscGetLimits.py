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


class TestGetLimits(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_get_limits(self):
        """

        :return:
        """
        list_1 = [1, 2, 3, 4]
        list_2 = [2, 3, 4, 5]

        limits = misc.get_limits(list_1, list_2)

        return self.assertIs(limits[0], 1) and self.assertIs(limits[1], 5)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
