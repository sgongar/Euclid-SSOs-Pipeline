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
from errors import AllSameException
import misc

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestAllSame(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_all_sources_are_equal(self):
        """

        :return: (True, 3)
        """
        sources_list = [1, 1, 1]
        output = misc.all_same(sources_list)

        statement_1 = self.assertIs(type(output), tuple)
        statement_2 = self.assertTrue(output[0])
        statement_3 = self.assertIs(output[1], 3)

        return statement_1 and statement_2 and statement_3

    def test_all_sources_are_false(self):
        """

        :return: (False, 0)
        """
        sources_list = ['False', 'False', 'False']
        output = misc.all_same(sources_list)

        statement_1 = self.assertIs(type(output), tuple)
        statement_2 = self.assertFalse(output[0])
        statement_3 = self.assertIs(output[1], 0)

        return statement_1 and statement_2 and statement_3

    def test_there_are_four_values_three_equal_sources_one_false(self):
        """

        :return: (True, 3)
        """
        sources_list = [1, 1, 1, 'False']
        output = misc.all_same(sources_list)

        statement_1 = self.assertIs(type(output), tuple)
        statement_2 = self.assertTrue(output[0])
        statement_3 = self.assertIs(output[1], 3)

        return statement_1 and statement_2 and statement_3

    def test_there_are_four_values_three_false_one_source(self):
        """

        :return: (False, 0)
        """
        sources_list = [1, 'False', 'False', 'False']
        output = misc.all_same(sources_list)

        statement_1 = self.assertIs(type(output), tuple)
        statement_2 = self.assertFalse(output[0])
        statement_3 = self.assertIs(output[1], 0)

        return statement_1 and statement_2 and statement_3

    def test_there_are_more_than_two_different_outputs(self):
        """

        :return: (False, 0)
        """
        sources_list = [1, 2, 'False', 'False']
        output = misc.all_same(sources_list)

        statement_1 = self.assertIs(type(output), tuple)
        statement_2 = self.assertFalse(output[0])
        statement_3 = self.assertIs(output[1], 0)

        return statement_1 and statement_2 and statement_3

    def test_there_are_three_values_sources_and_false_values(self):
        """

        :return: (False, 0)
        """
        sources_list = [1, 2, 'False']
        output = misc.all_same(sources_list)

        statement_1 = self.assertIs(type(output), tuple)
        statement_2 = self.assertFalse(output[0])
        statement_3 = self.assertIs(output[1], 0)

        return statement_1 and statement_2 and statement_3

    def test_there_is_an_empty_list(self):
        """

        :return: AllSameException
        """
        sources_list = []

        return self.assertRaises(AllSameException,
                                 misc.all_same, sources_list)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
