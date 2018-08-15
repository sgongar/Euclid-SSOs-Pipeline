#!/usr/bin/python
# -*- coding: utf-8 -*-


class UnsupportedFile(Warning):

    """
    One of the following situations may raise this error:
        1. A file not valid for sextracting process has been found
    """


class TooManyDetections(Warning):

    """
    One of the following situations may raise this warning:
        1. A multiple detections for the same object has been done
    """
