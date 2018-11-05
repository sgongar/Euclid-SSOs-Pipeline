#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.

Todo:
    * Improve documentation

*GNU Terry Pratchett*
"""
from subprocess import Popen

import misc

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class ScampELViS:

    def __init__(self, logger, scmp_d):
        """

        :param logger:
        :param scmp_d:
        """
        self.prfs_d = misc.extract_settings_elvis()

        self.logger = logger
        self.scmp_d = scmp_d

        self.scamp_process()

    def scamp_process(self):
        """

        :return:
        """
        self.logger.info('Scamp process')

        scmp_1 = 'scamp -c {}'.format(self.prfs_d['conf_scamp'])
        scmp_2 = ' {}/*t.cat'.format(self.prfs_d['fits_dir'])
        scmp_3 = ' -ASTREFCAT_NAME {}/stars_catalogue.cat'.format(self.prfs_d['references'])
        scmp_4 = ' -PIXSCALE_MAXERR {}'.format(self.scmp_d['pixscale_maxerr'])
        scmp_5 = ' -POSANGLE_MAXERR {}'.format(self.scmp_d['posangle_maxerr'])
        scmp_6 = ' -POSITION_MAXERR {}'.format(self.scmp_d['position_maxerr'])
        scmp_7 = ' -CROSSID_RADIUS {}'.format(self.scmp_d['crossid_radius'])
        # Output catalogs location
        merged_cat_n = '{}/merged.cat'.format(self.prfs_d['output_cats'])
        scmp_8 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat_n)
        full_cat_n = '{}/full.cat'.format(self.prfs_d['output_cats'])
        scmp_9 = ' -FULLOUTCAT_NAME {}'.format(full_cat_n)
        scmp_p = scmp_1 + scmp_2 + scmp_3 + scmp_4 + scmp_5
        scmp_p = scmp_p + scmp_6 + scmp_7 + scmp_8 + scmp_9
        scmp_p = scmp_p

        process_scamp = Popen(scmp_p, shell=True)
        process_scamp.wait()

        self.logger.info('Scamp process finished.')

        return True
