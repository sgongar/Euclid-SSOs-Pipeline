#!/usr/bin/python
# -*- coding: utf-8 -*-

""" #TODO

Versions:
* 0.1 - 

In order to improve #TODO
* c = catalog
* n = name
* d = dictionary
* loc = location
* cmp = comparation

Todo:
    * Improve log messages
    * Improve docstring
    * Create versions history
    * License??
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pandas import read_csv

from misc import extract_settings, get_fits

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class DrawHistograms:

    def __init__(self):
        """

        """
        prfs_d = extract_settings()
        distance_d = self.read_stats(prfs_d)

        if not self.draw_histograms(distance_d):
            raise Exception

    def read_stats(self, prfs_d):
        """ Reads stats.

        @param prfs_d:

        @return distance_d:
        """
        # Creates dict.
        distance_d = {}
        # Gets fits files from directory.
        fits_files = get_fits(unique=True)

        for fits_n in fits_files:
            csv_n = '{}/sources_{}.csv'.format(prfs_d['tmp_out'],
                                               fits_n[-13:-5])
            cat = read_csv(csv_n, index_col=0)
            distance_d[fits_n[-13:-5]] = cat['distance'].tolist()

        return distance_d

    def draw_histograms(self, distance_d):
        """

        @param distance_d:

        """
        fig, ax_l = plt.subplots(ncols=3, nrows=3,
                                 figsize=(16.53, 11.69), dpi=100)

        ((ax0, ax1, ax2), (ax3, ax4, ax5), (ax6, ax7, ax8)) = ax_l

        keys = ['x0_y0_d1', 'x0_y1_d1', 'x0_y2_d1',
                'x1_y0_d1', 'x1_y1_d1', 'x1_y2_d1',
                'x2_y0_d1', 'x2_y1_d1', 'x2_y2_d1']

        idx_l = 0
        for idx_g in range(0, 7, 3):
            for idx_ax in range(0, 3, 1):
                key_ = keys[idx_g + idx_ax]
                ticks = [0.001, 0.002, 0.003, 0.004, 0.005]

                ax_l[idx_l][idx_ax].hist(distance_d[key_])
                ax_l[idx_l][idx_ax].set_title(key_)
                ax_l[idx_l][idx_ax].grid(True)
                # ax_l[idx_l][idx_ax].set_xticks(ticks)
            idx_l += 1

        fig.tight_layout()

        with PdfPages('test.pdf') as pdf:
            pdf.savefig()


if __name__ == '__main__':

    DrawHistograms()
