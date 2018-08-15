#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Plots pm and magnitude against A/B/elongation
Galaxies data

Versions:
- 0.1 Initial release

Todo:
    * Improve usability

"""
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import pyplot
import numpy as np
from numpy import arange
from pandas import read_csv


class PlotMagVSPMStars:

    def __init__(self):
        stars_df = read_csv('stars_df.csv', index_col=0)
        self.stars_d = stars_df.to_dict()

        self.data_d = {}
        self.manage_stars_dict()
        self.plot_figure()

    def manage_stars_dict(self):
        """

        :return:
        """
        # Stars
        stars_mag = []
        for key_ in self.stars_d['stars_mag'].keys():
            stars_mag.append(self.stars_d['stars_mag'][key_])
        stars_pm = []
        for key_ in self.stars_d['stars_pm'].keys():
            stars_pm.append(self.stars_d['stars_pm'][key_])

        self.data_d['x_stars'] = np.array(stars_mag)
        self.data_d['y_stars'] = np.array(stars_pm)

    def plot_figure(self):
        """

        :return:
        """
        fig = pyplot.figure(figsize=(16.53, 11.69), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        # ax.semilogy(self.data_d['x_galaxies'], self.data_d['y_galaxies'],
        #             'bs', ms=1)
        ax.plot(self.data_d['x_stars'], self.data_d['y_stars'], 'bs', ms=1)
        ax.set_ylim([0.001, 1])
        # ax.grid(b=True, which='major', ls='-', lw=2)
        # ax.grid(b=True, which='minor', ls='--', lw=1)

        pyplot.grid(True)
        pyplot.show()


class PlotMagVSPMGalaxies:

    def __init__(self):
        self.galaxies_df = read_csv('galaxies_df.csv', index_col=0)
        self.galaxies_d = self.galaxies_df.to_dict()

        self.data_d = {}
        self.manage_galaxies_dict()

        mag_gaps = [[16, 18], [18, 20], [20, 22],
                    [22, 24], [24, 26], [26, 28]]

        pdf_name = 'histograms.pdf'
        with PdfPages(pdf_name) as pdf:
            fig = self.plot_figure()
            pdf.savefig()
            pyplot.clf()
            pyplot.close(fig)
            for gap_ in mag_gaps:
                fig = self.plot_histograms(gap_)
                pdf.savefig()
                pyplot.clf()
                pyplot.close(fig)

    def manage_galaxies_dict(self):
        """

        :return:
        """
        # Galaxies
        galaxies_mag = []
        for key_ in self.galaxies_d['galaxies_mag'].keys():
            galaxies_mag.append(self.galaxies_d['galaxies_mag'][key_])
        galaxies_pm = []
        for key_ in self.galaxies_d['galaxies_pm'].keys():
            galaxies_pm.append(self.galaxies_d['galaxies_pm'][key_])

        self.data_d['x_galaxies'] = np.array(galaxies_mag)
        self.data_d['y_galaxies'] = np.array(galaxies_pm)

    def plot_figure(self):
        """

        :return:
        """
        fig = pyplot.figure(figsize=(16.53, 11.69), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        # ax.semilogy(self.data_d['x_galaxies'], self.data_d['y_galaxies'],
        #             'bs', ms=1)
        ax.plot(self.data_d['x_galaxies'], self.data_d['y_galaxies'], 'bs',
                ms=1)
        ax.set_ylim([0.001, 1])
        # ax.grid(b=True, which='major', ls='-', lw=2)
        # ax.grid(b=True, which='minor', ls='--', lw=1)

        pyplot.grid(True)

        return fig

    def plot_histograms(self, gap):
        """

        :param gap:
        :return:
        """
        fig = pyplot.figure(figsize=(16.53, 11.69), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        galaxies_df = self.galaxies_df[self.galaxies_df['galaxies_mag'] < gap[1]]
        galaxies_df = galaxies_df[galaxies_df['galaxies_mag'] > gap[0]]

        pyplot.hist(galaxies_df['galaxies_pm_err'], arange(0, 0.5, 0.005),
                    normed=1, facecolor='green', alpha=0.75)

        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 15.0)

        ax.set_xticks(arange(0, 0.6, 0.1), minor=False)
        ax.set_xticks(arange(0, 0.5, 0.005), minor=True)
        ax.set_yticks(arange(0, 15, 2), minor=False)
        ax.set_yticks(arange(0, 15, 1), minor=True)

        ax.grid(b=True, which='major', ls='-', lw=2)
        ax.grid(b=True, which='minor', ls='--', lw=1)

        return fig




if __name__ == "__main__":
    test = PlotMagVSPMGalaxies()
