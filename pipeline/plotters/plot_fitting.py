#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1 Initial release

Todo:
    *

"""
from datetime import datetime, timedelta
from math import modf

from astropy import units
from astropy.coordinates import Angle

from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from pyds9 import DS9

from misc import extract_settings
from plot_aux import create_ticks


class PlotFitting:

    def __init__(self, star, data_d, output_path, fits_d):
        """

        :param plot_d:
        """
        self.star = star
        self.prfs_d = extract_settings()
        self.data_d = data_d
        self.output_path = output_path
        self.fits_d = fits_d

        # Page size
        self.plot_size = [16.53, 11.69]
        self.plot_dpi = 100

        self.plot()

    def odr_regression(self):
        """

        :return: i_alpha_seconds, i_delta_seconds
        """
        x_odr = []
        tmp_hour = []
        tmp_minute = []
        for alpha_ in self.fits_d['x_odr']:
            a = Angle(alpha_, units.degree)
            dms = a.dms
            degree = int(dms[0])
            tmp_hour.append(degree)
            minute = int(dms[1])
            tmp_minute.append(minute)
            second = float("{0:.6f}".format(dms[2]))
            x_odr.append('{}'.format(second))

        # delta coordinates
        y_odr = []
        tmp_degree = []
        tmp_minute = []
        for delta_ in self.fits_d['y_odr']:
            d = Angle(delta_, units.degree)
            dms = d.dms
            degree = int(dms[0])
            tmp_degree.append(degree)
            minute = int(dms[1])
            tmp_minute.append(minute)
            second = float("{0:.6f}".format(dms[2]))
            y_odr.append('{}'.format(second))

        x_odr = [float(i) for i in x_odr]  # needed?
        x_odr = [float("{0:.6f}".format(i)) for i in x_odr]
        y_odr = [float(i) for i in y_odr]  # needed?
        y_odr = [float("{0:.6f}".format(i)) for i in y_odr]

        return x_odr, y_odr

    def i_coordinates(self):
        """

        :return: i_alpha_seconds, i_delta_seconds
        """
        i_alpha_tmpstmp = []
        i_alpha_seconds = []
        i_alpha_arcseconds = []
        tmp_hour = []
        tmp_minute = []
        for alpha_ in self.data_d['i_alpha']:
            a = Angle(alpha_, units.degree)
            dms = a.dms
            degree = int(dms[0])
            tmp_hour.append(degree)
            minute = int(dms[1])
            tmp_minute.append(minute)
            second = float("{0:.6f}".format(dms[2]))
            arcsecond = float("{0:.6f}".format(a.arcsecond))
            i_alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute,
                                                     second))
            i_alpha_seconds.append('{}'.format(second))
            i_alpha_arcseconds.append(arcsecond)

        # delta coordinates
        i_delta_tmpstmp = []
        i_delta_seconds = []
        tmp_degree = []
        tmp_minute = []
        for delta_ in self.data_d['i_delta']:
            d = Angle(delta_, units.degree)
            dms = d.dms
            degree = int(dms[0])
            tmp_degree.append(degree)
            minute = int(dms[1])
            tmp_minute.append(minute)
            second = float("{0:.6f}".format(dms[2]))
            i_delta_tmpstmp.append(
                '{}:{}.{}'.format(degree, minute, second))
            i_delta_seconds.append('{}'.format(second))

        i_alpha_seconds = [float(i) for i in i_alpha_seconds]  # needed?
        i_alpha_seconds = [float("{0:.6f}".format(i)) for i in i_alpha_seconds]
        i_delta_seconds = [float(i) for i in i_delta_seconds]  # needed?
        i_delta_seconds = [float("{0:.6f}".format(i)) for i in i_delta_seconds]

        return i_alpha_seconds, i_delta_seconds

    def o_coordinates(self):
        """

        :return: o_alpha_seconds, o_delta_seconds
        """
        o_alpha_tmpstmp = []
        o_alpha_seconds = []
        tmp_hour = []
        tmp_minute = []
        for alpha_ in self.data_d['o_alpha']:
            if alpha_ is not False:
                a = Angle(alpha_, units.degree)
                dms = a.dms
                degree = int(dms[0])
                tmp_hour.append(degree)
                minute = int(dms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(dms[2]))
                o_alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute,
                                                         second))
                o_alpha_seconds.append('{}'.format(second))

        o_delta_tmpstmp = []
        o_delta_seconds = []
        tmp_degree = []
        tmp_minute = []
        for delta_ in self.data_d['o_delta']:
            if delta_ is not False:
                d = Angle(delta_, units.degree)
                dms = d.dms
                degree = int(dms[0])
                tmp_degree.append(degree)
                minute = int(dms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(dms[2]))
                o_delta_tmpstmp.append('{}:{}.{}'.format(degree, minute,
                                                         second))
                o_delta_seconds.append('{}'.format(second))

        o_alpha_seconds = [float(i) for i in o_alpha_seconds]  # needed?
        o_alpha_seconds = [float("{0:.6f}".format(i)) for i in o_alpha_seconds]
        o_delta_seconds = [float(i) for i in o_delta_seconds]  # needed?
        o_delta_seconds = [float("{0:.6f}".format(i)) for i in o_delta_seconds]

        return o_alpha_seconds, o_delta_seconds

    def general_view(self, pdf):
        """ shows both input and output extractions

        :param pdf:
        :return:
        """
        fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
        ax_1 = plt.subplot2grid((1, 5), (0, 0), colspan=5)

        source_str = 'source {}'.format(self.data_d['source'][0])
        pm_str = 'input_pm {}'.format(self.data_d['i_pm'][0])
        # ok_str = 'ok {}'.format(self.ok)
        ax_1.set_title('{}\n{}'.format(source_str, pm_str))

        i_alpha_seconds, i_delta_seconds = self.i_coordinates()
        o_alpha_seconds, o_delta_seconds = self.o_coordinates()

        # alpha coordinates

        alpha_seconds = i_alpha_seconds + o_alpha_seconds
        delta_seconds = i_delta_seconds + o_delta_seconds

        alpha_seconds = [float(i) for i in alpha_seconds]  # needed?
        alpha_seconds = [float("{0:.6f}".format(i)) for i in alpha_seconds]
        delta_seconds = [float(i) for i in delta_seconds]  # needed?
        delta_seconds = [float("{0:.6f}".format(i)) for i in delta_seconds]

        # Plots data
        blue_colors = cm.Blues(np.linspace(0, 1, len(i_alpha_seconds) + 1))
        for idx in range(0, len(i_alpha_seconds), 1):
            ax_1.scatter(i_alpha_seconds[idx],
                         i_delta_seconds[idx],
                         c=blue_colors[idx + 1], s=32)

        red_colors = cm.Reds(np.linspace(0, 1, len(o_alpha_seconds) + 1))
        for idx in range(0, len(o_alpha_seconds), 1):
            ax_1.scatter(o_alpha_seconds[idx],
                         o_delta_seconds[idx],
                         c=red_colors[idx + 1], s=32)

        try:
            x_ticks = create_ticks(alpha_seconds)
            y_ticks = create_ticks(delta_seconds)
        except ZeroDivisionError:
            print('alpha: {}'.format(self.data_d['i_alpha']))
            print('delta: {}'.format(self.data_d['i_delta']))
            print('source {}'.format(self.data_d['source'][0]))
            raise Exception

        # x-ticks assignation
        ax_1.set_xticks(x_ticks['major_t'], minor=False)
        ax_1.set_xticklabels(x_ticks['major_t'])
        ax_1.set_xticks(x_ticks['minor_t'], minor=True)

        # y-ticks assignation
        ax_1.set_yticks(y_ticks['major_t'], minor=False)
        ax_1.set_yticklabels(y_ticks['major_t'])
        ax_1.set_yticks(y_ticks['minor_t'], minor=True)

        # Formats grids
        ax_1.grid(b=True, which='major', linestyle='-', linewidth=2)
        ax_1.grid(b=True, which='minor', linestyle='--', linewidth=1)

        # x-axis
        x_label_ra = 'Right ascension \n'
        major_s = x_ticks['major_s']
        x_label_major_step = 'major step size {}"\n'.format(major_s)
        minor_s = x_ticks['minor_s']
        x_label_minor_step = 'minor step size {}"'.format(minor_s)
        ax_1.set_xlabel('{}{}{}'.format(x_label_ra, x_label_major_step,
                                        x_label_minor_step))
        # x-axis
        y_label_ra = 'Declination \n'
        major_s = y_ticks['major_s']
        y_label_major_step = 'major step size {}"\n'.format(major_s)
        minor_s = y_ticks['minor_s']
        y_label_minor_step = 'minor step size {}"'.format(minor_s)
        ax_1.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                        y_label_minor_step))

        # Plots a legend
        plt.legend(loc=0, ncol=2, borderaxespad=0.)
        plt.setp(ax_1.get_xticklabels(), visible=True)
        plt.setp(ax_1.get_yticklabels(), visible=True)
        plt.draw()

        pdf.savefig()
        plt.close(fig)

    def fit_view(self, pdf):
        """ shows extractions and fittings

        :rtype: object
        :param pdf:
        :return:
        """

        # Second page
        fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
        ax_1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)

        source_str = 'source {}'.format(self.data_d['source'][0])
        if self.star:
            pm_str = 'input_pm 0'
        else:
            pm_str = 'input_pm {}'.format(self.data_d['i_pm'][0])
        # ok_str = 'ok {}'.format(self.ok)
        ax_1.set_title('{}\n{}'.format(source_str, pm_str))

        alpha_seconds, delta_seconds = self.o_coordinates()

        # Plots data
        red_colors = cm.Reds(np.linspace(0, 1, len(alpha_seconds) + 1))
        for idx in range(0, len(alpha_seconds), 1):
            ax_1.scatter(alpha_seconds[idx], delta_seconds[idx],
                         c=red_colors[idx + 1], s=36)

        # Plots fitting lines
        x_odr, y_odr = self.odr_regression()
        ax_1.plot(x_odr, y_odr, 'g', label='fit_odr')

        try:
            x_ticks = create_ticks(alpha_seconds)
            y_ticks = create_ticks(delta_seconds)
        except ZeroDivisionError:
            print('alpha: {}'.format(self.data_d['i_alpha']))
            print('delta: {}'.format(self.data_d['i_delta']))
            print('source {}'.format(self.data_d['source'][0]))
            # raise Exception

        # x-ticks assignation
        ax_1.set_xticks(x_ticks['major_t'], minor=False)
        ax_1.set_xticklabels(x_ticks['major_t'])
        ax_1.set_xticks(x_ticks['minor_t'], minor=True)

        # y-ticks assignation
        ax_1.set_yticks(y_ticks['major_t'], minor=False)
        ax_1.set_yticklabels(y_ticks['major_t'])
        ax_1.set_yticks(y_ticks['minor_t'], minor=True)

        # Formats grids
        ax_1.grid(b=True, which='major', linestyle='-', linewidth=2)
        ax_1.grid(b=True, which='minor', linestyle='--', linewidth=1)

        # x-axis
        x_label_ra = 'Right ascension \n'
        major_s = x_ticks['major_s']
        x_label_major_step = 'major step size {}"\n'.format(major_s)
        minor_s = x_ticks['minor_s']
        x_label_minor_step = 'minor step size {}"'.format(minor_s)
        ax_1.set_xlabel('{}{}{}'.format(x_label_ra, x_label_major_step,
                                        x_label_minor_step))
        # x-axis
        y_label_ra = 'Declination \n'
        major_s = y_ticks['major_s']
        y_label_major_step = 'major step size {}"\n'.format(major_s)
        minor_s = y_ticks['minor_s']
        y_label_minor_step = 'minor step size {}"'.format(minor_s)
        ax_1.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                        y_label_minor_step))

        # Plots a legend
        plt.legend(loc=0, ncol=2, borderaxespad=0.)
        plt.setp(ax_1.get_xticklabels(), visible=True)
        plt.setp(ax_1.get_yticklabels(), visible=True)
        plt.draw()

        ax_2 = plt.subplot2grid((1, 5), (0, 4))

        if self.star:
            i_pm = 0.0
            i_pm_alpha = 0.0
            i_pm_delta = 0.0
        else:
            i_pm = float("{0:.6f}".format(self.data_d['i_pm'][1]))
            i_pm_alpha = float("{0:.6f}".format(self.data_d['i_pm_alpha'][1]))
            i_pm_delta = float("{0:.6f}".format(self.data_d['i_pm_delta'][1]))
        o_pm = float("{0:.6f}".format(self.data_d['o_pm'][1]))
        o_pm_alpha = float("{0:.6f}".format(self.data_d['o_pm_alpha'][1]))
        o_pm_delta = float("{0:.6f}".format(self.data_d['o_pm_delta'][1]))
        fit_odr = float("{0:.6f}".format(self.fits_d['fit_odr']))
        # print('fit_odr {} - format {}'.format(fit_odr, type(fit_odr)))

        table_ = [['', 'catalog values'],
                  ['cat_pm', i_pm],
                  ['cat_pm_alpha', i_pm_alpha],
                  ['cat_pm_delta', i_pm_delta],
                  ['', 'extracted values'],
                  ['ext_pm', o_pm],
                  ['ext_pm_alpha', o_pm_alpha],
                  ['ext_pm_delta', o_pm_delta],
                  ['', 'goodness of fitting'],
                  ['fit_odr', fit_odr]]

        # ax_2.axis('tight')
        ax_2.axis('off')
        # ax_2.axis('on')
        data_table = ax_2.table(cellText=table_, colLabels=None,
                                loc='center')
        data_table.set_fontsize(14)

        pdf.savefig()
        plt.close(fig)

    # def fits_images(self, pdf):
    #     """
    #
    #     :param pdf:
    #     :return:
    #     """
    #     # Gets limits
    #     alpha_center = 0
    #     delta_center = 0
    #     if len(self.data_d['i_alpha']) == 2:
    #         alpha_sum = self.data_d['i_alpha'][0] + self.data_d['i_alpha'][1]
    #         alpha_center = alpha_sum / 2
    #         delta_sum = self.data_d['i_delta'][0] + self.data_d['i_delta'][1]
    #         delta_center = delta_sum / 2
    #     elif len(self.data_d['i_alpha']) == 3:
    #         alpha_center = self.data_d['i_alpha'][1]
    #         delta_center = self.data_d['i_delta'][1]
    #     elif len(self.data_d['i_alpha']) == 4:
    #         alpha_sum = self.data_d['i_alpha'][1] + self.data_d['i_alpha'][2]
    #         alpha_center = alpha_sum / 2
    #         delta_sum = self.data_d['i_delta'][1] + self.data_d['i_delta'][2]
    #         delta_center = delta_sum / 2
    #
    #     if alpha_center is 0 or delta_center is 0:
    #         print('ERROR')
    #         print("self.err_d['i_alpha'] {}".format(self.data_d['i_alpha']))
    #         print("self.err_d['i_delta'] {}".format(self.data_d['i_delta']))
    #
    #         raise Exception
    #
    #     for idx in range(0, len(self.data_d['i_alpha']), 1):
    #         # Get regions for desired fits file
    #         regions_file = get_regions(self.fits_files[idx],
    #                                    self.data_d['sex_cf'][0],
    #                                    self.prfs_d, self.mag)
    #
    #         # Creates image
    #         dither = self.fits_files[idx][-6:-5]
    #         i_regs = '{}/dither_{}.reg'.format(self.prfs_d['dithers_out'],
    #                                            dither)
    #         if idx == 0:
    #             p_alpha = 0
    #             p_delta = 0
    #         else:
    #             p_alpha = self.data_d['i_alpha'][idx - 1]
    #             p_delta = self.data_d['i_delta'][idx - 1]
    #
    #         # Zoomed image
    #         zoom = True
    #         img = image(i_regs, self.fits_files[idx], regions_file,
    #                     alpha_center, delta_center, idx,
    #                     float(self.data_d['i_pm'][0]), p_alpha, p_delta,
    #                     zoom)
    #
    #         img = img[1:-50, :]  # Removes ds9's legend
    #
    #         fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
    #         ax_1 = fig.add_subplot(1, 1, 1)
    #         ax_1.set_title('Dither {} - Enlarged view '.format(dither))
    #
    #         plt.imshow(img)
    #
    #         labels = [item.get_text() for item in ax_1.get_xticklabels()]
    #
    #         empty_string_labels = [''] * len(labels)
    #         ax_1.set_xticklabels(empty_string_labels)
    #
    #         labels = [item.get_text() for item in ax_1.get_yticklabels()]
    #
    #         empty_string_labels = [''] * len(labels)
    #         ax_1.set_yticklabels(empty_string_labels)
    #
    #         ax_1.set_xlabel('right ascension')
    #         ax_1.set_ylabel('declination')
    #
    #         pdf.savefig()
    #         plt.close(fig)
    #
    #         # Wide image
    #         zoom = False
    #         img = image(i_regs, self.fits_files[idx], regions_file,
    #                     alpha_center, delta_center, idx,
    #                     float(self.err_d['i_pm'][0]), p_alpha, p_delta,
    #                     zoom)
    #
    #         img = img[1:-50, :]  # Removes ds9's legend
    #
    #         fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
    #         ax_1 = fig.add_subplot(1, 1, 1)
    #         ax_1.set_title('Dither {} - Wide view '.format(dither))
    #
    #         plt.imshow(img)
    #
    #         labels = [item.get_text() for item in ax_1.get_xticklabels()]
    #
    #         empty_string_labels = [''] * len(labels)
    #         ax_1.set_xticklabels(empty_string_labels)
    #
    #         labels = [item.get_text() for item in ax_1.get_yticklabels()]
    #
    #         empty_string_labels = [''] * len(labels)
    #         ax_1.set_yticklabels(empty_string_labels)
    #
    #         ax_1.set_xlabel('right ascension')
    #         ax_1.set_ylabel('declination')
    #
    #         pdf.savefig()
    #         plt.close(fig)

    def plot(self):
        """

        :return:
        """
        pdf_name = '{}/{}.pdf'.format(self.output_path,
                                      self.data_d['source'][0])

        with PdfPages(pdf_name) as pdf:
            if self.star:
                self.fit_view(pdf)
            else:
                self.general_view(pdf)
                self.fit_view(pdf)

    # def gets_ccd_names(self, fits_):
    #     """
    #
    #     :return:
    #     """
    #     ccds = []
    #     for ccd in fits_:
    #         ccds.append(ccd[-13:-5])
    #
    #     return ccds
