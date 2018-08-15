#!/usr/bin/python
# -*- coding: utf-8 -*-

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import arange
from pandas import read_csv, DataFrame

from misc import setting_logger, extract_settings, create_configurations
from misc import create_scamp_dict, get_limits


class PMs_vs_Detections:

    def __init__(self):
        """

        """
        logger = setting_logger()
        prfs_d = extract_settings()
        mode = {'type': 'scamp'}
        confs, total_confs = create_configurations(mode)

        self.pms_vs_dect(logger, prfs_d, confs)

    def pms_vs_dect(self, logger, prfs_d, confs):
        """

        @param logger:
        @param prfs_d:
        @param confs:
        """
        with PdfPages('propermotions_vs_detections.pdf') as pdf:
            for mag_ in prfs_d['mags']:
                for conf in confs:
                    (scmp_d, len_confs) = create_scamp_dict(logger, prfs_d,
                                                            confs.index(conf))
                    for confidence_ in prfs_d['confidences']:
                        if type(self.open_file(logger, mag_, conf,
                                               confidence_)) is DataFrame:
                            table_ = self.open_file(logger, mag_, conf,
                                                    confidence_)
                            print table_
                            self.create_fig(mag_, conf, table_, prfs_d,
                                            confidence_)
                            # Saves current figure
                            pdf.savefig()
                            plt.clf()  # Clear current figure
                        elif self.open_file(logger, mag_, conf,
                                            confidence_) is str:
                            file_ = self.open_file(logger, mag_, conf,
                                                   confidence_)
                            logger.error('file {} not present'.format(file_))

    def open_file(self, logger, mag_, conf, confidence_):
        """

        @param logger:
        @param mag_:
        @param conf:
        @param pm_:

        @return data_table: a DataFrame object ready to be use
        """
        file_ = "stats_{}_{}_{}_{}_['{}'].csv".format(conf[0], conf[1],
                                                      conf[2], conf[3], mag_)

        try:
            table_ = read_csv(file_)
            table_ = table_[table_['confidence'].isin([confidence_])]

            return table_
        except IOError:
            return file_

    def create_fig(self, mag, conf, data_table, prfs_d, confidence_):
        """

        @param mag:
        @param conf:
        @param data_table:
        @param prfs_d:
        @param pm_:

        @return fig:
        """
        fig = plt.figure(figsize=(11.7, 16.5 / 2), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for pm_ in prfs_d['pms']:
            df_1 = data_table[data_table['pm'].isin([pm_])]
            ax.scatter(pm_, df_1['right_scamp'].iloc[0], c='g')
            ax.scatter(pm_, df_1['right_filter'].iloc[0], c='r')

        t = [confidence_, conf[0], conf[2], conf[3]]
        suptitle_1 = 'confidence {} crossid {} '.format(t[0], t[1])
        suptitle_2 = 'posangle {} position {}'.format(t[2], t[3])
        suptitle_t = suptitle_1 + suptitle_2
        fig.suptitle(suptitle_t)

        ax.set_xlabel('proper motions')
        ax.set_ylabel('right detections')

        # x-scale
        ax.set_xlim(0.1, 1)
        # ax.set_xscale('symlog', subsx=[1, 2, 3, 4, 5, 6, 7, 8, 9])
        x_major_ticks = arange(0, 1, 0.5)
        x_minor_ticks = arange(0, 1, 0.1)
        ax.set_xticks(x_major_ticks, minor=False)
        ax.set_xticks(x_minor_ticks, minor=True)

        # y-scale
        ax.set_ylim(0, 20)
        y_major_ticks = arange(0, 20, 5)
        y_minor_ticks = arange(0, 20, 1)
        ax.set_yticks(y_major_ticks, minor=True)
        ax.set_yticks(y_minor_ticks, minor=True)

        ax.grid(b=True, which='major', color='b', linewidth=0.1)
        ax.grid(b=True, which='minor', color='b', linewidth=0.01)

        return fig


class Detections_vs_Tolerance:

    def __init__(self):
        """

        """
        logger = setting_logger()
        prfs_d = extract_settings()
        mode = {'type': 'scamp'}
        confs, total_confs = create_configurations(mode)

        # Launchs instance
        # self.dect_vs_tol(logger, prfs_d, confs)

    def dect_vs_tol(self, logger, prfs_d, confs):
        """

        @param logger:
        @param prfs_d:
        @param confs:
        """
        with PdfPages('detections_vs_tolerance.pdf') as pdf:
            for mag_ in prfs_d['mags']:
                for conf in confs:
                    (scmp_d, len_confs) = create_scamp_dict(logger, prfs_d,
                                                            confs.index(conf))
                    for pm_ in prfs_d['pms']:
                        if type(self.open_file(logger, mag_,
                                               conf, pm_)) is DataFrame:
                            data_table = self.open_file(logger, mag_,
                                                        conf, pm_)
                            self.create_fig(mag_, conf, data_table,
                                            prfs_d, pm_)
                            # Saves current figure
                            pdf.savefig()
                            plt.clf()  # Clear current figure
                        elif self.open_file(logger, mag_, conf, pm_) is str:
                            file_ = self.open_file(logger, mag_, conf, pm_)
                            logger.error('file {} not present'.format(file_))

    def open_file(self, logger, mag_, conf, pm_):
        """

        @param logger:
        @param mag_:
        @param conf:
        @param pm_:

        @return data_table: a DataFrame object ready to be use
        """
        data_file = 'stats_{}_{}_{}_{}_{}.csv'.format(conf[0], conf[1],
                                                      conf[2], conf[3], mag_)

        try:
            data_table = read_csv(data_file)
            data_table = data_table[data_table['pm'].isin([pm_])]

            return data_table
        except IOError:
            return data_file

    def create_fig(self, mag, conf, data_table, prfs_d, pm_):
        """

        @param mag:
        @param conf:
        @param data_table:
        @param prfs_d:
        @param pm_:

        @return fig:
        """
        fig = plt.figure(figsize=(11.7, 16.5 / 2), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for confidence_ in prfs_d['confidences']:
            df_1 = data_table[data_table['confidence'].isin([confidence_])]
            ax.scatter(confidence_, df_1['right_scamp'].iloc[0], c='g')
            ax.scatter(confidence_, df_1['right_filter'].iloc[0], c='r')

        t = [mag, conf[0], conf[2], conf[3]]
        fig.suptitle('pm {} crossid {} posangle {} position {}'.format(t[0],
                                                                       t[1],
                                                                       t[2],
                                                                       t[3]))
        ax.set_xlabel('tolerance levels')
        ax.set_ylabel('right detections')

        log_scale = True

        if log_scale:
            ax.set_xlim(1, 500)
            ax.set_ylim(0, 20)
            ax.set_xscale('symlog', subsx=[1, 2, 3, 4, 5, 6, 7, 8, 9])

            ax.grid(True)
        else:
            ax.set_xlim(1, 500)
            ax.set_ylim(0, 20)
            x_major_ticks = arange(0, 500, 10)
            x_minor_ticks = arange(0, 500, 2)
            y_major_ticks = arange(0, 20, 0.1)
            y_minor_ticks = arange(0, 20, 0.02)

            ax.set_xticks(x_major_ticks)
            ax.set_xticks(x_minor_ticks, minor=True)
            ax.set_yticks(y_major_ticks)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(True)

        return fig


class Detections_vs_Created:

    def __init__(self):
        """

        """
        logger = setting_logger()
        prfs_d = extract_settings()
        mode = {'type': 'scamp'}
        confs, total_confs = create_configurations(mode)

        # Launchs instance
        self.dect_vs_created(logger, prfs_d, confs)

    def dect_vs_created(self, logger, prfs_d, confs):
        with PdfPages('detections_vs_created.pdf') as pdf:
            for mag_ in prfs_d['mags']:
                for conf in confs:
                    (scmp_d, len_confs) = create_scamp_dict(logger, prfs_d,
                                                            confs.index(conf))
                    for pm_ in prfs_d['pms']:
                        for confidence_ in prfs_d['confidences']:
                            self.create_fig(logger, mag_, conf,
                                            prfs_d, pm_, confidence_)
                            # Saves current figure
                            pdf.savefig()

    def open_file(self, logger, mag_, conf, pm_, crossid_):
        """

        @param logger:
        @param mag_:
        @param conf:
        @param pm_:

        @return data_table: a DataFrame object ready to be use
        """
        data_file = 'stats_{}_{}_{}_{}_{}.csv'.format(crossid_, conf[1],
                                                      conf[2], conf[3], mag_)

        try:
            data_table = read_csv(data_file, index_col=0)
            data_table = data_table[data_table['pm'].isin([pm_])]

            return data_table
        except IOError:
            return data_file

    def create_fig(self, logger, mag_, conf, prfs_d, pm_, confidence_):
        """

        @param mag:
        @param conf:
        @param data_table:
        @param prfs_d:
        @param pm_:

        @return fig:
        """
        # Creates a new fig for all crossid values
        fig = plt.figure(figsize=(11.7, 16.5 / 2), dpi=100)
        ax = fig.add_subplot(1, 1, 1)
        # Loops over crossid values
        for crossid_ in prfs_d['cross_ids']:
            if type(self.open_file(logger, mag_, conf,
                                   pm_, crossid_)) is DataFrame:
                data_table = self.open_file(logger, mag_,
                                            conf, pm_, crossid_)
                df_1 = data_table[data_table['confidence'].isin([confidence_])]

                # Gets the name of filter catalog for desired cross_id
                mrgd_n = '{}_{}_{}_{}_{}'.format(crossid_, conf[1], conf[2],
                                                 conf[3], mag_)
                mrgd_n = '{}/filt_{}__5.csv'.format(prfs_d['results_dir'],
                                                    mrgd_n)

                merged_cat = read_csv(mrgd_n, index_col=0)
                merged_list = merged_cat['SOURCE_NUMBER'].tolist()
                # Gets unique elements of filter catalog
                unique_list_len = len(list(set(merged_list)))

                print df_1

                ax.scatter(unique_list_len,
                           df_1['right_filter'].iloc[0], c='r')


            # if there is no file, 
            elif self.open_file(logger, mag_, conf, pm_, crossid_) is str:
                file_ = self.open_file(logger, mag_, conf, pm_)
                logger.error('file {} not present'.format(file_))

                return None

        return fig


        # df_1['right_camp'].iloc[0], 
        """
        ax.scatter(confidence_, df_1['right_scamp'].iloc[0], c='g')
        ax.scatter(confidence_, df_1['right_filter'].iloc[0], c='r')
        """
        """
        fig.suptitle('crossid: {} posangle: {} position: {}'.format(conf[0],
                                                                    conf[2],
                                                                    conf[3]))
        ax.set_xlabel('pm: {}'.format(pm_))
        ax.set_ylabel('right detections')

        log_scale = True

        if log_scale:
            ax.set_xlim(1, 500)
            ax.set_ylim(0, 20)
            ax.set_xscale('symlog', subsx=[1, 2, 3, 4, 5, 6, 7, 8, 9])

            ax.grid(True)
        else:
            ax.set_xlim(1, 500)
            ax.set_ylim(0, 20)
            x_major_ticks = arange(0, 500, 10)
            x_minor_ticks = arange(0, 500, 2)
            y_major_ticks = arange(0, 20, 0.1)
            y_minor_ticks = arange(0, 20, 0.02)

            ax.set_xticks(x_major_ticks)
            ax.set_xticks(x_minor_ticks, minor=True)
            ax.set_yticks(y_major_ticks)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(True)

        """
        return fig


class Curves:

    def __init__(self, logger, prfs_d):
        self.plot_curves(logger, prfs_d)


    def plot_curves(logger, prfs_d):
        """

        @param logger:
        @param prefs_dict:
        """

        # df_1, df_2 = read_output(logger, prefs_dict, 0)

        df = read_csv('full_stats_2.csv', index_col=0)

        with PdfPages('output.pdf') as pdf:
            markers = ['o', 's', '*', '^', 'v', '<', '>']
            purity_level = 0.9
            alpha_v = 0.7
            size = 35
            margin = 0.1

            x_legend = np.arange(17)
            y_legend = np.empty(17)
            print y_legend
            y_legend = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1]

            t = np.arange(17)
            colors = ['b', 'g', 'r', 'c', 'm']

            filters = ['models/gauss_2.0_5x5.conv',
                       'models/tophat_2.0_3x3.conv']
            mags = [0, 1, 2]
            speeds = [0.001, 1, 3, 5, 10]

            # Looping over magnitudes
            for mag in mags:
                d_1 = df[df['mag'].isin([mag])]
                logger.debug('analysis for magnitude {} launched!'.format(mag))
                # Looping over speeds
                for speed in speeds:
                    logger.debug('analysis for speed {} launched!'.format(speed))
                    x_range = {}
                    y_range = {}
                    # Create figure
                    fig = plt.figure(figsize=(16.5, 11.7), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)
                    d_2 = d_1[d_1['speed'].isin([speed])]
                    # Looping over filters
                    for filt in filters:
                        logger.debug('filter {} starts!'.format(filt[7:]))
                        d_3 = d_2[d_2['filter'].isin([filt])]
                        f_com = d_3['completeness']
                        f_pur = d_3['purity']

                        ax.scatter(f_com, f_pur, c=t, s=size,
                                   marker=markers[filters.index(filt)],
                                   label='s {} f {}'.format(speed,
                                                            filt[7:-5]))

                        ax.plot(f_com, f_pur, color=colors[filters.index(filt)],
                                label='filter {}'.format(filt[7:-5]),
                                linewidth=1)

                        ax.plot([0.9, 0.9], [0, 1.1], color='b')
                        ax.plot([-1, 1.1], [0.9, 0.9], color='b')
                        ax.plot([0.8, 0.8], [0, 1.1], color='r')
                        ax.plot([-1, 1.1], [0.8, 0.8], color='r')

                    x_min = 0
                    x_max = 1.1
                    y_min = 0
                    y_max = 1.1
                    ax.axis([x_min, x_max, y_min, y_max])
                    ax.legend(bbox_to_anchor=(1.125, 1))

                    ax.set_xlabel('f_com')
                    ax.set_ylabel('f_pur')
                    mags_t = ['24-25', '25-26', '26-27']
                    fig.suptitle('magnitude {} speed {}'.format(mags_t[mag],
                                                                speed))

                    x_major_ticks = np.arange(x_min, x_max, 0.10)
                    x_minor_ticks = np.arange(x_min, x_max, 0.02)
                    y_major_ticks = np.arange(y_min, y_max, 0.10)
                    y_minor_ticks = np.arange(y_min, y_max, 0.02)

                    ax.set_xticks(x_major_ticks)
                    ax.set_xticks(x_minor_ticks, minor=True)
                    ax.set_yticks(y_major_ticks)
                    ax.set_yticks(y_minor_ticks, minor=True)

                    ax.grid(which='minor', alpha=0.5)
                    ax.grid(which='major', alpha=0.9)

                    pdf.savefig()

            logger.debug('n_seas vs.')
            for mag in mags:
                d_1 = df[df['mag'].isin([mag])]
                logger.debug('analysis for magnitude {} launched!'.format(mag))
                # Looping over speeds
                for speed in speeds:
                    logger.debug('analysis for speed {} launched!'.format(speed))
                    x_range = {}
                    y_range = {}
                    # Create figure
                    fig = plt.figure(figsize=(16.5, 11.7), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)
                    d_2 = d_1[d_1['speed'].isin([speed])]
                    # Looping over filters
                    for filt in filters:
                        logger.debug('filter {} starts!'.format(filt[7:]))
                        d_3 = d_2[d_2['filter'].isin([filt])]
                        det = d_3['detected']
                        ok = d_3['right']
                        # x range analysis
                        x_range[filters.index(filt)] = det.tolist()
                        # y range analysis
                        y_range[filters.index(filt)] = ok.tolist()
                        ax.scatter(det, ok, c=t, s=size,
                                   marker=markers[filters.index(filt)],
                                   label='s {} f {}'.format(speed,
                                                            filt[7:-5]))

                        ax.plot(det, ok, color=colors[filters.index(filt)],
                                label='filter {}'.format(filt[7:-5]),
                                linewidth=1)

                    x_min = 0
                    x_max = 425
                    y_min = 0
                    y_max = 120

                    ax.plot([0, 125], [0, 125], color='b',
                            label='f_pur=100')
                    ax.plot([0, 125 / 0.9], [0, 125], color='g',
                            label='f_pur=0.9')
                    ax.plot([0, 125 / 0.8], [0, 125], color='r',
                            label='f_pur=0.8')
                    ax.plot([0, 125 / 0.7], [0, 125], color='c',
                            label='f_pur=0.7')
                    ax.plot([0, 125 / 0.6], [0, 125], color='m',
                            label='f_pur=0.7')
                    ax.plot([100, 100], [0, 125], color='m')

                    ax.axis([x_min, x_max, y_min, y_max])
                    ax.legend(bbox_to_anchor=(1.125, 1))

                    ax.set_xlabel('n_meas')
                    ax.set_ylabel('n_seas')
                    mags_t = ['24-25', '25-26', '26-27']
                    fig.suptitle('magnitude {} speed {}'.format(mags_t[mag],
                                                                speed))

                    x_major_ticks = np.arange(x_min, x_max, 25)
                    x_minor_ticks = np.arange(x_min, x_max, 5)
                    y_major_ticks = np.arange(y_min, y_max, 25)
                    y_minor_ticks = np.arange(y_min, y_max, 5)

                    ax.set_xticks(x_major_ticks)
                    ax.set_xticks(x_minor_ticks, minor=True)
                    ax.set_yticks(y_major_ticks)
                    ax.set_yticks(y_minor_ticks, minor=True)

                    ax.grid(which='minor', alpha=0.5)
                    ax.grid(which='major', alpha=0.9)

                    pdf.savefig()

            """
            logger.debug('legend plotting!')
            fig = plt.figure(figsize=(16.5, 11.7), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            l_threshold = [1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45,
                           1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625,
                           1.65, 1.675, 1.7]

            for thres in l_threshold:
                print t[l_threshold.index(thres)]
                ax.scatter(x_legend[l_threshold.index(thres)],
                           y_legend[l_threshold.index(thres)],
                           color=t[l_threshold.index(thres)], s=size,
                           marker=markers[filters.index(filt)],
                           label='threshold {}'.format(thres))

            ax.legend(bbox_to_anchor=(1.125, 1))
            fig.suptitle('threshold legend')

            pdf.savefig()
            """ 