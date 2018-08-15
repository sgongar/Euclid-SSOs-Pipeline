#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from cats_management import get_max_mag


def plot_tables(logger, prfs_d, configurations, conf_num):
    """

    @param logger: a logger object
    @param prfs_d: a dictionary with all configuration information
    @param configurations:
    @param conf_num:

    @return fig:
    """
    speeds = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1,
              0.3, 0.5, 1, 3, 5, 10, 30, 50, 100, 130, 150]

    fig, ax = plt.subplots(figsize=(11.69, 8.27), dpi=100)

    # Hide axes
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.autoscale_view('tight')

    max_values = get_max_mag(logger, prfs_d, conf_num, speeds)
    std = 0.5

    table_speeds = []

    # Table from Ed Smith answer
    for speed in range(len(speeds)):
        table_speeds_part = []
        table_speeds_part.insert(0, conf_num)

        table_speeds_part.insert(1, speeds[speed])

        detect_minarea = configurations[conf_num][0]
        table_speeds_part.insert(2, detect_minarea)

        detect_thresh = configurations[conf_num][1]
        table_speeds_part.insert(3, detect_thresh)

        analysis_thresh = configurations[conf_num][2]
        table_speeds_part.insert(4, analysis_thresh)

        deblend_nthresh = configurations[conf_num][3]
        table_speeds_part.insert(5, deblend_nthresh)

        deblend_mincount = configurations[conf_num][4]
        table_speeds_part.insert(6, deblend_mincount)

        table_speeds_part.insert(7, max_values[speed][0])

        std = get_std(logger, prfs_d, conf_num)
        table_speeds_part.insert(8, std)

        table_speeds.append(table_speeds_part)

    colums_label = ("analysis", "speed", "detect_minarea",
                    "detect_thresh", "analysis_thresh",
                    "deblend_nthresh", "deblend_mincount",
                    "Max. magnitude", "Standard deviation")
    ax.table(cellText=table_speeds, colLabels=colums_label,
             loc='center')

    return fig
