#!/usr/bin/python
# -*- coding: utf-8 -*-

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import array, arange, mean, std
from scipy import stats

# TODO read automatically
first_galaxy_x = [99.4183578491, 199.506011963, 299.507415771, 399.524291992,
                  499.418365479, 599.378479004, 699.528015137, 799.580505371,
                  899.462768555, 999.611877441, 1099.47290039, 1199.52978516]
first_galaxy_y = [599.74609375, 599.810424805, 599.979187012, 599.985656738,
                  600.090393066, 600.268676758, 600.105102539, 600.241516113,
                  600.316040039, 600.132446289, 600.706481934, 600.829528809]

second_galaxy_x = [99.8847122192, 199.859207153, 299.961517334, 399.94732666,
                   499.832000732, 599.973388672, 699.793945312, 799.907531738,
                   899.992126465, 999.84765625, 1099.8548584, 1200.10827637]
second_galaxy_y = [1100.11486816, 1100.30175781, 1100.4263916, 1100.54382324,
                   1100.38085938, 1100.63879395, 1100.78198242, 1100.67871094,
                   1100.88256836, 1100.88110352, 1101.11779785, 1101.13586426]

third_galaxy_y = [1600.81323242, 1600.99536133, 1600.80651855, 1600.92785645,
                  1600.91833496, 1601.06518555, 1601.17004395, 1601.09436035,
                  1601.24560547, 1601.31970215, 1601.79345703, 1601.8503418]
third_galaxy_x = [99.3178939819, 199.40536499, 299.470031738, 399.261077881,
                  499.432373047, 599.323669434, 699.460083008, 799.548156738,
                  899.383117676, 999.465026855, 1099.22949219, 1199.37475586]


def get_stats():
    """

    :return:
    """
    d_stats = {}

    # First galaxy stats
    first_input_galaxy_y_a = arange(599.5, 600.6, 0.1)
    first_input_galaxy_y_b = arange(599.6, 600.7, 0.1)
    first_input_galaxy_y_c = arange(599.7, 600.9, 0.1)

    difference_a_l = []
    for idx, output_a in enumerate(first_galaxy_y):
        # print(output_a)
        difference_a = abs(first_input_galaxy_y_a[idx] - output_a)
        difference_a_l.append(difference_a)
    difference_b_l = []
    for idx, output_b in enumerate(first_galaxy_y):
        # print(output_a)
        difference_b = abs(first_input_galaxy_y_b[idx] - output_b)
        difference_b_l.append(difference_b)
    difference_c_l = []
    for idx, output_c in enumerate(first_galaxy_y):
        # print(output_a)
        difference_c = abs(first_input_galaxy_y_c[idx] - output_c)
        difference_c_l.append(difference_c)

    d_stats['mean_1st_a'] = mean(difference_a_l)
    d_stats['std_1st_a'] = std(difference_a_l)
    d_stats['mean_1st_b'] = mean(difference_b_l)
    d_stats['std_1st_b'] = std(difference_b_l)
    d_stats['mean_1st_c'] = mean(difference_c_l)
    d_stats['std_1st_c'] = std(difference_c_l)

    # Second galaxy stats
    second_input_galaxy_y_a = arange(1100.1, 1101.2, 0.1)
    second_input_galaxy_y_b = arange(1100.15, 1101.35, 0.1)
    second_input_galaxy_y_c = arange(1100.2, 1101.4, 0.1)

    difference_a_l = []
    for idx, output_a in enumerate(second_galaxy_y):
        # print(output_a)
        difference_a = abs(second_input_galaxy_y_a[idx] - output_a)
        difference_a_l.append(difference_a)
    difference_b_l = []
    for idx, output_b in enumerate(second_galaxy_y):
        # print(output_a)
        difference_b = abs(second_input_galaxy_y_b[idx] - output_b)
        difference_b_l.append(difference_b)
    difference_c_l = []
    for idx, output_c in enumerate(second_galaxy_y):
        # print(output_a)
        difference_c = abs(second_input_galaxy_y_c[idx] - output_c)
        difference_c_l.append(difference_c)

    d_stats['mean_2nd_a'] = mean(difference_a_l)
    d_stats['std_2nd_a'] = std(difference_a_l)
    d_stats['mean_2nd_b'] = mean(difference_b_l)
    d_stats['std_2nd_b'] = std(difference_b_l)
    d_stats['mean_2nd_c'] = mean(difference_c_l)
    d_stats['std_2nd_c'] = std(difference_c_l)

    """
    # Third galaxy stats
    third_input_galaxy_y_a = arange(1100.1, 1101.2, 0.1)
    second_input_galaxy_y_b = arange(1100.15, 1101.35, 0.1)
    second_input_galaxy_y_c = arange(1100.2, 1101.4, 0.1)

    difference_a_l = []
    for idx, output_a in enumerate(second_galaxy_y):
        # print(output_a)
        difference_a = abs(second_input_galaxy_y_a[idx] - output_a)
        difference_a_l.append(difference_a)
    difference_b_l = []
    for idx, output_b in enumerate(second_galaxy_y):
        # print(output_a)
        difference_b = abs(second_input_galaxy_y_b[idx] - output_b)
        difference_b_l.append(difference_b)
    difference_c_l = []
    for idx, output_c in enumerate(second_galaxy_y):
        # print(output_a)
        difference_c = abs(second_input_galaxy_y_c[idx] - output_c)
        difference_c_l.append(difference_c)

    d_stats['mean_2nd_a'] = mean(difference_a_l)
    d_stats['std_2nd_a'] = std(difference_a_l)
    d_stats['mean_2nd_b'] = mean(difference_b_l)
    d_stats['std_2nd_b'] = std(difference_b_l)
    d_stats['mean_2nd_c'] = mean(difference_c_l)
    d_stats['std_2nd_c'] = std(difference_c_l)
    """

    return d_stats


def plot(d_stats):
    """

    :param d_stats:
    :return:
    """
    # TODO generate automatically
    # x1 ticks
    x1_ticks = {'major_ticks': arange(100, 1300, 100)}
    y1_ticks = {'major_ticks': arange(599.5, 601, 0.1),
                'minor_ticks': arange(599.5, 600.92, 0.02)}

    # x2 ticks
    x2_ticks = {'major_ticks': arange(100, 1300, 100)}
    y2_ticks = {'major_ticks': arange(1099.9, 1101.4, 0.1),
                'minor_ticks': arange(1099.9, 1101.32, 0.02)}

    # x3 ticks
    x3_ticks = {'major_ticks': arange(100, 1300, 100)}
    y3_ticks = {'major_ticks': arange(1600.5, 1602.02, 0.1),
                'minor_ticks': arange(1600.5, 1602.02, 0.02)}

    plot_size = [16.53, 11.69]
    plot_dpi = 100

    pdf_name = 'movement.pdf'

    with PdfPages(pdf_name) as pdf:
        # First galaxy curve
        fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
        ax_1 = fig.add_subplot(1, 1, 1)
        ax_1.set_title('First galaxy mean-std')

        step = 0.005
        n_estimated_l = arange(599, 600, step)

        estimated_means = []
        estimated_stds = []
        for idx_n, n_ in enumerate(n_estimated_l):
            y_estimated_l = []
            for idx_x, x_ in enumerate(first_galaxy_x ):
                y_estimated = x_ * 0.001 + n_
                y_difference = abs(first_galaxy_y[idx_x] - y_estimated)
                y_estimated_l.append(y_difference)
            estimated_means.append(mean(y_estimated_l))
            estimated_stds.append(std(y_estimated_l))

        # Minimum value for estimated means
        idx_mean = estimated_means.index(min(estimated_means))
        n_min_mean = n_estimated_l[idx_mean]
        # Minimum value for estimated standard deviations
        idx_std = estimated_stds.index(min(estimated_stds))
        n_min_std = n_estimated_l[idx_std]

        ax_1.plot(n_estimated_l, estimated_means, label='mean value')
        ax_1.plot(n_estimated_l, estimated_stds, label='std value')

        ax_1.set_xticks([599.0, 599.1, 599.2, 599.3, 599.4, 599.5, 599.6,
                         599.7, 599.8, 599.9, 600], minor=False)
        ax_1.set_xticks(arange(598.95, 600.05, 0.01), minor=True)
        ax_1.set_yticks(arange(0.1, 0.6, 0.1), minor=False)
        ax_1.set_yticks(arange(0.05, 0.55, 0.01), minor=True)

        # Labels
        x_label = 'y-value (pixels) for x=0\n'
        x_label_major_step = 'major step size 0.1 pix\n'
        x_label_minor_step = 'minor step size 0.01 pix'
        ax_1.set_xlabel('{}{}{}'.format(x_label, x_label_major_step,
                                        x_label_minor_step))
        y_major_step = 'major step size 0.1\n'
        y_minor_step = 'minor step size 0.01'
        ax_1.set_ylabel('{}{}'.format(y_major_step, y_minor_step))

        # Grid configuration
        ax_1.grid(b=True, which='major', linestyle='-', linewidth=3)
        ax_1.grid(b=True, which='minor', linestyle='-', linewidth=1)

        ax_1.legend(loc=4)
        pdf.savefig()  # saves figure
        plt.clf()  # clear current figure
        plt.close(fig)  # removes figure

        # First galaxy
        fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
        ax_2 = fig.add_subplot(1, 1, 1)
        ax_2.set_title('First galaxy positions')
        ax_2.scatter(first_galaxy_x, first_galaxy_y)
        # Fitting for minimum mean
        ax_2.plot([100, 1200], [(0.001 * 100) + n_min_mean,
                                (0.001 * 1200) + n_min_mean],
                  label='minimum mean')
        # Fitting for minimum std
        ax_2.plot([100, 1200], [(0.001 * 100) + n_min_std,
                                (0.001 * 1200) + n_min_std],
                  label='minimum std')

        # Ticks positions
        ax_2.set_xticks(x1_ticks['major_ticks'], minor=False)
        ax_2.set_yticks(y1_ticks['major_ticks'], minor=False)
        ax_2.set_yticks(y1_ticks['minor_ticks'], minor=True)

        # Labels
        x_label_major_step = 'major step size 100 pix'
        ax_2.set_xlabel('y-axis pixels\n{}'.format(x_label_major_step))
        y_label_major_step = 'major step size 0.1 pix\n'
        y_label_minor_step = 'minor step size 0.02 pix'
        ax_2.set_ylabel('x-axis pixels\n{}{}'.format(y_label_major_step,
                                                     y_label_minor_step))

        # Grid configuration
        ax_2.grid(b=True, which='major', linestyle='-', linewidth=2)
        ax_2.grid(b=True, which='minor', linestyle='--', linewidth=1)

        ax_2.legend(loc=4)
        pdf.savefig()  # saves figure
        plt.clf()  # clear current figure
        plt.close(fig)  # removes figure

        # Second galaxy curve
        fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
        ax_3 = fig.add_subplot(1, 1, 1)
        ax_3.set_title('Second galaxy mean-std')

        step = 0.005
        n_estimated_l = arange(1099.8, 1100.1, step)

        estimated_means = []
        estimated_stds = []
        for idx_n, n_ in enumerate(n_estimated_l):
            y_estimated_l = []
            for idx_x, x_ in enumerate(second_galaxy_x ):
                y_estimated = x_ * 0.001 + n_
                y_difference = abs(second_galaxy_y[idx_x] - y_estimated)
                y_estimated_l.append(y_difference)
            estimated_means.append(mean(y_estimated_l))
            estimated_stds.append(std(y_estimated_l))

        # Minimum value for estimated means
        idx_mean = estimated_means.index(min(estimated_means))
        n_min_mean = n_estimated_l[idx_mean]
        # Minimum value for estimated standard deviations
        idx_std = estimated_stds.index(min(estimated_stds))
        n_min_std = n_estimated_l[idx_std]

        ax_3.plot(n_estimated_l, estimated_means, label='mean value', c='r')
        ax_3.plot(n_estimated_l, estimated_stds, label='std value', c='y')

        ax_3.set_xticks([1099.8, 1099.85, 1099.9, 1099.95, 1100, 1100.05,
                         1100.1], minor=False)
        ax_3.set_xticks(arange(1099.785, 1100.11, 0.005), minor=True)
        ax_3.set_yticks(arange(0.04, 0.22, 0.02), minor=False)
        ax_3.set_yticks(arange(0.04, 0.22, 0.005), minor=True)

        ax_3.set_xticklabels(['1099.8', '1099.85', '1099.9', '1099.95', '1100',
                              '1100.05', '1100.1'])

        # Labels
        x_label = 'y-value (pixels) for x=0\n'
        x_label_major_step = 'major step size 0.05 pix\n'
        x_label_minor_step = 'minor step size 0.005 pix'
        ax_3.set_xlabel('{}{}{}'.format(x_label, x_label_major_step,
                                        x_label_minor_step))
        y_major_step = 'major step size 0.02\n'
        y_minor_step = 'minor step size 0.005'
        ax_3.set_ylabel('{}{}'.format(y_major_step, y_minor_step))

        # Grid configuration
        ax_3.grid(b=True, which='major', linestyle='-', linewidth=3)
        ax_3.grid(b=True, which='minor', linestyle='-', linewidth=1)

        ax_3.legend(loc=4)
        pdf.savefig()  # saves figure
        plt.clf()  # clear current figure
        plt.close(fig)  # removes figure

        # Second galaxy
        fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
        ax_4 = fig.add_subplot(1, 1, 1)
        ax_4.set_title('Second galaxy positions')
        ax_4.scatter(second_galaxy_x, second_galaxy_y)
        # Fitting for minimum mean
        ax_4.plot([100, 1200], [(0.001 * 100) + n_min_mean,
                                (0.001 * 1200) + n_min_mean],
                  label='minimum mean', c='r')
        # Fitting for minimum std
        ax_4.plot([100, 1200], [(0.001 * 100) + n_min_std,
                                (0.001 * 1200) + n_min_std],
                  label='minimum std', c='y')

        ax_4.set_xticks(x2_ticks['major_ticks'], minor=False)
        ax_4.set_yticks(y2_ticks['major_ticks'], minor=False)
        ax_4.set_yticks(y2_ticks['minor_ticks'], minor=True)

        # Labels
        x_label_major_step = 'major step size 100 pix'
        ax_4.set_xlabel('y-axis pixels\n{}'.format(x_label_major_step))
        y_label_major_step = 'major step size 0.1 pix\n'
        y_label_minor_step = 'minor step size 0.02 pix'
        ax_4.set_ylabel('x-axis pixels\n{}{}'.format(y_label_major_step,
                                                     y_label_minor_step))

        ax_4.grid(b=True, which='major', linestyle='-', linewidth=2)
        ax_4.grid(b=True, which='minor', linestyle='--', linewidth=1)

        ax_4.legend(loc=4)

        pdf.savefig()  # saves figure
        plt.clf()  # clear current figure
        plt.close(fig)  # removes figure

        # Third galaxy curve
        fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
        ax_5 = fig.add_subplot(1, 1, 1)
        ax_5.set_title('Third galaxy mean-std')

        step = 0.005
        n_estimated_l = arange(1600.4, 1600.7, step)

        estimated_means = []
        estimated_stds = []
        for idx_n, n_ in enumerate(n_estimated_l):
            y_estimated_l = []
            for idx_x, x_ in enumerate(third_galaxy_x):
                y_estimated = x_ * 0.001 + n_
                y_difference = abs(third_galaxy_y[idx_x] - y_estimated)
                y_estimated_l.append(y_difference)
            estimated_means.append(mean(y_estimated_l))
            estimated_stds.append(std(y_estimated_l))

        # Minimum value for estimated means
        idx_mean = estimated_means.index(min(estimated_means))
        n_min_mean = n_estimated_l[idx_mean]
        # Minimum value for estimated standard deviations
        idx_std = estimated_stds.index(min(estimated_stds))
        n_min_std = n_estimated_l[idx_std]

        ax_5.plot(n_estimated_l, estimated_means, label='mean value', c='g')
        ax_5.plot(n_estimated_l, estimated_stds, label='std value', c='k')

        ax_5.set_xticks([1600.4, 1600.45, 1600.5, 1600.55, 1600.6,
                         1600.65, 1600.7], minor=False)
        ax_5.set_xticks(arange(1600.385, 1600.71, 0.005), minor=True)
        ax_5.set_yticks(arange(0.06, 0.20, 0.02), minor=False)
        ax_5.set_yticks(arange(0.06, 0.20, 0.005), minor=True)

        ax_5.set_xticklabels(['1600.4', '1600.45', '1600.5', '1600.55',
                              '1600.6', '1600.65', '1600.7'])

        # Labels
        x_label = 'y-value (pixels) for x=0\n'
        x_label_major_step = 'major step size 0.05 pix\n'
        x_label_minor_step = 'minor step size 0.005 pix'
        ax_5.set_xlabel('{}{}{}'.format(x_label, x_label_major_step,
                                        x_label_minor_step))
        y_major_step = 'major step size 0.02\n'
        y_minor_step = 'minor step size 0.005'
        ax_5.set_ylabel('{}{}'.format(y_major_step, y_minor_step))

        # Grid configuration
        ax_5.grid(b=True, which='major', linestyle='-', linewidth=3)
        ax_5.grid(b=True, which='minor', linestyle='-', linewidth=1)

        ax_5.legend(loc=4)
        pdf.savefig()  # saves figure
        plt.clf()  # clear current figure
        plt.close(fig)  # removes figure

        # Third galaxy
        fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
        ax_6 = fig.add_subplot(1, 1, 1)
        ax_6.set_title('Third galaxy positions')

        ax_6.scatter(third_galaxy_x, third_galaxy_y)
        # Fitting for minimum mean
        ax_6.plot([100, 1200], [(0.001 * 100) + n_min_mean,
                                (0.001 * 1200) + n_min_mean],
                  label='minimum mean', c='g')
        # Fitting for minimum std
        ax_6.plot([100, 1200], [(0.001 * 100) + n_min_std,
                                (0.001 * 1200) + n_min_std],
                  label='minimum std', c='k')

        ax_6.set_xticks(x3_ticks['major_ticks'], minor=False)
        ax_6.set_yticks(y3_ticks['major_ticks'], minor=False)
        ax_6.set_yticks(y3_ticks['minor_ticks'], minor=True)

        # Labels
        x_label_major_step = 'major step size 100 pix'
        ax_6.set_xlabel('y-axis pixels\n{}'.format(x_label_major_step))
        y_label_major_step = 'major step size 0.1 pix\n'
        y_label_minor_step = 'minor step size 0.02 pix'
        ax_6.set_ylabel('x-axis pixels\n{}{}'.format(y_label_major_step,
                                                     y_label_minor_step))

        ax_6.grid(b=True, which='major', linestyle='-', linewidth=2)
        ax_6.grid(b=True, which='minor', linestyle='--', linewidth=1)

        ax_6.legend(loc=4)
        pdf.savefig()
        plt.clf()  # clear current figure
        plt.close(fig)


if __name__ == "__main__":
    data_stats = get_stats()
    plot(data_stats)
