#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Plots pm and magnitude against A/B/elongation

Versions:
- 0.1 Initial release

Todo:
    * Improve usability

"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv


class Plot3D:

    def __init__(self, argv):
        self.arguments = argv
        stars_df = read_csv('stars_df.csv', index_col=0)
        self.stars_d = stars_df.to_dict()
        galaxies_df = read_csv('galaxies_df.csv', index_col=0)
        self.galaxies_d = galaxies_df.to_dict()
        ssos_df = read_csv('ssos_df.csv', index_col=0)
        self.ssos_d = ssos_df.to_dict()

        self.data_d = {}
        self.manage_dict()
        self.plot_3d_figure()

    def manage_dict(self):
        """

        :return:
        """
        stars_mag = []
        for key_ in self.stars_d['stars_mag'].keys():
            stars_mag.append(self.stars_d['stars_mag'][key_])
        stars_pm = []
        for key_ in self.stars_d['stars_pm'].keys():
            stars_pm.append(self.stars_d['stars_pm'][key_])
        stars_elongation = []
        for key_ in self.stars_d['stars_elongation'].keys():
            stars_elongation.append(self.stars_d['stars_elongation'][key_])

        self.data_d['x_stars'] = np.array(stars_mag)
        self.data_d['y_stars'] = np.array(stars_pm)
        self.data_d['z_stars'] = np.array(stars_elongation)

        # Galaxies
        galaxies_mag = []
        for key_ in self.galaxies_d['galaxies_mag'].keys():
            galaxies_mag.append(self.galaxies_d['galaxies_mag'][key_])
        galaxies_pm = []
        for key_ in self.galaxies_d['galaxies_pm'].keys():
            galaxies_pm.append(self.galaxies_d['galaxies_pm'][key_])
        galaxies_a = []
        for key_ in self.galaxies_d['galaxies_a'].keys():
            galaxies_a.append(self.galaxies_d['galaxies_a'][key_])
        galaxies_b = []
        for key_ in self.galaxies_d['galaxies_b'].keys():
            galaxies_b.append(self.galaxies_d['galaxies_b'][key_])
        galaxies_elongation = []
        for key_ in self.galaxies_d['galaxies_elongation'].keys():
            galaxies_elongation.append(self.galaxies_d['galaxies_elongation'][key_])

        self.data_d['x_galaxies'] = np.array(galaxies_mag)
        self.data_d['y_galaxies'] = np.array(galaxies_pm)
        self.data_d['z_galaxies'] = np.array(galaxies_elongation)

        # SSOs
        ssos_mag = []
        for key_ in self.ssos_d['ssos_mag'].keys():
            ssos_mag.append(self.ssos_d['ssos_mag'][key_])
        ssos_pm = []
        for key_ in self.ssos_d['ssos_pm'].keys():
            ssos_pm.append(self.ssos_d['ssos_pm'][key_])
        ssos_a = []
        for key_ in self.ssos_d['ssos_a'].keys():
            ssos_a.append(self.ssos_d['ssos_a'][key_])
        ssos_b = []
        for key_ in self.ssos_d['ssos_b'].keys():
            ssos_b.append(self.ssos_d['ssos_b'][key_])
        ssos_elongation = []
        for key_ in self.ssos_d['ssos_elongation'].keys():
            ssos_elongation.append(self.ssos_d['ssos_elongation'][key_])

        self.data_d['x_ssos'] = np.array(ssos_mag)
        self.data_d['y_ssos'] = np.array(ssos_pm)
        self.data_d['z_ssos'] = np.array(ssos_elongation)

    def plot_3d_figure(self):
        """

        :return:
        """
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlabel('X - Magnitude')
        ax.set_ylabel('Y - Proper Motion "/h')
        ax.set_ylim3d([0, 5])
        ax.set_zlabel('Z - Elongation (A/B)')
        if self.arguments[1] == '-all':
            ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                         self.data_d['z_stars'], c='b')
            ax.scatter3D(self.data_d['x_galaxies'], self.data_d['y_galaxies'],
                         self.data_d['z_galaxies'], c='r')
            ax.scatter3D(self.data_d['x_ssos'], self.data_d['y_ssos'],
                         self.data_d['z_ssos'], c='g')
        elif self.arguments[1] == '-stars':
            if self.arguments[2] == '-galaxies':
                ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                             self.data_d['z_stars'], c='b')
                ax.scatter3D(self.data_d['x_galaxies'],
                             self.data_d['y_galaxies'],
                             self.data_d['z_galaxies'], c='r')
            elif self.arguments[2] == '-ssos':
                ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                             self.data_d['z_stars'], c='b')
                ax.scatter3D(self.data_d['x_ssos'], self.data_d['y_ssos'],
                             self.data_d['z_ssos'], c='g')
            elif self.arguments[2] == '-':
                ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                             self.data_d['z_stars'], c='b')
        elif self.arguments[1] == '-galaxies':
            if self.arguments[2] == '-stars':
                ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                             self.data_d['z_stars'], c='b')
                ax.scatter3D(self.data_d['x_galaxies'],
                             self.data_d['y_galaxies'],
                             self.data_d['z_galaxies'], c='r')
            elif self.arguments[2] == '-ssos':
                ax.scatter3D(self.data_d['x_galaxies'],
                             self.data_d['y_galaxies'],
                             self.data_d['z_galaxies'], c='r')
                ax.scatter3D(self.data_d['x_ssos'], self.data_d['y_ssos'],
                             self.data_d['z_ssos'], c='g')
            elif self.arguments[2] == '-':
                ax.scatter3D(self.data_d['x_galaxies'],
                             self.data_d['y_galaxies'],
                             self.data_d['z_galaxies'], c='r')
        elif self.arguments[1] == '-ssos':
            if self.arguments[2] == '-stars':
                ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                             self.data_d['z_stars'], c='b')
                ax.scatter3D(self.data_d['x_ssos'],
                             self.data_d['y_ssos'],
                             self.data_d['z_ssos'], c='g')
            elif self.arguments[2] == '-galaxies':
                ax.scatter3D(self.data_d['x_galaxies'],
                             self.data_d['y_galaxies'],
                             self.data_d['z_galaxies'], c='r')
                ax.scatter3D(self.data_d['x_ssos'], self.data_d['y_ssos'],
                             self.data_d['z_ssos'], c='g')
            elif self.arguments[2] == '-':
                ax.scatter3D(self.data_d['x_ssos'],
                             self.data_d['y_ssos'],
                             self.data_d['z_ssos'], c='g')
        else:
            raise Exception
        plt.show()


if __name__ == "__main__":
    from sys import argv

    test = Plot3D(argv)