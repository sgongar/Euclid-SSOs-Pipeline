#!/usr/bin/python
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
from pandas import read_csv

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"

"""
catalog_n = input('Catalog number: ')
i_alpha = input('ALPHA_J2000: ')
i_delta = input('DELTA_J2000: ')

o_file = fits.open('/home/sgongora/Documents/CarpetaCompartida/full_10_1.1_0.5_0.04_20-21_1.cat')
o_cat = Table(o_file[2].data).to_pandas()

tolerance = 0.0002

o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

print(o_df)

"""
# 0.3
# sources_list = [28193, 33556, 23878, 42982, 29232, 3892, 27867, 34952, 35368,
#                 19912, 30751, 27685, 35752, 25018, 9818, 36145, 45461, 22077,
#                 26159, 7933, 41378, 39627]
# 0.03
sources_list = [39209, 2685, 42050, 14999, 19603, 21143, 21977, 10403, 20254,
                3833, 26337, 37751, 29414, 12291, 28938, 20752, 29808]

dir_ = '/home/sgongora/Documents/CarpetaCompartida'
file_ = 'filt_10_1.1_0.5_0.04_20-21_3.csv'
o_cat = read_csv('{}/{}'.format(dir_, file_), index_col=0)

b_image_l = []

for idx, source_ in enumerate(sources_list):
    o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
    # print(o_df['CLASS_STAR'])
    for b_image in o_df['CLASS_STAR'].tolist():
        b_image_l.append(b_image)
    # print(o_df['B_IMAGE'])


from numpy import mean, std, median

print(mean(b_image_l))
print(std(b_image_l))
print(min(b_image_l))
print(max(b_image_l))
