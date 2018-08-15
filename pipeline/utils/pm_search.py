#!/usr/bin/python
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from math import sqrt, cos, radians, hypot
from numpy import genfromtxt, float64
from pandas import concat, read_csv, Series
from itertools import groupby
from operator import itemgetter

# """
# source_n = input('Catalog number: ')
#
# f_file = read_csv('/home/sgongora/Documents/CarpetaCompartida/filt_10_1.2_2.5_0.64_20-21_4.csv', index_col=0)

# f_pm = f_file[f_file['SOURCE_NUMBER'].isin([source_n])]
#
# alpha_1 = float(f_pm['ALPHA_J2000'].iloc[0])
# alpha_2 = float(f_pm['ALPHA_J2000'].iloc[1])
# alpha_3 = float(f_pm['ALPHA_J2000'].iloc[2])
#
# delta_1 = float(f_pm['DELTA_J2000'].iloc[0])
# delta_2 = float(f_pm['DELTA_J2000'].iloc[1])
# delta_3 = float(f_pm['DELTA_J2000'].iloc[2])
#
# m_32 = hypot(alpha_3 - alpha_2, delta_3 - delta_2)
# m_21 = hypot(alpha_2 - alpha_1, delta_2 - delta_1)
#
# pm_21 = (3600*m_21)/1003
# pm_21 = pm_21*3600
# pm_32 = (3600*m_21)/1003
# pm_32 = pm_32*3600
#
# f_pm = f_pm['PM'].iloc[0]
#
# print pm_21, pm_32, f_pm
# """
"""
source_n = input('Catalog number: ')

f_file = read_csv('/mnt/e/Documentos/CarpetaCompartida/input_sources.csv', index_col=0)

f_pm = f_file[f_file['source'].isin([source_n])]

alpha_1 = float(f_pm['alpha_j2000'].iloc[0])
alpha_2 = float(f_pm['alpha_j2000'].iloc[1])
alpha_3 = float(f_pm['alpha_j2000'].iloc[2])
alpha_4 = float(f_pm['alpha_j2000'].iloc[3])

delta_1 = float(f_pm['delta_j2000'].iloc[0])
delta_2 = float(f_pm['delta_j2000'].iloc[1])
delta_3 = float(f_pm['delta_j2000'].iloc[2])
delta_4 = float(f_pm['delta_j2000'].iloc[3])

m_43 = hypot(alpha_4 - alpha_3, delta_4 - delta_3)
m_32 = hypot(alpha_3 - alpha_2, delta_3 - delta_2)
m_21 = hypot(alpha_2 - alpha_1, delta_2 - delta_1)

pm_21 = (3600*m_21)/1003
pm_21 = pm_21*3600

pm_32 = (3600*m_32)/1003
pm_32 = pm_32*3600

pm_43 = (3600*m_43)/1003
pm_43 = pm_43*3600

#
# f_pm = f_pm['PM'].iloc[0]
#
print pm_21, pm_32, pm_43
"""
#
# input_ref = '/media/sf_CarpetaCompartida/luca_data/v14/Catalogs'
# first_sso = 134895
#
# input_d = {}
# for d in range(1, 5, 1):
#     input_d[d] = '{}/Cat_20-21_d{}.dat'.format(input_ref, d)
#
# for dither_ in input_d.keys():
#     catalog = genfromtxt(input_d[dither_])
#
#     list_x = catalog[:, 0]
#     list_y = catalog[:, 1]
#     list_mag = catalog[:, 2]
#     list_pm = catalog[:, 3]
#
#     x_values = []
#     y_values = []
#
#     speed_0_0003 = range(first_sso, 137446, 75)
#     speed_0_001 = range(first_sso + 10, 137456, 75)
#     speed_0_003 = range(first_sso + 20, 137466, 75)
#     speed_0_01 = range(first_sso + 30, 137476, 75)
#     speed_0_03 = range(first_sso + 40, 137486, 75)
#     speed_0_1 = range(first_sso + 50, 137496, 75)
#     speed_0_3 = range(first_sso + 60, 137506, 75)
#
#     speed_3 = range(first_sso + 67, 137513, 75)
#     print speed_3
#     speed_10 = range(first_sso + 68, 137514, 75)
#     speed_30 = range(first_sso + 69, 137515, 75)
#     speed_100 = range(first_sso + 70, 137516, 75)
#     speed_300 = range(first_sso + 71, 137517, 75)
#
#     for index in speed_0_0003:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.0003
#     for index in speed_0_001:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.001
#     for index in speed_0_003:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.003
#     for index in speed_0_01:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.01
#     for index in speed_0_03:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.03
#     for index in speed_0_1:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.1
#     for index in speed_0_3:
#         list_mag[index] = list_mag[index] - 2.5
#         list_pm[index] = 0.3
#     for index in speed_3:
#         list_pm[index] = list_pm[index] - 1000
#     for index in speed_10:
#         list_pm[index] = list_pm[index] - 1000
#     for index in speed_30:
#         list_pm[index] = list_pm[index] - 1000
#     for index in speed_100:
#         list_pm[index] = list_pm[index] - 1000
#     for index in speed_300:
#         list_pm[index] = list_pm[index] - 1000
#
#     indexes = (speed_0_0003 + speed_0_001 + speed_0_003 + speed_0_01 +
#                speed_0_03 + speed_0_1 + speed_0_3 + speed_3 + speed_10 +
#                speed_30 + speed_100 + speed_300)
#     indexes = sorted(indexes)
#
#     s1 = Series(list_x, name='X_IMAGE', dtype=float64)
#     s2 = Series(list_y, name='Y_IMAGE', dtype=float64)
#     s3 = Series(list_mag, name='MAG_VALUES', dtype=float64)
#     s4 = Series(list_pm, name='PM_INPUT', dtype=float64)
#
#     sources_df = concat([s1, s2, s3, s4], axis=1)
#     sources_df = sources_df.iloc[indexes, :]
#
#     fits_dir = '/media/sf_CarpetaCompartida/luca_data/v14/CCDs'
#     ccd_loc = 'mag_20-21_CCD_x0_y0_d1.fits'
#     fits_loc = '{}/{}'.format(fits_dir, ccd_loc)
#     hdulist = fits.open(fits_loc)
#     w = WCS(hdulist[0].header)
#
#     regions_list = []
#     for source_num in range(sources_df['X_IMAGE'].as_matrix().size):
#         x_value = sources_df['X_IMAGE'].as_matrix()[source_num]
#         y_value = sources_df['Y_IMAGE'].as_matrix()[source_num]
#         regions_list.append([x_value, y_value])
#
#     input_regions = w.all_pix2world(regions_list, 1)
#
#     print input_regions
#     """
#     alpha_list = []
#     delta_list = []
#     for idx, regions in enumerate(input_regions):
#         alpha_list.append(regions[0])
#         delta_list.append(regions[1])
#         x_values.append(regions_list[idx][0])
#         y_values.append(regions_list[idx][1])
#
#     fits_files_all = get_fits_d(dither=dither_)
#
#     fits_dict = {}
#     for fits_ in fits_files_all:
#         CCD = fits_[-13:-8]
#         fits_file = self.prfs_d['fits_dir'] + '/' + fits_
#         fits_dict[CCD] = get_fits_limits(fits_file)
#
#     i = 0
#     CCD_list = []
#     for alpha_, delta_ in zip(alpha_list, delta_list):
#         i += 1
#         flag = True
#         for key_ in fits_dict.keys():
#             below_ra = fits_dict[key_]['below_ra']
#             above_ra = fits_dict[key_]['above_ra']
#             below_dec = fits_dict[key_]['below_dec']
#             above_dec = fits_dict[key_]['above_dec']
#             alpha_comp = below_ra < alpha_ < above_ra
#             delta_comp = below_dec < delta_ < above_dec
#             if alpha_comp and delta_comp:
#                 CCD_list.append(key_)
#                 flag = False
#         if flag:
#             CCD_list.append('False')
#
#     # Creates a list for all sources
#     source_list = range(0, len(alpha_list), 1)
#     # Populates a list for all sources with the dither number
#     dither_list = []
#     for dither_idx in range(len(alpha_list)):
#         dither_list.append(dither_)
#
#     cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
#             ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
#             ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
#             ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
#             ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
#             ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
#             ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
#             ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
#             ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
#             ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
#             ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
#             ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]
#
#     cats_list = []
#     for dither_, CCD_ in zip(dither_list, CCD_list):
#         flag = True
#         for cat_ in cats:
#             if cat_[1] == dither_ and cat_[0] == CCD_:
#                 flag = False
#                 cats_list.append(cat_[2])
#
#         if flag:
#             cats_list.append(False)
#
#     # Creates a serie of Pandas Series
#     source = Series(source_list, name='source')
#
#     alpha_j2000 = Series(alpha_list, name='alpha_j2000')
#     delta_j2000 = Series(delta_list, name='delta_j2000')
#     mag = Series(sources_df['MAG_VALUES'].tolist(), name='mag_values')
#     pm = Series(sources_df['PM_INPUT'].tolist(), name='pm_values')
#     dither = Series(dither_list, name='dither_values')
#     CCD = Series(CCD_list, name='CCD')
#     cat = Series(cats_list, name='catalog')
#
#     if complete:
#         sources_df = concat([source, cat, alpha_j2000, delta_j2000,
#                              mag, pm, dither, CCD], axis=1)
#         sources_df = sources_df[~sources_df['CCD'].isin(['False'])]
#
#     else:
#         sources_df = concat([alpha_j2000, delta_j2000], axis=1)
#
#     dither_output = '{}/dither_{}'.format(self.prfs_d['dithers_out'],
#                                           dither_)
#
#     if save and not path.isfile(dither_output):
#         sources_df.to_csv(dither_output)
#
#     input_d[dither_] = sources_df

# source_i = input('Input number: ')
"""
i_file = read_csv('/home/sgongora/Documents/CarpetaCompartida/input_sources.csv', index_col=0)

unique_sources = list(set(i_file['source'].tolist()))

for source_i in unique_sources:
    # creates a dataframe for selected source
    i_pm = i_file[i_file['source'].isin([source_i])]

    # groups consecutive dithers
    tmp_groups = []
    data =  i_pm['dither_values'].tolist()
    for k, g in groupby(enumerate(data), lambda (i, x): i - x):
        tmp_groups.append(map(itemgetter(1), g))

    # look for a group of two consecutive dithers
    for group_ in tmp_groups:
        if len(group_) >= 2:
            group = group_

    i_pm_1 = i_pm[i_pm['dither_values'].isin([group[0]])]
    i_pm_2 = i_pm[i_pm['dither_values'].isin([group[1]])]

    i_alpha_1 = float(i_pm_1['alpha_j2000'].iloc[0])
    i_alpha_2 = float(i_pm_2['alpha_j2000'].iloc[0])

    i_delta_1 = float(i_pm_1['delta_j2000'].iloc[0])
    i_delta_2 = float(i_pm_2['delta_j2000'].iloc[0])

    mu_alpha = i_alpha_2 - i_alpha_1
    mu_delta = i_delta_2 - i_delta_1

    i_m_21 = mu_alpha ** 2 + (mu_delta ** 2) * cos(radians(i_delta_1))
    i_m_21 = sqrt(i_m_21)

    i_pm_21 = (3600 * i_m_21) / 1003
    i_pm_21 = i_pm_21 * 3600

    pm = i_pm['pm_values'].iloc[0]

    print pm, i_pm_21, (pm / i_pm_21)
"""

source_n = input('Catalog number: ')

o_file = fits.open('/mnt/e/Documentos/CarpetaCompartida/full_10_1.2_2.5_0.64_20-21_1.cat')
o_cat = Table(o_file[2].data).to_pandas()

o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_n])]
o_df = o_df[o_df['CATALOG_NUMBER'].isin([25])]
alpha_1 = float(o_df['ALPHA_J2000'].iloc[0])
delta_1 = float(o_df['DELTA_J2000'].iloc[0])
time_1 = float(o_df['EPOCH'].iloc[0])

o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_n])]
o_df = o_df[o_df['CATALOG_NUMBER'].isin([26])]
alpha_2 = float(o_df['ALPHA_J2000'].iloc[0])
delta_2 = float(o_df['DELTA_J2000'].iloc[0])
time_2 = float(o_df['EPOCH'].iloc[0])

o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_n])]
o_df = o_df[o_df['CATALOG_NUMBER'].isin([27])]
alpha_3 = float(o_df['ALPHA_J2000'].iloc[0])
delta_3 = float(o_df['DELTA_J2000'].iloc[0])
time_3 = float(o_df['EPOCH'].iloc[0])

# alpha_4 = float(o_df['ALPHA_J2000'].iloc[3])
# delta_4 = float(o_df['DELTA_J2000'].iloc[3])


print(time_1, time_2, time_3)

print(time_3 - time_2)
print(time_2 - time_1)

m_41 = hypot(abs(alpha_3 - alpha_1), abs(delta_3 - delta_1))
# m_43 = hypot(abs(alpha_4 - alpha_3), abs(delta_4 - delta_3))
m_32 = hypot(abs(alpha_3 - alpha_2), abs(delta_3 - delta_2))
m_21 = hypot(abs(alpha_2 - alpha_1), abs(delta_2 - delta_1))

# m_41 = hypot(alpha_4 - alpha_)

pm_21 = (3600*m_21)/1003
pm_21 = pm_21*3600

pm_32 = (3600*m_32)/1003
pm_32 = pm_32*3600

# pm_43 = (3600*m_43)/1003
# pm_43 = pm_43*3600

pm_41 = (3600*m_41)/2006
pm_41 = pm_41*3600

#
# f_pm = f_pm['PM'].iloc[0]
#
print pm_21, pm_32, pm_41
