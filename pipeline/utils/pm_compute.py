#!/usr/bin/python
# -*- coding: utf-8 -*-

# From ra/dec coordinates
"""
from math import hypot

i_alpha_1 = 99.6819699747
i_alpha_2 = 99.681967793
i_alpha_3 = 99.6819656113

i_delta_1 = 0.0311318211236
i_delta_2 = 0.0311344802527
i_delta_3 = 0.0311371393818

i_m_32 = hypot(i_alpha_3 - i_alpha_2, 
               i_delta_3 - i_delta_2)
i_m_21 = hypot(i_alpha_2 - i_alpha_1,
               i_delta_2 - i_delta_1)

i_pm_32 = (3600 * i_m_32) / 1003
i_pm_32 = i_pm_32 * 3600

print "in", i_pm_32

o_alpha_1 = 99.68200883
o_alpha_2 = 99.6819902
o_alpha_3 = 99.68198154

o_delta_1 = 0.03110524517
o_delta_2 = 0.03111189309
o_delta_3 = 0.03112198935

o_m_32 = hypot(o_alpha_3 - o_alpha_2,
               o_delta_3 - o_delta_2)
o_m_21 = hypot(o_alpha_2 - o_alpha_1,
               o_delta_2 - o_delta_1)

o_pm_32 = (3600 * o_m_32) / 1003
o_pm_32 = o_pm_32 * 3600

print "out", o_pm_32
"""
# From x/y coordinates
from math import hypot, sqrt
from math import cos, radians

from astropy.io import fits
from astropy.wcs import WCS

"""
# first object 0.0003
i_x_1 = 3639.09734476982
i_x_2 = 3639.0975329218
i_x_3 = 3639.09772107378

i_y_1 = 9306.29708939264
i_y_2 = 9306.29799883719
i_y_3 = 9306.29890828175

# second object 0.001
i_x_1 = 11445.6597063931
i_x_2 = 11445.6615949739
i_x_3 = 11445.6634835547

i_y_1 = 10906.1932248926
i_y_2 = 10906.195273228
i_y_3 = 10906.1973215633

# third object 0.003
i_x_1 = 7912.01135466143
i_x_2 = 7912.02063305155
i_x_3 = 7912.02991144166

i_y_1 = 1205.23029694468
i_y_2 = 1205.23069761134
i_y_3 = 1205.231098278

# forth object 0.01
i_x_1 = 6360.44023580468
i_x_2 = 6360.46686009872
i_x_3 = 6360.49348439276

i_y_1 = 7191.09837273731
i_y_2 = 7191.10658178593
i_y_3 = 7191.11479083455

# fifth object 0.03
i_x_1 = 13285.2894546902
i_x_2 = 13285.3551935217
i_x_3 = 13285.4209323531

i_y_1 = 3702.61412812138
i_y_2 = 3702.67972775349
i_y_3 = 3702.74532738559

# sixth object 0.1
i_x_1 = 5920.70195442806
i_x_2 = 5920.92451048715
i_x_3 = 5921.14706654624

i_y_1 = 11640.2784258782
i_y_2 = 11640.4460354
i_y_3 = 11640.6136449218

# seventh object 0.3
i_x_1 = 4986.99569245097
i_x_2 = 4987.07025047706
i_x_3 = 4987.14480850314

i_y_1 = 2929.31893847003
i_y_2 = 2930.24464450869
i_y_3 = 2931.17035054736

# eighth object 1
i_x_1 = 5280.94136043445
i_x_2 = 5287.53489685115
i_x_3 = 5294.12843326784

i_y_1 = 9383.98181934815
i_y_2 = 9389.11864849003
i_y_3 = 9394.25547763191
"""
fits_dir = '/home/sgongora/Documents/CarpetaCompartida/luca_data/v14/CCDs/'
ccd_loc = 'mag_20-21_CCD_x0_y0_d1.fits'
fits_loc = '{}/{}'.format(fits_dir, ccd_loc)
hdulist = fits.open(fits_loc)
w = WCS(hdulist[0].header)

regions_list = [[i_x_1, i_y_1], [i_x_2, i_y_2], [i_x_3, i_y_3]]
input_regions = w.all_pix2world(regions_list, 1)

# Coordinates
i_alpha_1 = float(input_regions[0][0])
i_alpha_2 = float(input_regions[1][0])
i_alpha_3 = float(input_regions[2][0])

i_delta_1 = float(input_regions[0][1])
i_delta_2 = float(input_regions[1][1])
i_delta_3 = float(input_regions[2][1])

# First method
i_m1_32 = hypot(i_alpha_3 - i_alpha_2,
                i_delta_3 - i_delta_2)
i_m1_21 = hypot(i_alpha_2 - i_alpha_1,
                i_delta_2 - i_delta_1)

i_pm1_21 = (3600 * i_m1_21) / 1003
i_pm1_21 = i_pm1_21 * 3600
print("%.2f" % i_pm1_21)

# Second method
mu_alpha = i_alpha_2 - i_alpha_1
mu_delta = i_delta_2 - i_delta_1

i_m2_21 = mu_alpha**2 + (mu_delta**2)*cos(radians(i_delta_1))
i_m2_21 = sqrt(i_m2_21)
i_pm2_21 = (3600 * i_m2_21) / 1003
i_pm2_21 = i_pm2_21 * 3600
print("%.2f" % i_pm2_21)

# Differences
# print "difference", i_pm1_21 - i_pm2_21
