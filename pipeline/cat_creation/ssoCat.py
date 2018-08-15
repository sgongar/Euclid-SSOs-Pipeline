# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 12:21:59 2016

@author: lconversi
"""
from sources import createObjectCatalogue
import numpy as np
import os

os.chdir('/pcdisk/holly/sgongora/Documents/Euclid/Euclid-tests/method_creation/')

catName = 'NoSSO'
outPath = '/pcdisk/holly/sgongora/Documents/Euclid/Euclid-tests/method_creation/cats/'

# Generate catalogue at 30 degrees ecliptic latitude to cover
# 1st quadrant and 4 dithers

# Set area to simulate: [4000,4000] pixels per CCD + 1000 in x and 3500 in y
# Set area to simulate: [4000,4000] pixels per CCD + 1000 in x and 3500 in y
# nxmax = 5000
# nymax = 8500
nxmax = 35000
nymax = 35000

# Adding stars + objects of type [8-15] (the same present in catalog30deg.dat)


settings = dict(besancon=False, deg=30, nx=nxmax, ny=nymax,
                outputprefix=outPath + catName, types=np.arange(8, 15),
                cutoff=25.5)
createObjectCatalogue.generateCatalog(**settings)


print "Catalogue created"

# Create SSO catalogue
nSso = 10
x0 = np.random.uniform(0, nxmax, nSso)
y0 = np.random.uniform(0, nymax, nSso)
mags = np.random.uniform(16, 26, nSso)
speeds = np.random.uniform(1, 10, nSso)
angles = np.random.uniform(0, 360, nSso)

# Convertion of speeds into indexes. Available speed files:
#   1-10 arcsec/h in steps of  1 arcsec/h as indexes 1001-1010
# 20-200 arcsec/h in steps of 10 arcsec/h as indexes 1020-1200

# Now rounding for indexes 1001-1010:
speeds[speeds < 10] = np.round(speeds[speeds < 10])

# Now rounding for indexes 1020-1200:
speeds[speeds >= 10] = np.round(speeds[speeds >= 10], -1)

# Save to catalogue
tbl = np.loadtxt(outPath + catName + '0.dat')

# print "before"
# print nSso + len(tbl[:, 0])

# print tbl[:, 0]

print "Create a new matrix using old objects and SSOs"
tblNew = np.zeros([nSso + len(tbl[:, 0]), 5])

tblNew[:, 0] = np.append(tbl[:, 0], x0)
tblNew[:, 1] = np.append(tbl[:, 1], y0)
tblNew[:, 2] = np.append(tbl[:, 2], mags)
tblNew[:, 3] = np.append(tbl[:, 3], speeds + 1000)
tblNew[:, 4] = np.append(tbl[:, 4], -angles)
np.savetxt(outPath + 'Cat_d1.dat', tblNew)


# Now computing SSOs position in the other dithers and save to catalogues

# 10 pixels/arcsec * delay between 2 exposure starts
# (1003 seconds from MOCD-B v6.8)

deltas = speeds / 3600. * 10. * 1003.


# Append to x and y position de new positions for the last 10 objects
# which are SSOs
tblNew[:, 0] = np.append(tbl[:, 0],
                         x0 + np.cos(angles / 180. * np.pi) * deltas)
tblNew[:, 1] = np.append(tbl[:, 1],
                         y0 + np.sin(angles / 180. * np.pi) * deltas)
np.savetxt(outPath + 'Cat_d2.dat', tblNew)

tblNew[:, 0] = np.append(tbl[:, 0],
                         x0 + np.cos(angles / 180. * np.pi) * deltas * 2)
tblNew[:, 1] = np.append(tbl[:, 1],
                         y0 + np.sin(angles / 180. * np.pi) * deltas * 2)
np.savetxt(outPath + 'Cat_d3.dat', tblNew)

tblNew[:, 0] = np.append(tbl[:, 0],
                         x0 + np.cos(angles / 180. * np.pi) * deltas * 3)
tblNew[:, 1] = np.append(tbl[:, 1],
                         y0 + np.sin(angles / 180. * np.pi) * deltas * 3)
np.savetxt(outPath + 'Cat_d4.dat', tblNew)

