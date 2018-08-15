#!/usr/bin/python
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
from numpy import sqrt

source = input('Type source number: ')

hdu_list = fits.open('/mnt/e/Documentos/CarpetaCompartida/merged_10_1.2_2.5_0.64_20-21_1.cat')
db = Table(hdu_list[2].data).to_pandas()

pm_alpha = db.loc[db['SOURCE_NUMBER'] == source, 'PMALPHA_J2000'].tolist()
pm_delta = db.loc[db['SOURCE_NUMBER'] == source, 'PMDELTA_J2000'].tolist()
pme_alpha = db.loc[db['SOURCE_NUMBER'] == source, 'PMALPHAERR_J2000'].tolist()
pme_delta = db.loc[db['SOURCE_NUMBER'] == source, 'PMDELTAERR_J2000'].tolist()

pm_alpha = float(pm_alpha[0]) / 8.75e6
pm_delta = float(pm_delta[0]) / 8.75e6
pme_alpha = float(pme_alpha[0]) / 8.75e6
pme_delta = float(pme_delta[0]) / 8.75e6

pm = sqrt(pm_alpha**2 + pm_delta**2)
pme = sqrt(pme_alpha**2 + pme_delta**2)

print('pm', pm)
print('pme', pme)
print('SN', pm / pme)

# pmalpha = Series(merged_db.field('PMALPHA_J2000') / 8.75e6)  # 8.75e6

"""
epoch = db.loc[db['SOURCE_NUMBER'] == source_, 'EPOCH'].tolist()

r = 0.95

edims = []
sigmas = []
epochs = []
dimensions = []

for dimension in [ra, dec]:
    x = np.array(epoch)
    epochs.append(x)  # epochs list
    y = np.array(dimension)
    dimensions.append(y)  # dimensions list
    if dimension == ra:
        sigma = db.loc[db['SOURCE_NUMBER'] == source_,
                       'ERRA_WORLD'].tolist()
        coordinate = 'ra'
    if dimension == dec:
        sigma = db.loc[db['SOURCE_NUMBER'] == source_,
                       'ERRB_WORLD'].tolist()
        coordinate = 'dec'
    edim = np.array([1 / var for var in sigma])
    edims.append(edim)
    # print("edim", edim)
    x = sm.add_constant(x)
    model = sm.WLS(y, x, weigths=edim)
    fitted = model.fit()
    # print fitted.summary()
"""