from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

source = input('Type source number: ')
sources = [source]
# sources = [16067, 2091, 2013, 9637, 16902, 12179, 7340]
# sources = [13453, 4164, 4149, 21013, 13959, 2013, 19532,
#            18054, 20130, 3984, 12146, 4891, 11401]

fig = plt.figure(figsize=(11.7, 16.5 / 2), dpi=100)

with PdfPages('curves') as pdf:
    for source_ in sources:
        print('source {}'.format(source_))

        passed = []

        hdu_list = fits.open('/mnt/e/Documentos/CarpetaCompartida/full_10_1.2_2.5_0.64_20-21_1.cat')
        db = Table(hdu_list[2].data).to_pandas()

        ra = db.loc[db['SOURCE_NUMBER'] == source_, 'ALPHA_J2000'].tolist()
        dec = db.loc[db['SOURCE_NUMBER'] == source_, 'DELTA_J2000'].tolist()
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
            if fitted.rsquared <= float(r):
                passed.append(source_)
                params = fitted.params
                fitted = fitted.rsquared
                sigma_e = sigma
            """
            print "antes", dimension
            for dimension_ in dimension:
               dimension[dimension.index(dimension_)] = dimension_ * 3600
            print "despues", dimension
            plt.title('Source {} - {} - {}'.format(source_, coordinate, fitted.rsquared))
            plt.plot(epoch, dimension, 'bs')
            plt.xlabel('epoch')
            plt.ylabel('{}'.format(coordinate))
            plt.grid(True)
            # plt.show()
            # print("params", list(params))
            pdf.savefig()
            plt.clf()  # Clear current figure
            print('fitted {} - {}'.format(fitted.rsquared, coordinate))

        # for dec_ in dec:
            # epoch_ = epoch[dec.index(dec_)]
            # out = epoch_*params[1] + params[0]
            # print("input_dec", dec_)
            # print("output_dec", out)
            # print("sigma", sigma_e[dec.index(dec_)])
            # print("difference", out - dec_)

        # print(sigma_e)
        plt.title('Source {}'.format(source_))
        plt.plot(ra, dec, 'rs')
        plt.xlabel('ra')
        plt.ylabel('dec')
        plt.grid(True)
        # plt.show()
        pdf.savefig()
        plt.clf()  # Clear current figure
