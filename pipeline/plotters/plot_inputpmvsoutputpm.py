#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Plots pm and magnitude against A/B/elongation
Galaxies data

Versions:
- 0.1 Initial release

Todo:
    * Improve usability

"""
from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot


class PlotPMVersus:

    def __init__(self):
        self.false, self.true = self.create_data()
        self.plot_figure()

    def create_data(self):
        """

        :return:
        """
        false_ = [0.00708149410755, 0.00784190518844, 0.00945733734449,
                  0.0137489812909, 0.00557554709042, 0.0124727475937,
                  0.00571015712355, 0.00852627635795, 0.00527698084051,
                  0.00802728287399, 0.00922154678764, 0.0085320674053,
                  0.0146300205049, 0.00593097686412, 0.00606252435655,
                  0.0144730629764, 0.00615715950526, 0.00861821586604,
                  0.00763674405868, 0.0118552272692, 0.00589868619168,
                  0.00634269136197, 0.00609790782804, 0.00771253678077,
                  0.00543890437423, 0.00578728142084, 0.00535127521748,
                  0.00642188706193, 0.00586426364591, 0.00500576499253,
                  0.00671738943994, 0.00546860503107, 0.00561586990595,
                  0.0103717188655, 0.0052149187623, 0.00919017219935,
                  0.00527864118443, 0.00584207311081, 0.0069568577138,
                  0.00579263224269, 0.00683196039396, 0.00962968776938,
                  0.00630890786881, 0.00683049934008, 0.00820426542722,
                  0.00539424867254, 0.0087387545894, 0.00639375040589,
                  0.00854927786425, 0.00924485448394, 0.00551936734568,
                  0.00547425655286, 0.00523399246182, 0.00558087309167,
                  0.00618575187195, 0.00662178815341, 0.00624789536885,
                  0.00954451092041, 0.00780889215461, 0.00503281758277,
                  0.00689335745189, 0.00591808604065, 0.00597137642521,
                  0.0132085130872, 0.00707768634155, 0.00842168616788,
                  0.00759389510065, 0.0130973535043, 0.00912823259977]

        true_ = [0.00375248921686, 0.0134678393445, 0.0111898173707,
                 0.00930122232523, 0.00794271958452, 0.0109284472308,
                 0.0022504889285, 0.0100229905735, 0.00967448253243,
                 0.0143252564094, 0.0128308548971]

        return false_, true_

    def plot_figure(self):
        """

        :return:
        """
        with PdfPages('0.01".pdf') as pdf:
            fig = pyplot.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            ax.plot(self.false, [1] * len(self.false), 'bs')
            ax.plot(self.true, [2] * len(self.true), 'rs')

            pyplot.grid(True)

            pdf.savefig()


class PlotCHIVersus:

    def __init__(self):
        self.false, self.true = self.create_data()
        self.plot_figure()

    def get_chi(self, sources):
        """

        :param sources:
        :return:
        """
        merged_cat_loc = '/home/sgongora/Documents/CarpetaCompartida'
        merged_cat_file = 'merged_10_1.1_0.5_0.04_20-21_1.cat'
        merged_cat_loc = '{}/{}'.format(merged_cat_loc, merged_cat_file)
        hdu_list = fits.open(merged_cat_loc)
        hdu = Table(hdu_list[2].data).to_pandas()

        values_list = []
        for source_ in sources:
            chi_table = hdu[hdu['SOURCE_NUMBER'].isin([source_])]
            chi_value = chi_table['CHI2_ASTROM'].iloc[0]
            values_list.append(chi_value)

        return values_list

    def create_data(self):
        """

        :return:
        """
        false_sources = [34842, 30765, 30772, 30782, 41027, 30797, 10386,
                         22738, 20705, 43291, 18773, 2440, 2454, 10662, 10736,
                         16883, 14848, 14891, 15115, 10849, 23147, 4758, 45733,
                         30856, 27492, 37761, 35718, 37814, 19394, 43972, 25542,
                         19399, 35799, 15338, 37909, 27161, 46358, 25942, 38256,
                         5511, 46511, 3669, 46620, 11821, 11823, 9795, 36429,
                         9837, 30324, 15295, 20100, 38508, 11949, 38589, 20199,
                         3821, 36614, 30495, 44861, 14153, 40793, 26476, 20339,
                         30587, 20410, 30655, 47082, 12268,18419]
        false_chi = self.get_chi(false_sources)

        true_sources = [14444, 20753, 24969, 45495, 14785, 14791,
                        20994, 23176, 8979, 42189, 3889]
        true_chi = self.get_chi(true_sources)

        return false_chi, true_chi

    def plot_figure(self):
        """

        :return:
        """
        with PdfPages('chi_0.01".pdf') as pdf:
            fig = pyplot.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            ax.plot(self.false, [1] * len(self.false), 'bs')
            ax.plot(self.true, [2] * len(self.true), 'rs')

            pyplot.grid(True)

            pdf.savefig()


if __name__ == "__main__":
    PlotCHIVersus()
