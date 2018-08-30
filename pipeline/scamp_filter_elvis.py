#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
For now intermediate versions are stored in physical memory, this will not
be longer needed in the future but in this moment could be useful for
pipeline development reasons.


Versions:
- 0.1: Initial release.

Todo:
    * Improve documentation

*GNU Terry Pratchett*
"""
from multiprocessing import Process

from astropy.io import fits
from astropy.table import Table
from numpy import mean, median, poly1d
from pandas import concat, read_csv, Series

from misc import pm_compute, extract_settings_elvis
from misc import check_source_elvis, confidence_filter

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class ScampFilterELViS:  # TODO Split scamp_filter method into single methods

    def __init__(self, logger, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_d:
        """
        # Analysis variables
        self.prfs_d = extract_settings_elvis()
        self.logger = logger

        self.save = True

        # Filtered catalog dir
        self.filter_dir = self.prfs_d['filtered']
        self.filt_n = 'filt_'
        self.filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], self.filt_n)

        # Saves _1.csv
        (merged_db, full_db) = self.scamp_filter()
        # Saves _2.csv
        full_df = self.compute_pm(merged_db, full_db)
        # Saves _3.csv
        full_df = self.get_areas(full_df)
        # Saves _4.csv
        full_df = full_df[full_df['PM'] > 0.01]
        # Saves _5.csv
        full_df = self.filter_class(full_df)

        fast_df = full_df[full_df['PM'] > 2]
        slow_df = full_df[full_df['PM'] < 2]

        fast_df = self.filter_coherence(fast_df)
        slow_df = self.filter_b_image(slow_df)  # 8th version

        full_df = concat([fast_df, slow_df])

        if self.save:
            self.save_message('9')
            full_df.to_csv('{}_9.csv'.format(self.filter_o_n))

    def save_message(self, order):
        """

        :param order:
        :return:
        """
        self.logger.debug('Saves data to: ')
        self.logger.debug('Dir: {}'.format(self.filter_dir))
        self.logger.debug('Name: {}_{}.csv'.format(self.filt_n, order))

    def get_cat(self, cat_n):
        """

        :param cat_n:
        :return: cat_file
        """
        cats = ['empty_cat']

        for x_ in range(1, 7, 1):
            for y_ in range(1, 7, 1):
                for d_ in range(1, 5, 1):
                    cat_name = 'CCD_x{}_y{}_d{}.cat'.format(x_, y_, d_)
                    cats.append(cat_name)

        return cats[cat_n]

    def scamp_filter(self):
        """

        :return:
        """
        self.logger.info("Filtering scamp's output")

        # Getting full catalog
        # Full catalog name
        full_cat_n = 'full_1.cat'  # todo - hardoced
        full_n = '{}/{}'.format(self.prfs_d['output_cats'], full_cat_n)

        # Just for logging reasons
        self.logger.debug('Opens full catalog')
        self.logger.debug('Dir: {}'.format(self.prfs_d['output_cats']))
        self.logger.debug('Name: {}'.format(full_cat_n))
        full_cat = fits.open(full_n)  # Opens full catalog
        full_db = Table(full_cat[2].data)  # Converts it to Astropy Table
        self.logger.debug('Converts full Astropy catalog to Pandas format')
        full_db = full_db.to_pandas()  # Converts it to Pandas format

        # Getting merge catalog
        # Merged catalog name
        mrgd_cat_n = 'merged_1.cat'
        mrgd_n = '{}/{}'.format(self.prfs_d['output_cats'], mrgd_cat_n)

        # Just for logging reasons
        self.logger.debug('Opens merged catalog')
        self.logger.debug('Dir: {}'.format(self.prfs_d['output_cats']))
        self.logger.debug('Name: merged_1.cat')  # Be careful, hardcoded!
        merged_cat = fits.open(mrgd_n)  # Opens merged catalog
        self.logger.debug('Converts merged Astropy catalog to Pandas format')
        # todo - checks if Pandas format gonna be useful!
        merged_db = Table(merged_cat[2].data)

        merged_db = merged_db.to_pandas()

        # Removing 0 catalog detections
        self.logger.debug('Removes 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        # Filter by 3 or more detections
        full_db, merged_db = self.filter_detections(full_db, merged_db, 3)

        if self.save:
            self.save_message('1')
            full_db.to_csv('{}_full_1.csv'.format(self.filter_o_n))
            merged_db.to_csv('{}_merged_1.csv'.format(self.filter_o_n))

        return merged_db, full_db

    def compute_pm(self, merged_db, full_db):
        """

        :param merged_db:
        :param full_db:
        :return: full_db
        """
        self.logger.debug('Computes proper motion')
        full_db = pm_compute(self.logger, merged_db, full_db)
        if self.save:
            self.save_message('2')
            full_db.to_csv('{}_2.csv'.format(self.filter_o_n))

        return full_db

    def get_areas(self, full_df):
        """

        :return: full_df
        """
        self.logger.debug('Populates filtered catalog with Sextractor data')

        # Gets unique sources from filtered file
        unique_sources = list(set(full_df['SOURCE_NUMBER'].tolist()))
        l_sourcs = len(unique_sources)  # Just to not break 79 characters
        self.logger.debug('Unique sources to be analysed {}'.format(l_sourcs))

        dict_keys = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                     'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE',
                     'ISOAREA_IMAGE', 'A_IMAGE', 'MEDIAN_A_IMAGE',
                     'MEAN_A_IMAGE', 'ERRA_IMAGE',  'MEDIAN_ERRA_IMAGE',
                     'MEAN_ERRA_IMAGE', 'B_IMAGE', 'MEDIAN_B_IMAGE',
                     'MEAN_B_IMAGE', 'ERRB_IMAGE', 'MEDIAN_ERRB_IMAGE',
                     'MEAN_ERRB_IMAGE', 'THETA_IMAGE', 'ERRTHETA_IMAGE',
                     'ALPHA_J2000', 'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD',
                     'ERRTHETA_WORLD', 'EPOCH', 'FWHM_IMAGE', 'CLASS_STAR',
                     'MEDIAN_CLASS_STAR', 'MEAN_CLASS_STAR', 'FLUX_ISO',
                     'MEDIAN_FLUX_ISO', 'MEAN_FLUX_ISO', 'FLUXERR_ISO',
                     'MEDIAN_FLUXERR_ISO', 'MEAN_FLUXERR_ISO', 'FLUX_RADIUS',
                     'ELONGATION', 'ELLIPTICITY', 'MEDIAN_ELLIPTICITY',
                     'MEAN_ELLIPTICITY', 'MAG', 'MAGERR', 'MAG_ISO',
                     'MEDIAN_MAG_ISO', 'MEAN_MAG_ISO', 'MAGERR_ISO',
                     'MEDIAN_MAGERR_ISO', 'MEAN_MAGERR_ISO',
                     'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'PM',
                     'PMERR', 'PMALPHA', 'PMDELTA', 'PMALPHAERR', 'PMDELTAERR',
                     'MAG_AUTO', 'MEDIAN_MAG_AUTO', 'MEAN_MAG_AUTO',
                     'MAGERR_AUTO', 'MEDIAN_MAGERR_AUTO', 'MEAN_MAGERR_AUTO']

        sub_list_size = len(unique_sources) / self.prfs_d['cores_number']

        sub_list_l = []
        for idx_sub_list in range(0, self.prfs_d['cores_number'], 1):
            if idx_sub_list != (self.prfs_d['cores_number'] - 1):
                idx_down = sub_list_size * idx_sub_list
                idx_up = sub_list_size * (idx_sub_list + 1)
                sub_list_l.append(unique_sources[idx_down:idx_up])
            else:
                idx_down = sub_list_size * idx_sub_list
                sub_list_l.append(unique_sources[idx_down:])

        cats_number = 144
        cat_d = {}
        for cat_n in range(1, cats_number + 1, 1):
            cat_file = self.get_cat(cat_n)
            cat_data = fits.open('{}/{}'.format(self.prfs_d['fits_dir'],
                                                cat_file))

            ccd_df = Table(cat_data[2].data)
            # self.logger.debug('CCD catalog {} to Pandas'.format(cat_n))
            cat_d[cat_n] = ccd_df.to_pandas()

        areas_j = []
        for idx_l in range(0, self.prfs_d['cores_number'], 1):
            areas_p = Process(target=self.get_areas_thread,
                              args=(dict_keys, sub_list_l[idx_l],
                                    full_df, idx_l, cat_d,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges areas
        # Merges catalogs
        csv_list = []
        for idx_csv in range(0, self.prfs_d['cores_number'], 1):
            csv_ = read_csv('{}_3_{}.csv'.format(self.filter_o_n, idx_csv),
                            index_col=0)
            csv_list.append(csv_)

        full_df = concat(csv_list)

        if self.save:
            self.save_message('3')
            full_df.to_csv('{}_3.csv'.format(self.filter_o_n))

        return full_df

    def get_areas_thread(self, dict_keys, unique_sources_thread,
                         filter_cat, idx_l, cat_d):
        """

        :param dict_keys:
        :param unique_sources_thread:
        :param filter_cat:
        :param idx_l:
        :param cat_d:
        :return:
        """
        tmp_d_keys = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                      'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE',
                      'ERRA_IMAGE', 'ERRB_IMAGE', 'ERRTHETA_IMAGE',
                      'ALPHA_J2000', 'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD',
                      'ERRTHETA_WORLD', 'EPOCH', 'MAG', 'MAGERR',
                      'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'PM',
                      'PMERR', 'PMALPHA', 'PMDELTA', 'PMALPHAERR',
                      'PMDELTAERR', 'THETA_IMAGE', 'ISOAREA_IMAGE',
                      'FWHM_IMAGE', 'ELONGATION', 'CLASS_STAR', 'A_IMAGE',
                      'B_IMAGE', 'ERRA_IMAGE', 'ERRB_IMAGE', 'CLASS_STAR',
                      'FLUX_ISO', 'FLUXERR_ISO', 'FLUX_RADIUS', 'MAG_ISO',
                      'MAGERR_ISO', 'MAG_AUTO', 'MAGERR_AUTO', 'ELLIPTICITY',
                      'MEDIAN_A_IMAGE', 'MEDIAN_B_IMAGE', 'MEDIAN_ERRA_IMAGE',
                      'MEDIAN_ERRB_IMAGE', 'MEDIAN_CLASS_STAR',
                      'MEDIAN_FLUX_ISO', 'MEDIAN_FLUXERR_ISO',
                      'MEDIAN_MAG_ISO', 'MEDIAN_MAGERR_ISO',
                      'MEDIAN_MAG_AUTO', 'MEDIAN_MAGERR_AUTO',
                      'MEDIAN_ELLIPTICITY', 'MEAN_A_IMAGE', 'MEAN_B_IMAGE',
                      'MEAN_ERRA_IMAGE', 'MEAN_ERRB_IMAGE', 'MEAN_CLASS_STAR',
                      'MEAN_FLUX_ISO', 'MEAN_FLUXERR_ISO', 'MEAN_MAG_ISO',
                      'MEAN_MAGERR_ISO', 'MEAN_MAG_AUTO', 'MEAN_MAGERR_AUTO',
                      'MEAN_ELLIPTICITY']
        tmp_d = {}
        for key_ in tmp_d_keys:
            tmp_d[key_] = []

        cat_n_l = []
        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources_thread):
            source_d_keys = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                             'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE',
                             'Y_IMAGE', 'ERRA_IMAGE', 'ERRB_IMAGE',
                             'ERRTHETA_IMAGE', 'ALPHA_J2000', 'DELTA_J2000',
                             'ERRA_WORLD', 'ERRB_WORLD', 'ERRTHETA_WORLD',
                             'EPOCH', 'MAG', 'MAGERR', 'FLAGS_EXTRACTION',
                             'FLAGS_SCAMP', 'FLAGS_IMA', 'PM', 'PMERR',
                             'PMALPHA', 'PMDELTA', 'PMALPHAERR', 'PMDELTAERR',
                             'THETA_IMAGE', 'ISOAREA_IMAGE', 'FWHM_IMAGE',
                             'ELONGATION', 'CLASS_STAR', 'A_IMAGE', 'B_IMAGE',
                             'ERRA_IMAGE', 'ERRB_IMAGE', 'CLASS_STAR',
                             'FLUX_ISO', 'FLUXERR_ISO', 'FLUX_RADIUS',
                             'MAG_ISO', 'MAGERR_ISO', 'MAG_AUTO', 'MAGERR_AUTO',
                             'ELLIPTICITY']
            source_d = {}
            for key_ in source_d_keys:
                source_d[key_] = []

            o_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]
            length_source = o_df['SOURCE_NUMBER'].size
            right_source = 0
            for i, row in enumerate(o_df.itertuples(), 1):
                # Populates temporal dictionary
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                cat_n = row.CATALOG_NUMBER
                cat_n_l.append(cat_n)
                # self.logger.debug('opening CCD catalog {}'.format(cat_file))

                cat_df = check_source_elvis(cat_d[cat_n], o_alpha, o_delta)
                if cat_df.empty:
                    # FIXME problemas entre sextractor y scamp
                    # Como el objeto no esta en las imagenes originales de
                    # sextractor lo borro del catalogo de scamp
                    # print('cat_empty')
                    # for idx_k, key_ in enumerate(keys_l):
                    #     tmp_d[key_].append('nan')
                    # for idx_k, key_ in enumerate(extra_keys):
                    #     tmp_d[key_].append('nan')
                    pass
                else:
                    right_source += 1
                    # print('cat not empty')
                    source_d['SOURCE_NUMBER'].append(row.SOURCE_NUMBER)
                    source_d['CATALOG_NUMBER'].append(row.CATALOG_NUMBER)
                    source_d['EXTENSION'].append(row.EXTENSION)
                    source_d['ASTR_INSTRUM'].append(row.ASTR_INSTRUM)
                    source_d['PHOT_INSTRUM'].append(row.PHOT_INSTRUM)
                    source_d['X_IMAGE'].append(row.X_IMAGE)
                    source_d['Y_IMAGE'].append(row.Y_IMAGE)
                    source_d['ERRA_IMAGE'].append(row.ERRA_IMAGE)
                    source_d['ERRB_IMAGE'].append(row.ERRB_IMAGE)
                    source_d['ERRTHETA_IMAGE'].append(row.ERRTHETA_IMAGE)
                    source_d['ALPHA_J2000'].append(row.ALPHA_J2000)
                    source_d['DELTA_J2000'].append(row.DELTA_J2000)
                    source_d['ERRA_WORLD'].append(row.ERRA_WORLD)
                    source_d['ERRB_WORLD'].append(row.ERRB_WORLD)
                    source_d['ERRTHETA_WORLD'].append(row.ERRTHETA_WORLD)
                    source_d['EPOCH'].append(row.EPOCH)
                    source_d['MAG'].append(row.MAG)
                    source_d['MAGERR'].append(row.MAGERR)
                    source_d['FLAGS_EXTRACTION'].append(row.FLAGS_EXTRACTION)
                    source_d['FLAGS_SCAMP'].append(row.FLAGS_SCAMP)
                    source_d['FLAGS_IMA'].append(row.FLAGS_IMA)
                    source_d['PM'].append(row.PM)
                    source_d['PMERR'].append(row.PMERR)
                    source_d['PMALPHA'].append(row.PMALPHA)
                    source_d['PMDELTA'].append(row.PMDELTA)
                    source_d['PMALPHAERR'].append(row.PMALPHAERR)
                    source_d['PMDELTAERR'].append(row.PMDELTAERR)
                    source_d['THETA_IMAGE'].append(cat_df['THETA_IMAGE'].iloc[0])
                    source_d['ISOAREA_IMAGE'].append(cat_df['ISOAREA_IMAGE'].iloc[0])
                    source_d['FWHM_IMAGE'].append(cat_df['FWHM_IMAGE'].iloc[0])
                    source_d['ELONGATION'].append(cat_df['ELONGATION'].iloc[0])
                    source_d['CLASS_STAR'].append(cat_df['CLASS_STAR'].iloc[0])
                    source_d['A_IMAGE'].append(cat_df['A_IMAGE'].iloc[0])
                    source_d['B_IMAGE'].append(cat_df['B_IMAGE'].iloc[0])
                    source_d['ERRA_IMAGE'].append(cat_df['ERRA_IMAGE'].iloc[0])
                    source_d['ERRB_IMAGE'].append(cat_df['ERRB_IMAGE'].iloc[0])
                    source_d['CLASS_STAR'].append(cat_df['CLASS_STAR'].iloc[0])
                    source_d['FLUX_ISO'].append(cat_df['FLUX_ISO'].iloc[0])
                    source_d['FLUXERR_ISO'].append(cat_df['FLUXERR_ISO'].iloc[0])
                    source_d['FLUX_RADIUS'].append(cat_df['FLUX_RADIUS'].iloc[0])
                    source_d['MAG_ISO'].append(cat_df['MAG_ISO'].iloc[0])
                    source_d['MAGERR_ISO'].append(cat_df['MAGERR_ISO'].iloc[0])
                    source_d['MAG_AUTO'].append(cat_df['MAG_AUTO'].iloc[0])
                    source_d['MAGERR_AUTO'].append(cat_df['MAGERR_AUTO'].iloc[0])
                    source_d['ELLIPTICITY'].append(cat_df['ELLIPTICITY'].iloc[0])

            if length_source == right_source:
                mean_a_image = mean(source_d['A_IMAGE'])
                median_a_image = median(source_d['A_IMAGE'])

                mean_b_image = mean(source_d['B_IMAGE'])
                median_b_image = median(source_d['B_IMAGE'])

                mean_erra_image = mean(source_d['ERRA_IMAGE'])
                median_erra_image = median(source_d['ERRA_IMAGE'])

                mean_errb_image = mean(source_d['ERRB_IMAGE'])
                median_errb_image = median(source_d['ERRB_IMAGE'])

                mean_class_star = mean(source_d['CLASS_STAR'])
                median_class_star = median(source_d['CLASS_STAR'])

                mean_flux_iso = mean(source_d['FLUX_ISO'])
                median_flux_iso = median(source_d['FLUX_ISO'])

                mean_fluxerr_iso = mean(source_d['FLUXERR_ISO'])
                median_fluxerr_iso = median(source_d['FLUXERR_ISO'])

                mean_mag_iso = mean(source_d['MAG_ISO'])
                median_mag_iso = median(source_d['MAG_ISO'])

                mean_magerr_iso = mean(source_d['MAGERR_ISO'])
                median_magerr_iso = median(source_d['MAGERR_ISO'])

                mean_mag_auto = mean(source_d['MAG_AUTO'])
                median_mag_auto = median(source_d['MAG_AUTO'])

                mean_magerr_auto = mean(source_d['MAGERR_AUTO'])
                median_magerr_auto = median(source_d['MAGERR_AUTO'])

                mean_ellipticity = mean(source_d['ELLIPTICITY'])
                median_ellipticiy = median(source_d['ELLIPTICITY'])

                # Saves mean and median values from sources
                for i_stats in range(0, len(o_df['SOURCE_NUMBER']), 1):
                    tmp_d['SOURCE_NUMBER'].append(source_d['SOURCE_NUMBER'][i_stats])
                    tmp_d['CATALOG_NUMBER'].append(source_d['CATALOG_NUMBER'][i_stats])
                    tmp_d['EXTENSION'].append(source_d['EXTENSION'][i_stats])
                    tmp_d['ASTR_INSTRUM'].append(source_d['ASTR_INSTRUM'][i_stats])
                    tmp_d['PHOT_INSTRUM'].append(source_d['PHOT_INSTRUM'][i_stats])
                    tmp_d['X_IMAGE'].append(source_d['X_IMAGE'][i_stats])
                    tmp_d['Y_IMAGE'].append(source_d['Y_IMAGE'][i_stats])
                    tmp_d['ERRA_IMAGE'].append(source_d['ERRA_IMAGE'][i_stats])
                    tmp_d['ERRB_IMAGE'].append(source_d['ERRB_IMAGE'][i_stats])
                    tmp_d['ERRTHETA_IMAGE'].append(source_d['ERRTHETA_IMAGE'][i_stats])
                    tmp_d['ALPHA_J2000'].append(source_d['ALPHA_J2000'][i_stats])
                    tmp_d['DELTA_J2000'].append(source_d['DELTA_J2000'][i_stats])
                    tmp_d['ERRA_WORLD'].append(source_d['ERRA_WORLD'][i_stats])
                    tmp_d['ERRB_WORLD'].append(source_d['ERRB_WORLD'][i_stats])
                    tmp_d['ERRTHETA_WORLD'].append(source_d['ERRTHETA_WORLD'][i_stats])
                    tmp_d['EPOCH'].append(source_d['EPOCH'][i_stats])
                    tmp_d['MAG'].append(source_d['MAG'][i_stats])
                    tmp_d['MAGERR'].append(source_d['MAGERR'][i_stats])
                    tmp_d['FLAGS_EXTRACTION'].append(source_d['FLAGS_EXTRACTION'][i_stats])
                    tmp_d['FLAGS_SCAMP'].append(source_d['FLAGS_SCAMP'][i_stats])
                    tmp_d['FLAGS_IMA'].append(source_d['FLAGS_IMA'][i_stats])
                    tmp_d['PM'].append(source_d['PM'][i_stats])
                    tmp_d['PMERR'].append(source_d['PMERR'][i_stats])
                    tmp_d['PMALPHA'].append(source_d['PMALPHA'][i_stats])
                    tmp_d['PMDELTA'].append(source_d['PMDELTA'][i_stats])
                    tmp_d['PMALPHAERR'].append(source_d['PMALPHAERR'][i_stats])
                    tmp_d['PMDELTAERR'].append(source_d['PMDELTAERR'][i_stats])
                    tmp_d['A_IMAGE'].append(source_d['A_IMAGE'][i_stats])
                    tmp_d['B_IMAGE'].append(source_d['B_IMAGE'][i_stats])
                    tmp_d['THETA_IMAGE'].append(source_d['THETA_IMAGE'][i_stats])
                    tmp_d['ISOAREA_IMAGE'].append(source_d['ISOAREA_IMAGE'][i_stats])
                    tmp_d['FWHM_IMAGE'].append(source_d['FWHM_IMAGE'][i_stats])
                    tmp_d['FLUX_ISO'].append(source_d['FLUX_ISO'][i_stats])
                    tmp_d['FLUXERR_ISO'].append(source_d['FLUXERR_ISO'][i_stats])
                    tmp_d['FLUX_RADIUS'].append(source_d['FLUX_RADIUS'][i_stats])
                    tmp_d['MAG_ISO'].append(source_d['MAG_ISO'][i_stats])
                    tmp_d['MAGERR_ISO'].append(source_d['MAGERR_ISO'][i_stats])
                    tmp_d['ELONGATION'].append(source_d['ELONGATION'][i_stats])
                    tmp_d['ELLIPTICITY'].append(source_d['ELLIPTICITY'][i_stats])
                    tmp_d['CLASS_STAR'].append(source_d['CLASS_STAR'][i_stats])
                    tmp_d['MEDIAN_A_IMAGE'].append(median_a_image)
                    tmp_d['MEDIAN_B_IMAGE'].append(median_b_image)
                    tmp_d['MEDIAN_ERRA_IMAGE'].append(median_erra_image)
                    tmp_d['MEDIAN_ERRB_IMAGE'].append(median_errb_image)
                    tmp_d['MEDIAN_CLASS_STAR'].append(median_class_star)
                    tmp_d['MEDIAN_FLUX_ISO'].append(median_flux_iso)
                    tmp_d['MEDIAN_FLUXERR_ISO'].append(median_fluxerr_iso)
                    tmp_d['MEDIAN_MAG_ISO'].append(median_mag_iso)
                    tmp_d['MEDIAN_MAGERR_ISO'].append(median_magerr_iso)
                    tmp_d['MEDIAN_MAG_AUTO'].append(median_mag_auto)
                    tmp_d['MEDIAN_MAGERR_AUTO'].append(median_magerr_auto)
                    tmp_d['MEDIAN_ELLIPTICITY'].append(median_ellipticiy)
                    tmp_d['MEAN_A_IMAGE'].append(mean_a_image)
                    tmp_d['MEAN_B_IMAGE'].append(mean_b_image)
                    tmp_d['MEAN_ERRA_IMAGE'].append(mean_erra_image)
                    tmp_d['MEAN_ERRB_IMAGE'].append(mean_errb_image)
                    tmp_d['MEAN_CLASS_STAR'].append(mean_class_star)
                    tmp_d['MEAN_FLUX_ISO'].append(mean_flux_iso)
                    tmp_d['MEAN_FLUXERR_ISO'].append(mean_fluxerr_iso)
                    tmp_d['MEAN_MAG_ISO'].append(mean_mag_iso)
                    tmp_d['MEAN_MAGERR_ISO'].append(mean_magerr_iso)
                    tmp_d['MEAN_MAG_AUTO'].append(mean_mag_auto)
                    tmp_d['MEAN_MAGERR_AUTO'].append(mean_magerr_auto)
                    tmp_d['MEAN_ELLIPTICITY'].append(mean_ellipticity)
            else:
                # Wrong source
                pass

        series_l = []
        series_d = {}
        for key_ in tmp_d.keys():
            series_d[key_] = Series(tmp_d[key_], name=key_)
            # print('idx {} - size {}'.format(idx_l, len(series_d[key_])))
            series_l.append(series_d[key_])

        full_db = concat(series_l, axis=1)

        if self.save:
            self.save_message('3_{}'.format(idx_l))
            full_db.to_csv('{}_3_{}.csv'.format(self.filter_o_n, idx_l),
                           columns=dict_keys)

    def filter_class(self, full_df):
        """

        :param full_df:
        :return:
        """
        rejected = []
        accepted = []
        unique_sources = list(set(full_df['SOURCE_NUMBER'].tolist()))

        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources):
            o_df = full_df[full_df['SOURCE_NUMBER'].isin([source_])].iloc[0]
            # mag = float(o_df['MEDIAN_MAG_ISO'])

            # b test
            pm = float(o_df['PM'])
            class_star = float(o_df['MEAN_CLASS_STAR'])

            if pm < 1.0 and class_star < 0.65:
                rejected.append(source_)
            else:
                accepted.append(source_)

        full_df = full_df[full_df['SOURCE_NUMBER'].isin(accepted)]

        if self.save:
            self.save_message('4')
            full_df.to_csv('{}_4.csv'.format(self.filter_o_n))

        return full_df

    def filter_b_image(self, full_df):
        """

        :return: full_df
        """
        self.logger.debug('Runs B_Image size filter')

        # Gets unique sources from filtered file
        unique_sources = list(set(full_df['SOURCE_NUMBER'].tolist()))
        l_sourcs = len(unique_sources)  # Just to not break 79 characters
        self.logger.debug('Unique sources to be analysed {}'.format(l_sourcs))

        dict_keys = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                     'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE',
                     'ISOAREA_IMAGE', 'A_IMAGE', 'MEDIAN_A_IMAGE',
                     'MEAN_A_IMAGE', 'ERRA_IMAGE',  'MEDIAN_ERRA_IMAGE',
                     'MEAN_ERRA_IMAGE', 'B_IMAGE', 'MEDIAN_B_IMAGE',
                     'MEAN_B_IMAGE', 'ERRB_IMAGE', 'MEDIAN_ERRB_IMAGE',
                     'MEAN_ERRB_IMAGE', 'THETA_IMAGE', 'ERRTHETA_IMAGE',
                     'ALPHA_J2000', 'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD',
                     'ERRTHETA_WORLD', 'EPOCH', 'FWHM_IMAGE', 'CLASS_STAR',
                     'MEDIAN_CLASS_STAR', 'MEAN_CLASS_STAR', 'FLUX_ISO',
                     'MEDIAN_FLUX_ISO', 'MEAN_FLUX_ISO', 'FLUXERR_ISO',
                     'MEDIAN_FLUXERR_ISO', 'MEAN_FLUXERR_ISO', 'FLUX_RADIUS',
                     'ELONGATION', 'ELLIPTICITY', 'MEDIAN_ELLIPTICITY',
                     'MEAN_ELLIPTICITY', 'MAG', 'MAGERR', 'MAG_ISO',
                     'MEDIAN_MAG_ISO', 'MEAN_MAG_ISO', 'MAGERR_ISO',
                     'MEDIAN_MAGERR_ISO', 'MEAN_MAGERR_ISO',
                     'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'PM',
                     'PMERR', 'PMALPHA', 'PMDELTA', 'PMALPHAERR', 'PMDELTAERR',
                     'MAG_AUTO', 'MEDIAN_MAG_AUTO', 'MEAN_MAG_AUTO',
                     'MAGERR_AUTO', 'MEDIAN_MAGERR_AUTO', 'MEAN_MAGERR_AUTO']

        # pm-a-b relation without error
        # new sextractor configuration
        upr_coefs_bright = [5.730178e-03, -6.528091e-01, 2.970536e+01,
                            -6.748832e+02, 7.655361e+03, -3.468232e+04]
        upr_limit_bright = poly1d(upr_coefs_bright)
        lwr_coefs_bright = [5.735456e-03, -6.536302e-01, 2.975250e+01,
                            -6.761721e+02, 7.672408e+03, -3.477074e+04]
        lwr_limit_bright = poly1d(lwr_coefs_bright)

        upr_coefs_faint = [-1.769922e-02, 1.871337e+00, -7.400582e+01,
                           1.297378e+03, -8.505660e+03]
        upr_limit_faint = poly1d(upr_coefs_faint)
        lwr_coefs_faint = [-3.273145e-02, 3.398693e+00, -1.321923e+02,
                           2.282417e+03, -1.475840e+04]
        lwr_limit_faint = poly1d(lwr_coefs_faint)

        filter_tests = {'upr_limit_bright': upr_limit_bright,
                        'lwr_limit_bright': lwr_limit_bright,
                        'upr_limit_faint': upr_limit_faint,
                        'lwr_limit_faint': lwr_limit_faint}

        sub_list_size = len(unique_sources) / self.prfs_d['cores_number']

        sub_list_l = []
        for idx_sub_list in range(0, self.prfs_d['cores_number'], 1):
            if idx_sub_list != (self.prfs_d['cores_number'] - 1):
                idx_down = sub_list_size * idx_sub_list
                idx_up = sub_list_size * (idx_sub_list + 1)
                sub_list_l.append(unique_sources[idx_down:idx_up])
            else:
                idx_down = sub_list_size * idx_sub_list
                sub_list_l.append(unique_sources[idx_down:])

        areas_j = []
        for idx_l in range(0, self.prfs_d['cores_number'], 1):
            areas_p = Process(target=self.filter_b_image_thread,
                              args=(dict_keys, sub_list_l[idx_l],
                                    full_df, filter_tests, idx_l,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges areas
        # Merges catalogs
        csv_list = []
        for idx_csv in range(0, self.prfs_d['cores_number'], 1):
            csv_ = read_csv('{}_8_{}.csv'.format(self.filter_o_n, idx_csv),
                            index_col=0)
            csv_list.append(csv_)

        full_df = concat(csv_list)

        if self.save:
            self.save_message('8')
            full_df.to_csv('{}_8.csv'.format(self.filter_o_n))

        return full_df

    def filter_b_image_thread(self, dict_keys, unique_sources_thread, full_df,
                              filter_tests, idx_l):
        """

        :param dict_keys:
        :param unique_sources_thread:
        :param full_df:
        :param idx_l:
        :return:
        """
        accepted = []
        rejected = []

        print(full_df.columns)

        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources_thread):
            print('filter_pm - thread {} - source {}'.format(idx_l, idx))

            o_df = full_df[full_df['SOURCE_NUMBER'].isin([source_])]

            # b test
            mag_auto = float(o_df['MEDIAN_MAG_AUTO'].iloc[0])
            b_image = float(o_df['MEDIAN_B_IMAGE'].iloc[0])
            pm = float(o_df['PM'].iloc[0])


            if mag_auto < 24.5:
                b_low = filter_tests['lwr_limit_bright'](mag_auto)
                b_upr = filter_tests['upr_limit_bright'](mag_auto)

                if b_low < b_image < b_upr:
                    accepted.append(source_)

                else:
                    rejected.append(source_)
            elif mag_auto > 24.5:
                b_low = filter_tests['lwr_limit_faint'](mag_auto)
                b_upr = filter_tests['upr_limit_faint'](mag_auto)

                if b_low < b_image < b_upr:
                    accepted.append(source_)

                else:
                    rejected.append(source_)

        full_df = full_df[full_df['SOURCE_NUMBER'].isin(accepted)]

        if self.save:
            self.save_message('8_{}'.format(idx_l))
            full_df.to_csv('{}_8_{}.csv'.format(self.filter_o_n, idx_l),
                           columns=dict_keys)

    def filter_coherence(self, full_db):
        """

        :param full_db:
        :return: full_db
        """
        self.logger.debug('Runs coherence motion filter')
        full_db = confidence_filter(full_db, 0.60)  # was 0.97

        return full_db

    def filter_detections(self, full_db, merged_db, detections):
        """

        :param full_db:
        :param merged_db:
        :param detections:
        :return:
        """
        self.logger.debug('Filter by detections number')
        # De momento lo quito
        full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                         if len(g) >= int(detections))
        # Filter by astrometry
        # merged_db = merged_db[merged_db['NPOS_OK'] >= int(detections)]
        # Filter by photometry
        # merged_db = merged_db[merged_db['NMAG'] >= int(detections)]

        # source_list = merged_db['SOURCE_NUMBER'].tolist()
        # full_db = full_db[full_db['SOURCE_NUMBER'].isin(source_list)]

        return full_db, merged_db
