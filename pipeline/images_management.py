#!/usr/bin/python
# -*- coding: utf-8 -*-

from logging import getLogger
from os import listdir

from astropy.io import fits
from astropy.wcs import WCS
from pandas import concat, Series, read_csv
from multiprocessing import Process, Manager
from numpy import isclose, logical_and

# from cats_management import get_output_catalogue


def merge_results(logger, prfs_d, stats_dict):
    """

    @param logger:
    @param psfs_d:
    @param stats_dict:

    @return stats_dict:
    """

    full_stats = []
    for key in stats_dict.keys():
        full_stats.append(stats_dict['{}'.format(key)])

    full_stats = concat(full_stats)

    list_mag = []
    list_dither = []
    list_sources = []
    list_detected = []
    list_right = []
    list_false = []
    list_f_pur = []
    list_f_com = []
    list_mincount = []
    list_threshold = []
    list_filter = []

    mags = ['20-21', '21-22']  # test reasons

    # for mag in prfs_d['mags']:
    for mag in mags:
        # logger.debug('looking for magnitude {}'.format(mag))
        df_1 = full_stats[full_stats['mag'].isin([mag])]
        for dither in range(1, 5, 1):
            # logger.debug('looking for dither {}'.format(dither))
            df_2 = df_1[df_1['dither'].isin(['d' + str(dither)])]

            list_mag.append(df_2['mag'].at[0][0])
            list_dither.append(df_2['dither'].at[0][0])
            if df_2['ccd'].size != 4:
                print "df_2 es", df_2
                raise Exception
            sources = df_2['sources'].at[0][0]
            list_sources.append(df_2['sources'].at[0][0])
            detected = df_2['detected'].sum()  # detected sources of all CCDs
            list_detected.append(detected)
            right = df_2['right'].sum()  # right sources of all CCDs
            list_right.append(right)
            wrong = df_2['wrong'].sum()  # wrong sources of all CCDs
            list_false.append(wrong)
            f_dr = float(detected) / float(sources)
            f_pur = float(right) / float(detected)
            f_com = f_dr * f_pur
            list_f_pur.append(f_pur)
            list_f_com.append(f_com)
            list_mincount.append(df_2['mincount'].at[0][0])
            list_threshold.append(df_2['threshold'].at[0][0])
            list_filter.append(df_2['filter'].at[0][0][7:])

    s1 = Series(list_mag, name='mag')
    s2 = Series(list_dither, name='dither')
    s3 = Series(list_sources, name='sources')
    s4 = Series(list_detected, name='detected')
    s5 = Series(list_right, name='right')
    s6 = Series(list_false, name='wrong')
    s7 = Series(list_f_pur, name='purity')
    s8 = Series(list_f_com, name='completeness')
    s9 = Series(list_mincount, name='mincount')
    s10 = Series(list_threshold, name='threshold')
    s11 = Series(list_filter, name='filter')

    stats = concat([s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11], axis=1)
    print "final stats is", stats
    stats.to_csv('finals.csv')

    return stats_dict


def check(logger, prfs_d, analysis_d):
    """

    @param logger:
    @param prfs_d:
    @param analysis_d:

    """
    log_check = getLogger("main.check")

    log_check.debug('input catalogs from {}'.format(prfs_d['input_cats']))
    i_cats = []
    i_cat_files = listdir(prfs_d['input_cats'])
    # Read  input catalog files
    for cat in i_cat_files:
        # All dithers, catalog should start with 'mag' string
        if cat[-4:] == '.dat':
            for mag in prfs_d['mags']:
                if cat[-12:-7] == mag:
                    i_cats.append(cat)

    # logger.debug('input cats are {}'.format(i_cats))

    log_check.debug('output catalogs from {}'.format(prfs_d['output_cats']))
    o_cats = []
    o_cat_files = listdir(prfs_d['output_cats'])
    # Read output catalog files
    for cat in o_cat_files:
        # All catalogs are chosen
        if cat[-4:] == '.cat' and cat[-31:-28] == 'mag':
            o_cats.append(cat)

    # ok!
    # i_cat_dict = {}
    cat_number = prfs_d['version'][-1:]  # cat version got from settings

    core_number = 8
    for cat_idx in range(0, len(i_cats), core_number):
        try:
            cat_j = []
            for proc in range(0, core_number, 1):
                idx = cat_idx + proc  # index
                cat = i_cats[idx]

                cat_p = Process(target=check_thread_i,
                                args=(logger, prfs_d, cat, cat_number,))
                cat_j.append(cat_p)
                cat_p.start()

            active_cat = list([job.is_alive() for job in cat_j])
            while True in active_cat:
                active_cat = list([job.is_alive() for job in cat_j])
                pass
        except IndexError:
            logger.debug('extraction finished')

    # Once input data has been extracted from input catalogs data sextracted
    # from images should be analysed

    o_cat_dict = {}
    # Extracts data from output catalogs
    for cat in o_cats:
        # Creates a dict's entry for each mag, dither and position
        cat = prfs_d['output_cats'] + '/' + cat
        print "cat es", cat
        mag = cat[-27:-22]
        dither = cat[-11:-9]
        x_y = cat[-17:-12]
        cat_data = get_output_catalogue(logger, cat)
        o_cat_dict['{}_{}_{}'.format(mag, dither, x_y)] = cat_data

    manager = Manager()  # Creates a Manager object
    stats_dict = manager.dict()  # Creates a shared dict to store results
    for cat_idx in range(0, len(o_cat_dict.keys()), core_number):
        try:
            cat_j = []
            for proc in range(0, int(prfs_d['cores_number']), 1):
                idx = cat_idx + proc  # index
                # read input catalogue associated to output one
                print "catalogo es ", o_cat_dict.keys()[idx]
                mag = o_cat_dict.keys()[idx][0:5]
                dither = o_cat_dict.keys()[idx][6:8]
                ccd = o_cat_dict.keys()[idx][9:14]
                logger.debug('reading cat {}'.format(o_cat_dict.keys()[idx]))
                cat_p = Process(target=check_thread_o,
                                args=(logger, prfs_d, mag, dither,
                                      ccd, o_cat_dict, idx,
                                      analysis_d, stats_dict,))
                cat_j.append(cat_p)
                cat_p.start()

            active_cat = list([job.is_alive() for job in cat_j])
            while True in active_cat:
                active_cat = list([job.is_alive() for job in cat_j])
                pass
        except IndexError:
            logger.debug('extraction finished')

    stats_dict = merge_results(logger, prfs_d, stats_dict)

    return stats_dict


def check_thread_i(log_check, prfs_d, cat, cat_number):
    """

    @param logger:
    @param prfs_d:
    @param cat:
    @param cat_number:

    """
    cat = prfs_d['input_cats'] + '/' + cat
    mag = cat[-12:-7]  # referenced to the end of the file!
    dither = cat[-6:-4]

    complete_list = True

    log_check.debug('reading catalog {}'.format(cat))
    cat_data = get_input_cat(log_check, prfs_d, cat_number,
                             cat, complete_list)

    # Convert xy data to ra/dec
    # Get WCS data from first dither
    ccd_wcs = '{}/mag_{}_CCD_x0_y0_d1.fits'.format(prfs_d['fits_dir'], mag)
    hdulist = fits.open(ccd_wcs)
    w = WCS(hdulist[0].header)

    regions_list = []
    for source_num in range(cat_data['x_values'].as_matrix().size):
        x_value = cat_data['x_values'].as_matrix()[source_num]
        y_value = cat_data['y_values'].as_matrix()[source_num]
        regions_list.append([x_value, y_value])

    # Transform x/y positions related to first dither to ra/dec positions
    log_check.debug('getting ra/dec coordinates from {}'.format(ccd_wcs))
    input_regions = w.wcs_pix2world(regions_list, 1)

    # Creates a series of list with ra/dec positions, magnitudes and proper
    # motions and populates them
    ra_values = []
    dec_values = []
    mag_values = []
    pm_values = []

    for i in range(cat_data['x_values'].as_matrix().size):
        ra_values.append(input_regions[i][0])
        dec_values.append(input_regions[i][1])
        mag_values.append(cat_data['mag_values'].iloc[i])
        pm_values.append(cat_data['pm_values'].iloc[i])

    # Creates differente Series for ra/dec positions, mangnitudes and motions
    log_check.debug('creating series for ra/dec, mag and pm values')
    c1 = Series(ra_values, name='ra_values')
    c2 = Series(dec_values, name='dec_values')
    c3 = Series(mag_values, name='mag_values')
    c4 = Series(pm_values, name='pm_values')

    log_check.debug('concatenig data to pandas DataFrame')
    cat_data = concat([c1, c2, c3, c4], axis=1)  # Join series into a DataFrame

    # Read columns and transform them
    log_check.debug('saving data to .csv file')
    cat_data.to_csv('{}/{}_{}.csv'.format(prfs_d['fits_dir'],
                                          mag, dither))  # test reasons!

    return True


def check_thread_o(logger, prfs_d, mag, dither, ccd, o_cat_dict,
                   cat, analysis_d, stats_dict):
    """

    @param logger:
    @param prfs_d:
    @param mag:
    @param dither:
    @param o_cat_dict:
    @param cat:

    @return True: if everything goes alright
    """
    # Stuff need for comparation

    input_name = '{}/{}_{}.csv'.format(prfs_d['fits_dir'],
                                       mag, dither)
    logger.debug('looking to {} input cat'.format(input_name))
    input_cat = read_csv(input_name)
    logger.debug('looking to {} output cat'.format(o_cat_dict.keys()[cat]))
    output_cat = o_cat_dict[o_cat_dict.keys()[cat]]

    ok = 0

    # Creates lists for Series creation
    list_mag = []
    list_dither = []
    list_ccd = []
    list_sources = []
    list_detected = []
    list_ok = []
    list_no = []
    list_mincount = []
    list_threshold = []
    list_filter = []

    sources_number = input_cat['ra_values'].size
    logger.debug('{} sources to analyse in {}'.format(sources_number,
                                                      input_name))
    for source_num in range(sources_number):
        source_alpha = input_cat['ra_values'].iloc[source_num]
        source_delta = input_cat['dec_values'].iloc[source_num]

        alpha_mask = isclose(output_cat['ALPHA_J2000'].values[:, None],
                             [source_alpha],
                             atol=float(prfs_d['alpha_t'])).any(axis=1)
        delta_mask = isclose(output_cat['DELTA_J2000'].values[:, None],
                             [source_delta],
                             atol=float(prfs_d['delta_t'])).any(axis=1)

        mask = logical_and(alpha_mask, delta_mask)
        df = output_cat[mask]

        if not df.empty:
            ok += 1

    list_mag.append(mag)
    list_dither.append(dither)
    list_ccd.append(ccd)
    list_sources.append(sources_number)
    sources_detected = output_cat['ALPHA_J2000'].size
    list_detected.append(sources_detected)
    list_ok.append(ok)
    no = sources_detected - ok
    list_no.append(no)
    list_mincount.append(analysis_d['deblend_mincount'])
    list_threshold.append(analysis_d['analysis_thresh'])
    list_filter.append(analysis_d['filter'])

    s1 = Series(list_mag, name='mag')
    s2 = Series(list_dither, name='dither')
    s3 = Series(list_ccd, name='ccd')
    s4 = Series(list_sources, name='sources')
    s5 = Series(list_detected, name='detected')
    s6 = Series(list_ok, name='right')
    s7 = Series(list_no, name='wrong')
    s8 = Series(list_mincount, name='mincount')
    s9 = Series(list_threshold, name='threshold')
    s10 = Series(list_filter, name='filter')

    stats = concat([s1, s2, s3, s4, s5, s6, s7, s8, s9, s10], axis=1)

    logger.debug('saving file {}_{}_{}.csv'.format(mag, dither, ccd))
    stats.to_csv('{}/pipeline/stats_{}_{}_{}.csv'.format(prfs_d['home'],
                                                         mag, dither, ccd))
    stats_dict['{}_{}'.format(mag, dither)] = stats

    return True


def get_fits_limits(fits_image):
    """

    @param logger:
    @param fits_image: fits image

    @return limits: a dict with ra/dec limits above_ra, below_ra,
                    above_dec_, below_dec
    """
    # logger.info('getting limits of {} image'.format(fits_image))

    data, header = fits.getdata(fits_image, header=True)
    w = WCS(fits_image)

    above_x, above_y = header['NAXIS1'], header['NAXIS2']
    above_ra, above_dec = w.all_pix2world(above_x, above_y, 0)

    below_ra, below_dec = w.all_pix2world(0, 0, 0)

    limits = {'below_ra': float(above_ra), 'above_ra': float(below_ra),
              'below_dec': float(below_dec), 'above_dec': float(above_dec)}

    # check position
    # sometimes some values could be higher when are tagged as "lowest"
    return limits
