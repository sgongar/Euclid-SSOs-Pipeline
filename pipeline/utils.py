#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import Counter
from multiprocessing import Process

from os import listdir
from sys import argv

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from pandas import concat, DataFrame, Series, read_csv
from numpy import isclose, logical_and, genfromtxt, float64, sqrt

from cats_management import rebase_catalogue, rewriting_catalogue
from misc import setting_logger, extract_settings
from misc import all_same, compare_floats
from regions import Create_regions


def input_ssos(logger, prfs_d):
    """

    :param logger:
    :param prfs_d:
    :return:
    """
    # Gets all possible values of proper motions
    ssos = read_csv('ssos.csv', index_col=0)
    ssos = concat(g for _, g in ssos.groupby('source') if len(g) >= 3)

    unique_sources = list(set(ssos['source'].tolist()))
    total_pms = []

    for source_ in unique_sources:
        df_ssos = ssos[ssos['source'].isin([source_])]
        total_pms.append(float(df_ssos['pm_values'].iloc[0]))

    total_values = dict(Counter(total_pms))
    print('possible values are {}'.format(total_values))

    # Gets right proper motion values
    (input_pm_list, output_pm_list) = modify_table(logger, prfs_d)

    final_input_pm = []
    final_output_pm = []

    for pm_input in input_pm_list:
        pm_output = output_pm_list[input_pm_list.index(pm_input)]
        error = pm_input * 2
        if compare_floats(pm_input, pm_output, error):
            final_input_pm.append(pm_input)
            final_output_pm.append(pm_output)

    for value_ in final_input_pm:
        if value_ == 0.29999999999999999:
            final_input_pm[final_input_pm.index(value_)] = 0.3
        elif value_ == 0.029999999999999999:
            final_input_pm[final_input_pm.index(value_)] = 0.03
        elif value_ == 0.0030000000000000001:
            final_input_pm[final_input_pm.index(value_)] = 0.003

    test = dict(Counter(final_input_pm))

    # Gets obtained proper motion values
    file_sources = read_csv('filter_100_1.2_0.05_20-21_5.csv', index_col=0)

    sources = list(set(file_sources['SOURCE_NUMBER'].tolist()))
    total_pms = []

    """
    for source_ in sources:
        # print source_, type(source_)
        df_sources = file_sources[file_sources['SOURCE_NUMBER'].isin([source_])]
        # print df_sources
        total_pms.append(float(df_sources['PM'].iloc[0]))

    for value_ in total_pms:
        for key_ in total_values.keys():
            print value_, key_
            if isclose([value_], [key_]):
                print value_
    total_values = dict(Counter(total_pms))
    # print 'original values are {}'.format(total_values)
    """


def modify_table(logger, prfs_d):
    """

    @param logger:
    @param prfs_d:

    Important parameters:
    @data_input: a DataFrame object which contains all ssos created
    @data_output: a DataFrame object which contains all ssos detected
                  by sextractor and scamp. These sources were checked against
                  positions obtained from ssos.csv.

    @return pm_input_list, pm_output_list:
    """

    # Gets all sources detected by scamp after filtering process
    post_filter_scamp_cat = 'filter_100_1.2_0.05_20-21_5.csv'
    logger.debug('reading post-filter scamp catalog {}'.format(post_filter_scamp_cat))
    df_sources = read_csv(post_filter_scamp_cat)
    # Gets all unique sources after filter
    logger.debug('getting all unique sources from post-filter scamp catalog')
    valid_sources = list(set(df_sources['SOURCE_NUMBER'].tolist()))

    confidence_list = [2, 5, 10, 20, 50, 100, 200, 500]

    counter = 0
    pm_input_list = []
    pm_output_list = []

    # A dict object ready for contains valuable data
    pm_keys = ['CCD', 'i_alpha', 'i_delta', 'o_alpha', 'o_delta',
               'm_alpha', 'm_delta', 'dither', 'pm_input', 'pm_output',
               'source_scamp', 'source_merged', 'dispersion']
    pm_dict = {}
    for key_ in pm_keys:
        pm_dict[key_] = []

    stats_keys = ['confidence', 'mag', 'pm', 'input_number',
                  'total_scamp', 'right_scamp', 'total_filter', 'right_filter']
    stats_dict = {}
    for key_ in stats_keys:
        stats_dict[key_] = []

    # Gets all sources created
    sources_created_cat = 'ssos.csv'
    logger.debug('reading sources created cat {}'.format(sources_created_cat))
    data_input = read_csv(sources_created_cat, index_col=0)

    # Creates a new DataFrame just to get ssos' speed distribution
    ssos_created_number = concat(g for _, 
                                 g in data_input.groupby('source') if len(g) >= 3)

    unique_ssos = list(set(ssos_created_number['source'].tolist()))
    total_pms = []

    for source_ in unique_ssos:
        df_ssos = ssos_created_number[ssos_created_number['source'].isin([source_])]
        total_pms.append(float(df_ssos['pm_values'].iloc[0]))

    # Gets all sources detected by sextractor and scamp
    sources_detected_cat = 'sso_cat.csv'
    logger.debug('reading sources detected cat {}'.format(sources_detected_cat))
    data_output = read_csv(sources_detected_cat, index_col=0)
    # Gets all unique sources detected. Tip: sources values are refered to
    # sextractor catalog values, scamp detections for each dither shared same
    # catalog number.
    logger.debug('getting all unique sources from detected sources cat')
    unique_sources = list(set(data_output['sources'].tolist()))

    # Gets all sources detected by scamp before filtering process.
    # Useful for obtained scamp calculated proper motion.
    merged_cat_name = '/merged_1.cat'
    logger.debug('reading sources from scamp merged cat {}'.format(merged_cat_name))
    merged_cat = fits.open(prfs_d['results_dir'] + merged_cat_name)
    data_merged = Table(merged_cat[2].data).to_pandas()

    # Looks for sources one by one in sextractor/scamp catalog
    # Try to find if they're an SSO or not
    # source_ number refered to custom catalog, not sextractor/scamp one
    for confidence_ in confidence_list:
        total_values = dict(Counter(total_pms))
        values = sorted(total_values.keys())
        # Harcoded TODO
        mag = '20-21'
        empty_value = 0
        for value_ in values:
            stats_dict['confidence'].append(confidence_)
            stats_dict['mag'].append(mag)
            stats_dict['pm'].append(value_)
            stats_dict['input_number'].append(total_values[value_])
            stats_dict['total_scamp'].append(empty_value)
            stats_dict['right_scamp'].append(empty_value)
            stats_dict['total_filter'].append(empty_value)
            stats_dict['right_filter'].append(empty_value)
        for i, source_ in enumerate(unique_sources):
            flag_input = False
            flag_output = False

            # Gets the input data associated to selected source (luca's data)
            df = data_input[data_input['source'].isin([source_])]
            # Needs the total number of detections of an unique source
            dithers_number = df['dither_values'].size
            # If detections are lower than 3 that source is not interested
            # for our study
            if dithers_number >= 3:
                pm_input = df['pm_values'].iloc[0]  # expected pm
                flag_input = True
            # There are only 2 or less positions availables in our
            # focal plane array.
            else:
                flag_output = False

            # Gets data obtained by sextractor and scamp for selected source
            dt = data_output[data_output['sources'].isin([source_])]
            # Gets references in scamp catalogs for selected source, add them
            # to a list
            same_elements = all_same(dt['scamp_source_number'].tolist())
            # Gets references quantity for selected source
            number_elements = dt['scamp_source_number'].size

            if same_elements and number_elements >= 3 and flag_input:
                # Checks if selected source is present in scamp's post-filter
                # catalog
                if int(dt['scamp_source_number'].iloc[0]) in valid_sources:  # valid_sources = sources in post-filter
                    flag_output = True
                    source_scamp = int(dt['scamp_source_number'].iloc[0])
                    # All detections had the same scamp source number
                    source_merged = dt['scamp_source_number'].iloc[0]

                    # looks in merged catalog so there is only one row per
                    # source
                    dm = data_merged[data_merged['SOURCE_NUMBER'].isin([int(source_merged)])]
                    pm_alpha = dm['PMALPHA_J2000'].iloc[0] / 8.75e6
                    pm_delta = dm['PMDELTA_J2000'].iloc[0] / 8.75e6
                    pm_output = sqrt(pm_alpha**2 + pm_delta**2)
                    if pm_input < pm_output:
                        percentaje_1 = (abs(pm_output - pm_input)) / pm_input
                        dispersion = percentaje_1 * 100.0
                    elif pm_input > pm_output:
                        percentaje_1 = (abs(pm_output - pm_input)) / pm_output
                        dispersion = percentaje_1 * 100.0
                    else:
                        raise Exception

                    if dispersion < confidence_:
                        idx = stats_dict['pm'].index(pm_input)
                        idx = idx + (12 * confidence_list.index(confidence_))
                        stats_dict['right_scamp'][idx] += 1
                        stats_dict['right_filter'][idx] += 1
                else:
                    # Presents in scamp output but not in filtered catalog
                    # Appending to right_scamp
                    source_scamp = int(dt['scamp_source_number'].iloc[0])
                    # All detections had the same scamp source number
                    source_merged = dt['scamp_source_number'].iloc[0]

                    # looks in merged catalog so there is only one row per
                    # source
                    dm = data_merged[data_merged['SOURCE_NUMBER'].isin([int(source_merged)])]
                    pm_alpha = dm['PMALPHA_J2000'].iloc[0] / 8.75e6
                    pm_delta = dm['PMDELTA_J2000'].iloc[0] / 8.75e6
                    pm_output = sqrt(pm_alpha**2 + pm_delta**2)

                    if pm_input < pm_output:
                        percentaje_1 = (abs(pm_output - pm_input)) / pm_input
                        dispersion = percentaje_1 * 100.0
                    elif pm_input > pm_output:
                        percentaje_1 = (abs(pm_output - pm_input)) / pm_output
                        dispersion = percentaje_1 * 100.0

                    if dispersion < confidence_:
                        idx = stats_dict['pm'].index(pm_input)
                        idx = idx + (12 * confidence_list.index(confidence_))
                        stats_dict['right_scamp'][idx] += 1 

                    flag_output = False
            else:
                flag_output = False

            if flag_input and flag_output:
                pm_input_list.append(pm_input)
                pm_output_list.append(pm_output)
                counter = counter + 1
                for d in df['dither_values'].tolist():
                    df_2 = df[df['dither_values'].isin([d])]
                    dt_2 = dt[dt['dithers'].isin([d])]
                    dm_2 = dm
                    dispersion = abs((pm_input / pm_output) - 1) * 100

                    pm_dict['source_scamp'].append(source_scamp)
                    pm_dict['source_merged'].append(source_merged)
                    pm_dict['CCD'].append(str(dt_2['CCD'].iloc[0]))
                    pm_dict['i_alpha'].append(df_2['alpha_j2000'].iloc[0])
                    pm_dict['i_delta'].append(df_2['delta_j2000'].iloc[0])
                    pm_dict['o_alpha'].append(dt_2['scmp_alpha_j2000'].iloc[0])
                    pm_dict['o_delta'].append(dt_2['scmp_delta_j2000'].iloc[0])
                    pm_dict['m_alpha'].append(dm_2['ALPHA_J2000'].iloc[0])
                    pm_dict['m_delta'].append(dm_2['DELTA_J2000'].iloc[0])
                    pm_dict['dither'].append(d)
                    pm_dict['pm_input'].append(pm_input)
                    pm_dict['pm_output'].append(pm_output)
                    pm_dict['dispersion'].append(dispersion)
            else:
                pass

    dict_dataframe = DataFrame(pm_dict)
    dict_dataframe.to_csv('problematic.csv')

    stats_df = DataFrame(stats_dict)
    stats_df = stats_df[['confidence', 'mag', 'pm', 'input_number',
                         'right_scamp', 'total_scamp',
                         'right_filter', 'total_filter']]
    stats_df.to_csv('stats.csv')

    return stats_dict


def extract_trails(logger, prfs_d):
    """

    @param logger:
    @param prfs_d:
    """

    # junta catalogos del mismo tipo
    # agrupa por fuentes groupby source_number
    # mira la velocidad de cada una
    filter_list = ['filter_10_1.2_5_20-21_6.csv',
                   'filter_50_1.2_5_20-21_6.csv',
                   'filter_130_1.2_5_20-21_6.csv']

    for filt in filter_list:
        cat_list = []
        for dither in range(1, 5, 1):
            cat_test = read_csv('{}_{}_{}.csv'.format(dither, filt, dither))
            cat_list.append(cat_test)
        df = concat(cat_list)  # remove old index

        db = concat(g for _, g in df.groupby("SOURCE_NUMBER") if len(g) >= 2)
        unicos = db['SOURCE_NUMBER'].unique()

        df.to_csv(filt)

    return True


def get_input_cat(logger, prfs_d, cat_file):
    """

    @param logger:
    @param prfs_d:
    @param cat_file:

    @return True: if everything goes alright.
    """

    catalogue = genfromtxt(cat_file)

    list_x = catalogue[:, 0]
    list_y = catalogue[:, 1]
    list_mag = catalogue[:, 2]
    list_pm = catalogue[:, 3]

    speed_0_001 = range(prfs_d['first_sso'], 137447, 75)
    speed_0_003 = range(prfs_d['first_sso'] + 10, 137457, 75)
    speed_0_01 = range(prfs_d['first_sso'] + 20, 137467, 75)
    speed_0_03 = range(prfs_d['first_sso'] + 30, 137477, 75)
    speed_0_1 = range(prfs_d['first_sso'] + 40, 137487, 75)
    speed_0_3 = range(prfs_d['first_sso'] + 50, 137497, 75)
    speed_1 = range(prfs_d['first_sso'] + 60, 137507, 75)

    speed_3 = range(prfs_d['first_sso'] + 66, 137512, 75)
    speed_10 = range(prfs_d['first_sso'] + 67, 137513, 75)
    speed_30 = range(prfs_d['first_sso'] + 68, 137514, 75)
    speed_100 = range(prfs_d['first_sso'] + 69, 137515, 75)
    speed_300 = range(prfs_d['first_sso'] + 70, 137516, 75)

    for index in speed_0_001:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 0.001
    for index in speed_0_003:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 0.003
    for index in speed_0_01:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 0.01
    for index in speed_0_03:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 0.03
    for index in speed_0_1:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 0.1
    for index in speed_0_3:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 0.3
    for index in speed_1:
        list_mag[index] = list_mag[index] - 2.5
        list_pm[index] = 1
    for index in speed_3:
        list_pm[index] = list_pm[index] - 1000
    for index in speed_10:
        list_pm[index] = list_pm[index] - 1000
    for index in speed_30:
        list_pm[index] = list_pm[index] - 1000
    for index in speed_100:
        list_pm[index] = list_pm[index] - 1000
    for index in speed_300:
        list_pm[index] = list_pm[index] - 1000

    indexes = (speed_0_001 + speed_0_003 + speed_0_01 + speed_0_03 +
               speed_0_1 + speed_0_3 + speed_1 + speed_3 + speed_10 +
               speed_30 + speed_100 + speed_300)
    indexes = sorted(indexes)

    s1 = Series(list_x, name='X_IMAGE', dtype=float64)
    s2 = Series(list_y, name='Y_IMAGE', dtype=float64)
    s3 = Series(list_mag, name='mag_values', dtype=float64)
    s4 = Series(list_pm, name='pm_values', dtype=float64)

    sources_df = concat([s1, s2, s3, s4], axis=1)
    sources_df = sources_df.iloc[indexes, :]

    return sources_df


def check(logger, prfs_d):
    """

    @param logger:
    @param prfs_d:
    """

    cat_files = []
    cat_files_d = listdir(prfs_d['results_dir'])
    for cat_file in cat_files_d:
        if cat_file[-5:] == '6.csv':
            cat_files.append(prfs_d['results_dir'] + '/' + cat_file)

    i_cats = {}
    i_cat_files = listdir(prfs_d['input_cats'])
    # Read  input catalog files
    for cat in i_cat_files:
        # All dithers, catalog should start with 'mag' string
        if cat[-4:] == '.dat':
            for mag in prfs_d['mags']:
                if cat[-12:-7] == mag:
                    cat_name = prfs_d['input_cats'] + '/' + cat
                    print(cat_name)
                    """
                    catalog = read_csv(cat_name,
                                       names=['X_IMAGE', 'Y_IMAGE', 'ID_1',
                                              'ID_2', 'ID_3'], delimiter=" ",
                                       skiprows=range(1, 134893))
                    """
                    catalog = get_input_cat(logger, prfs_d, cat_name)
                    # It's reading everything!
                    i_cats['{}'.format(int(cat[-5:-4]))] = catalog

    fits_loc = '{}/m_{}_x0_y0_d1_copy.fits'.format(prfs_d['fits_dir'], mag)
    hdulist = fits.open(fits_loc)
    w = WCS(hdulist[0].header)

    i_cats_t = {}

    for dither in range(1, 5, 1):
        input_table = i_cats[str(dither)]
        regions_list = []
        for source_num in range(input_table['X_IMAGE'].as_matrix().size):
            x_value = input_table['X_IMAGE'].as_matrix()[source_num]
            y_value = input_table['Y_IMAGE'].as_matrix()[source_num]
            regions_list.append([x_value, y_value])

        logger.debug('transforming x/y values to ra/dec ones')
        input_regions = w.wcs_pix2world(regions_list, 1)

        logger.debug('creating new list with ra/dec values')
        alpha_list = []
        delta_list = []
        for regions in input_regions:
            alpha_list.append(regions[0])
            delta_list.append(regions[1])

        mag_values_list = input_table['mag_values'].tolist()
        pm_values_list = input_table['pm_values'].tolist()

        mag_values = Series(mag_values_list, name='mag_values')
        pm_values = Series(pm_values_list, name='pm_values')

        logger.debug('creating new Pandas series for ra/dec values')
        alpha_j2000 = Series(alpha_list, name='alpha_j2000')
        delta_j2000 = Series(delta_list, name='delta_j2000')

        logger.debug('concatenating series of ra and dec')
        coords_table = concat([alpha_j2000, delta_j2000,
                               mag_values, pm_values], axis=1)

        coords_table.to_csv('dither_{}.csv'.format(dither), index=False,
                            header=False, sep=" ")

        i_cats_t[dither] = coords_table

    cat_n = {1: [1, 2, 3, 4, 5, 6, 7, 8, 9],
             2: [10, 11, 12, 13, 14, 15, 16, 17, 18],
             3: [19, 20, 21, 22, 23, 24, 25, 26, 27],
             4: [28, 29, 30, 31, 32, 33, 34, 35, 36]}

    o_cats = {}
    for cat_file in cat_files:
        # cat_file[70:] id!
        filter_n = cat_file[70:]  # harcoded!
        cat = read_csv(cat_file, index_col=0)
        for dither in range(1, 5, 1):
            cat_d = cat[cat['CATALOG_NUMBER'].isin(cat_n[dither])]

            o_cats['{}_{}'.format(filter_n, dither)] = cat_d

    check_j = []
    for proc in range(1, 5, 1):
        # Each dither in a different thread

        temp_keys = []

        for key in o_cats.keys():
            if key[-1:] == str(proc):
                temp_keys.append(key)

        o_cat_t = {new_key: o_cats[new_key] for new_key in temp_keys}

        check_p = Process(target=check_thread,
                          args=(logger, prfs_d, i_cats_t[proc],
                                o_cat_t, proc,))
        check_j.append(check_p)
        check_p.start()

    active_check = list([job.is_alive() for job in check_j])
    while True in active_check:
        active_check = list([job.is_alive() for job in check_j])
        pass


def check_thread(logger, prfs_d, i_cats_t, o_cats, dither):
    """

    @param logger:
    @param prfs_d:
    @param i_cats_t:
    @param o_cats:

    @return True:
    """

    # counts
    stats_dict = {}
    for d in range(1, 5, 1):
        for o_cat in o_cats.keys():
            stats_dict['{}_{}'.format(dither, o_cat)] = {}

    columns = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION', 'ASTR_INSTRUM',
               'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE', 'ERRA_IMAGE',
               'ERRB_IMAGE', 'ERRTHETA_IMAGE', 'ALPHA_J2000', 'DELTA_J2000',
               'ERRA_WORLD', 'ERRB_WORLD', 'ERRTHETA_WORLD', 'EPOCH', 'MAG',
               'MAGERR', 'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'PM',
               'PMERR', 'PMALPHA', 'PMDELTA', 'PMALPHAERR', 'PMDELTAERR']
    for key_ in stats_dict.keys():
        for column in columns:
            stats_dict[key_][column] = []

    # print "o_cats", o_cats
    # print "o_cats.keys_()", o_cats.keys()

    sources_number = i_cats_t['alpha_j2000'].size  # same as delta_j2000.size
    for source_num in range(sources_number):
        # voy mirando una a una las fuentes de sextractor y si esta en el
        # original copio todo a una table nueva
        # getting alpha, delta for each source sextracted
        # FIXME no uso valores absolutos!
        source_alpha = i_cats_t['alpha_j2000'].iloc[source_num]
        source_delta = i_cats_t['delta_j2000'].iloc[source_num]

        for o_cat in o_cats.keys():
            alpha_mask = isclose(o_cats[o_cat]['ALPHA_J2000'].values[:, None],
                                 [source_alpha],
                                 atol=prfs_d['alpha_t']).any(axis=1)
            delta_mask = isclose(o_cats[o_cat]['DELTA_J2000'].values[:, None],
                                 [source_delta],
                                 atol=prfs_d['delta_t']).any(axis=1)

            mask = logical_and(alpha_mask, delta_mask)
            df = o_cats[o_cat][mask]
            if df.size == 27:
                dict_name = '{}_{}'.format(dither, o_cat)
                for column in columns:
                    stats_dict[dict_name][column].append(df[column].iloc[0])

    series_dict = {}
    for key_ in stats_dict.keys():
        series_dict[key_] = {}
        for column in columns:
            series_dict[key_][column] = Series(stats_dict[key_][column],
                                               name=column)

    df_dict = {}
    for key_ in series_dict.keys():
        tmp_list = []
        for column in columns:
            tmp_list.append(series_dict[key_][column])
        df_dict[key_] = concat(tmp_list, axis=1)

    for key_ in df_dict.keys():
        df_dict[key_].to_csv(key_ + '.csv', index=False)

    return True


def change_times(logger, prfs_d, mag):
    """
    hardcoded times:
    2021-06-26T09:00:00.00000
    2021-06-26T09:16:43.00000
    2021-06-26T09:33:26.00000
    2021-06-26T09:50:09.00000

    hardcoded dir:

    @param logger:
    @param prfs_d:

    """
    core_number = int(prfs_d['cores_number'])
    fits_dir = '{}/{}/CCDs'.format(prfs_d['fits_dir'], mag)
    files = listdir(fits_dir)
    fits_files = []

    for file in files:
        if file[-5:] == '.fits':
            fits_files.append(file)

    for idx_fits in range(0, len(fits_files), core_number):
        try:
            time_j = []
            for proc in range(0, core_number, 1):
                idx = idx_fits + proc
                time_p = Process(target=change_times_thread,
                                 args=(logger, prfs_d, fits_files[idx],))
                time_j.append(time_p)
                time_p.start()

            active_time = list([job.is_alive() for job in time_j])
            while True in active_time:
                active_time = list([job.is_alive() for job in time_j])
                pass
        except IndexError:
            logger.debug('extraction finished')

    return True


def change_times_thread(logger, prfs_d, fits_image):
    """

    @param logger:
    @param prfs_d:
    @param fits_image:

    @return True:
    """
    fits_dir = '{}/{}/CCDs'.format(prfs_d['fits_dir'], mag)

    data, header = fits.getdata('{}/{}'.format(fits_dir, fits_image),
                                header=True)
    dither = fits_image[-6:-5]

    if dither == '1':
        header['DATE-OBS'] = prfs_d['time_1']
        # Change format
        t = Time(prfs_d['time_1'])
        t.format = 'mjd'
        header['MJD-OBS'] = float(str(t))
    elif dither == '2':
        header['DATE-OBS'] = prfs_d['time_2']
        t = Time(prfs_d['time_2'])
        t.format = 'mjd'
        header['MJD-OBS'] = float(str(t))
    elif dither == '3':
        header['DATE-OBS'] = prfs_d['time_3']
        t = Time(prfs_d['time_3'])
        t.format = 'mjd'
        header['MJD-OBS'] = float(str(t))
    elif dither == '4':
        header['DATE-OBS'] = prfs_d['time_4']
        t = Time(prfs_d['time_4'])
        t.format = 'mjd'
        header['MJD-OBS'] = float(str(t))
    else:
        header['DATE-OBS'] = prfs_d['time_1']
        # Change format
        t = Time(prfs_d['time_1'])
        t.format = 'mjd'
        header['MJD-OBS'] = float(str(t))

    print('opens {} date-obs {} mjd-obs {} d {}'.format(fits_image,
                                                        header['DATE-OBS'],
                                                        header['MJD-OBS'],
                                                        dither))

    fits.writeto('{}/{}_copy'.format(fits_dir, fits_image),
                 data, header, clobber=True)

    return True


def rebase_thread(logger, prfs_d, mag):
    """

    @param logger:
    @param prfs_d:
    @param mag:

    """
    # TODO create start catalog
    input_cat = '{}/catalogue_{}.cat'.format(prfs_d['output_cats'], mag)
    logger.debug('input_catalogue is {}'.format(input_cat))
    remove_ssos = True
    final_catalogue = rebase_catalogue(logger, mag, prfs_d,
                                       remove_ssos, input_cat)

    if not rewriting_catalogue(logger, prfs_d,
                               final_catalogue, mag):
        raise Exception


if __name__ == '__main__':
    logger = setting_logger()
    prfs_d = extract_settings()

    try:
        if argv[1] == '-change':
            for mag in prfs_d['mags']:
                print('mag {}'.format(mag))
                change_times(logger, prfs_d, mag)
        elif argv[1] == '-rebase':
            cores_number = prfs_d['cores_number']
            if cores_number > len(prfs_d['mags']):
                cores_number = len(prfs_d['mags'])
            rebase_j = []
            for proc in range(0, cores_number, 1):
                rebase_p = Process(target=rebase_thread,
                                   args=(logger, prfs_d,
                                         prfs_d['mags'][proc]))
                rebase_j.append(rebase_p)
                rebase_p.start()

            active_rebase = list([job.is_alive() for job in rebase_j])
            while True in active_rebase:
                active_rebase = list([job.is_alive() for job in rebase_j])
                pass
        elif argv[1] == '-check':
            check(logger, prfs_d)
        elif argv[1] == '-full_cats':
            i_c1 = '/mnt/c/CCDs'
            i_c2 = '/full_10_1.2_5_1_20-21_1.cat'
            i_c = i_c1 + i_c2

            logger.debug('opening catalog file {}'.format(i_c))
            Create_regions(i_c).full_cats()
        elif argv[1] == '-full_regions':

            i_c1 = '/mnt/c/CCDs'
            i_c2 = '/full_10_1.2_5_1_20-21_1.cat'
            i_c = i_c1 + i_c2

            logger.debug('opening catalog file {}'.format(i_c))
            Create_regions(i_c).full_regions()
        elif argv[1] == '-regions':
            opts = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1],
                    [1, 2], [2, 0], [2, 1], [2, 2]]
            for opt_ in opts:
                mag = '20-21'

                sex_folders = ['30_1.5_1.5_0.01_4']
                cat_dir = prfs_d['fits_dir']
                for folder_ in sex_folders:
                    for d in range(1, 5, 1):
                        cat_name = 'mag_{}_CCD_x{}_y{}_d{}.cat'.format(mag,
                                                                       opt_[0],
                                                                       opt_[1],
                                                                       d)
                        i_c = '{}/{}/CCDs/{}/{}'.format(cat_dir, mag, folder_,
                                                        cat_name)
                        Create_regions(i_c).fits()
                        logger.debug('opening catalog file {}'.format(i_c))
        elif argv[1] == '-luca_check':
            input_d = {}
            for d in range(1, 5, 1):
                input_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                           d)
            Create_regions(input_d).check_luca(True, True)
        elif argv[1] == '-modify':
            list_1, list_2 = modify_table(logger, prfs_d)
        elif argv[1] == '-input':
            input_ssos(logger, prfs_d)
        else:
            "Wrong choice"
    except Exception as e:
        print(e)
        logger.error(e)
        logger.error('Select a valid option')


def input_regions():
    """

    :return:
    """
    prfs_d = extract_settings()

    input_d = {}
    for d in range(1, 5, 1):
        input_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'], d)
    input_d = Create_regions(input_d, prfs_d).check_luca(True, True)

    # Creates a DataFrame from an input dictionary
    input_l = []
    for key_ in input_d.keys():
        input_l.append(input_d[key_])

    i_df = concat(input_l, axis=0)
    # Look for < 3 coincidences
    i_df = concat(g for _, g in i_df.groupby('source')
                  if len(g) >= 3)
    i_df = i_df.reset_index(drop=True)

    # Saves sources
    i_df.to_csv('input_sources.csv')

# from astropy.io import fits
# from astropy.time import Time
#
# time_1 = '2021-06-26T09:00:00.00000'
# fits_file = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00.fits'
# for idx in range(0, 36, 1):
#     print(idx)
#     """
#     fits_name = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_{}.fits'.format(idx)
#     data, header = fits.getdata(fits_name, header=True)
#     header['DATE-OBS'] = time_1
#     # Change format
#     t = Time(time_1)
#     t.format = 'mjd'
#     header['MJD-OBS'] = float(str(t))
#
#     fits.writeto('{}_copy'.format(fits_name), data, header)
#
#     fits_name = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_f{}.fits'.format(idx)
#     data, header = fits.getdata(fits_name, header=True)
#     header['DATE-OBS'] = time_1
#     # Change format
#     t = Time(time_1)
#     t.format = 'mjd'
#     header['MJD-OBS'] = float(str(t))
#
#     fits.writeto('{}_copy'.format(fits_name), data, header)
#     """
#
#     fits_name = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_w{}.fits'.format(idx)
#     data, header = fits.getdata(fits_name, header=True)
#     header['DATE-OBS'] = time_1
#     # Change format
#     t = Time(time_1)
#     t.format = 'mjd'
#     header['MJD-OBS'] = float(str(t))
#
#     fits.writeto('{}_copy'.format(fits_name), data, header)