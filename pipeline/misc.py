#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements


Todo:
    * Improve log messages
    * Improve usability
"""

from collections import Counter
from decimal import Decimal
from math import hypot
from multiprocessing import cpu_count
import os
import platform

from ConfigParser import ConfigParser
import numpy as np
from pandas import concat, Series
import statsmodels.api as sm

from errors import BadSettings, AllSameException, WrongOS
from errors import InvalidScampConfiguration
from logging import getLogger, config


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_cats():
    """ Loops over the different folders.

    :return:
    """
    prfs_d = extract_settings()

    mode = {'type': 'sextractor'}
    confs, total_confs = create_configurations(mode)

    folders = []
    for idx, conf_ in enumerate(confs):
        analysis_d, len_dicts = create_sextractor_dict(idx, False)
        folder_n = '{}_{}_{}_{}_{}'.format(analysis_d['deblend_nthresh'],
                                           analysis_d['analysis_thresh'],
                                           analysis_d['detect_thresh'],
                                           analysis_d['deblend_mincount'],
                                           analysis_d['detect_minarea'])
        folders.append(folder_n)

    cat_dict = {}
    # FIXME Magnitude it's harcoded...change it!
    for folder_ in folders:
        cat_dict[folder_] = []
        files = os.listdir('{}/20-21/CCDs/{}'.format(prfs_d['fits_dir'],
                                                     folder_))
        for file_ in files:
            if file_[:1] == 'm' and file_[-4:] == '.cat':
                cat_dict[folder_].append(file_)

    return cat_dict


def get_os():
    """ a function that gets the current operative system
    for now works in Debian, Fedora and Ubuntu shell (Microsoft version)

    @return os_system: a string which contains the operative system name
    """

    if 'fedora-23' in platform.platform():
        os_system = 'test'
    elif 'Debian' in platform.platform():
        os_system = 'debian'
    elif 'Ubuntu' in platform.platform():
        os_system = 'ubuntu'
    elif 'fedora-26' in platform.platform():
        os_system = 'fedora'
    elif 'fedora-19' in platform.platform():
        os_system = 'cab'
    elif 'centos' in platform.platform():
        os_system = 'centos'
    else:
        raise WrongOS

    return os_system


def all_same(items):

    # If there are more than two different values by definition should
    # be wrong

    length_items = len(list(set(items)))
    items_w_o_false = [x for x in items if x != 'False']

    if length_items is 1 and 'False' not in dict(Counter(items)).keys():
        return True, len(items_w_o_false)
    elif length_items is 1 and 'False' in dict(Counter(items)).keys():
        return False, 0
    elif length_items is 2 and 'False' in dict(Counter(items)).keys():
        if len(items) == 4:
            if dict(Counter(items))['False'] > 1:
                return False, 0
            elif dict(Counter(items))['False'] is 1:
                return True, 3  # Harcoded?
        else:
            return False, 0
    elif length_items > 2:
        return False, 0
    else:
        raise AllSameException


def create_sextractor_dict(conf_num, cat_conf):
    """

    @param conf_num:
    @param cat_conf:

    @return analysis_d:
    """

    if cat_conf:
        configurations = [2, 0.1, 5, 4, 'models/gauss_2.0_5x5.conv']
        len_conf = 1
    else:
        mode = {'type': 'sextractor'}  # harcoded
        configurations, len_conf = create_configurations(mode)

    analysis_l = []

    if type(configurations[0]) is list:
        for configuration in range(len(configurations)):
            temp_list = [configurations[configuration][0],
                         configurations[configuration][1],
                         configurations[configuration][2],
                         configurations[configuration][2],
                         configurations[configuration][3],
                         configurations[configuration][4]]
            analysis_l.append(temp_list)
        analysis_d = {'deblend_nthresh': analysis_l[conf_num][0],
                      'deblend_mincount': analysis_l[conf_num][1],
                      'detect_thresh': analysis_l[conf_num][2],
                      'analysis_thresh': analysis_l[conf_num][3],
                      'detect_minarea': analysis_l[conf_num][4],
                      'filter': analysis_l[conf_num][5]}
    else:
        analysis_d = {'deblend_nthresh': configurations[0],
                      'deblend_mincount': configurations[1],
                      'detect_thresh': configurations[2],
                      'analysis_thresh': configurations[2],
                      'detect_minarea': configurations[3],
                      'filter': configurations[4]}

    return analysis_d, len_conf


def create_scamp_dict(conf_num):
    """

    :param conf_num:
    :return:
    """

    scamp_list = []
    mode = {'type': 'scamp'}
    configurations, len_conf = create_configurations(mode)

    for conf in configurations:
        temp_list = [conf[0], conf[1], conf[2], conf[3]]
        scamp_list.append(temp_list)

    try:
        scamp_dict = {'crossid_radius': scamp_list[conf_num][0],
                      'pixscale_maxerr': scamp_list[conf_num][1],
                      'posangle_maxerr': scamp_list[conf_num][2],
                      'position_maxerr': scamp_list[conf_num][3]}
    except IndexError:
        raise InvalidScampConfiguration

    return scamp_dict, len_conf


def create_configurations(mode):
    """ creates a list of configuration lists merging
        different input parameters
    :param mode: can be 'sextractor' or 'scamp'
    :return:
    """
    if mode['type'] == 'sextractor':
        l_deblending = [30]
        l_mincount = [0.01]
        l_threshold = [1.5]

        l_area = [4]
        l_filter_name = ['models/gauss_2.0_5x5.conv']

        configurations = []
        for deblending in l_deblending:
            for mincount in l_mincount:
                for threshold in l_threshold:
                    for area in l_area:
                        for filt in l_filter_name:
                            configurations.append([deblending, mincount,
                                                   threshold, area, filt])
        configurations_len = len(configurations)
        return configurations, configurations_len

    elif mode['type'] == 'scamp':
        l_crossid_radius = [10]  # [10] seconds
        l_pixscale_maxerr = [1.1]  # [1.2] scale-factor
        l_posangle_maxerr = [0.5]  # [0.5, 2.5] degrees
        l_position_maxerr = [0.04]

        configurations = []

        for crossid in l_crossid_radius:
            for pixscale in l_pixscale_maxerr:
                for posangle in l_posangle_maxerr:
                    for position in l_position_maxerr:
                        configurations.append([crossid, pixscale,
                                               posangle, position])

        configurations_len = len(configurations)

        return configurations, configurations_len


def confmap(config_, section):
    """

    @param config_:
    @param section:

    @return dict1:
    """
    dict1 = {}
    options = config_.options(section)
    for option in options:
        try:
            dict1[option] = config_.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except KeyError:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def extract_settings_elvis():
    """ creates a dictionary with all the configuration parameters
        at this moment configuration file location is fixed at main directory

    @return prfs_d: a dictionary which contains all valuable data
    """
    cf = ConfigParser()
    cf.read(".settings_ELViS.ini")

    prfs_d = {}
    os_version = get_os()

    if os_version == 'centos':
        prfs_d['version'] = confmap(cf, "Version")['centos_version']
    elif os_version == 'cab':
        prfs_d['version'] = confmap(cf, "Version")['cab_version']
    else:
        raise BadSettings('Operative system not chosen')

    if os_version == 'centos':
        prfs_d['home'] = confmap(cf, "HomeDirs")['centos_home']
    elif os_version == 'cab':
        prfs_d['home'] = confmap(cf, "HomeDirs")['cab_home']
    else:
        raise BadSettings('Operative system not chosen')

    prfs_d['fits_dir'] = confmap(cf, "ImagesDirs")['fits_dir']
    prfs_d['fits_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fits_dir'])
    prfs_d['fpas_dir'] = confmap(cf, "ImagesDirs")['fpas_dir']
    prfs_d['fpas_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fpas_dir'])

    # todo - comment!
    prfs_d['output_cats'] = confmap(cf, "CatsDirs")['output_cats']
    prfs_d['output_cats'] = prfs_d['version'] + prfs_d['output_cats']
    # todo - comment!
    prfs_d['references'] = confmap(cf, "CatsDirs")['references']
    prfs_d['references'] = prfs_d['version'] + prfs_d['references']
    # todo - comment!
    prfs_d['filtered'] = confmap(cf, "CatsDirs")['filtered']
    prfs_d['filtered'] = prfs_d['version'] + prfs_d['filtered']

    prfs_d['time_1'] = confmap(cf, "ImagesTime")['time_1']  # 1st dither time
    prfs_d['time_2'] = confmap(cf, "ImagesTime")['time_2']  # 2nd dither time
    prfs_d['time_3'] = confmap(cf, "ImagesTime")['time_3']  # 3nd dither time
    prfs_d['time_4'] = confmap(cf, "ImagesTime")['time_4']  # 4th dither time

    outputdirs_list = ['conf_scamp', 'conf_sex', 'params_sex', 'neural_sex',
                       'params_cat', 'logger_config']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = confmap(cf, "ConfigDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['detections'] = int(confmap(cf, "Misc")['detections'])
    prfs_d['pm_low'] = float(confmap(cf, "Misc")['pm_low'])
    prfs_d['pm_up'] = float(confmap(cf, "Misc")['pm_up'])
    prfs_d['pm_sn'] = float(confmap(cf, "Misc")['pm_sn'])
    prfs_d['r_fit'] = confmap(cf, "Misc")['r_fit']
    prfs_d['cores_number'] = confmap(cf, "Misc")['cores_number']
    if prfs_d['cores_number'] == '0':
        prfs_d['cores_number'] = int(str(cpu_count()))
        # TODO should leave free at least 20% of processors
    else:
        prfs_d['cores_number'] = int(prfs_d['cores_number'])
    prfs_d['tolerance'] = float(confmap(cf, "Misc")['tolerance'])

    return prfs_d


def extract_settings():
    """ creates a dictionary with all the configuration parameters
        at this moment configuration file location is fixed at main directory

    @return prfs_d: a dictionary which contains all valuable data
    """
    cf = ConfigParser()
    cf.read(".settings.ini")

    os_version = get_os()

    prfs_d = {'cat': confmap(cf, "Version")['cat_version']}

    if os_version == 'fedora':
        prfs_d['home'] = confmap(cf, "HomeDirs")['fed_home']
    elif os_version == 'ubuntu':
        prfs_d['home'] = confmap(cf, "HomeDirs")['ub_home']
    elif os_version == 'test':
        prfs_d['home'] = confmap(cf, "HomeDirs")['test_home']
    elif os_version == 'centos':
        prfs_d['home'] = confmap(cf, "HomeDirs")['centos_home']
    else:
        raise BadSettings('Operative system not chosen')

    if os_version == 'fedora':
        prfs_d['version'] = confmap(cf, "Version")['fed_version']
        prfs_d['version'] = prfs_d['version'] + prfs_d['cat']
    elif os_version == 'ubuntu':
        prfs_d['version'] = confmap(cf, "Version")['ub_version']
        prfs_d['version'] = prfs_d['version'] + prfs_d['cat']
    elif os_version == 'test':
        prfs_d['version'] = confmap(cf, "Version")['test_version']
        prfs_d['version'] = prfs_d['version'] + prfs_d['cat']
    elif os_version == 'centos':
        prfs_d['version'] = confmap(cf, "Version")['centos_version']
        prfs_d['version'] = prfs_d['version'] + prfs_d['cat']
    else:
        raise BadSettings('Operative system not chosen')

    prfs_d['fits_dir'] = confmap(cf, "ImagesDirs")['fits_dir']
    # prfs_d['fits_dir'] = prfs_d['version'] + prfs_d['fits_dir']
    prfs_d['fits_dir'] = prfs_d['version']  # hardcoded

    prfs_d['fpas_dir'] = confmap(cf, "ImagesDirs")['fpas_dir']
    prfs_d['fpas_dir'] = prfs_d['version'] + prfs_d['fpas_dir']
    # TODO This hardcoded lines should be improve
    prfs_d['fits_ref'] = confmap(cf, "ImagesDirs")['fits_ref']

    prfs_d['time_1'] = confmap(cf, "ImagesTime")['time_1']
    prfs_d['time_2'] = confmap(cf, "ImagesTime")['time_2']
    prfs_d['time_3'] = confmap(cf, "ImagesTime")['time_3']
    prfs_d['time_4'] = confmap(cf, "ImagesTime")['time_4']

    outputdirs_list = ['conf_scamp', 'conf_sex', 'params_sex', 'neural_sex',
                       'params_cat', 'logger_config']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = confmap(cf, "ConfigDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['output_cats'] = confmap(cf, "CatsDirs")['output_cats']
    prfs_d['output_cats'] = prfs_d['version'] + prfs_d['output_cats']
    prfs_d['input_cats'] = confmap(cf, "CatsDirs")['input_cats']
    prfs_d['input_cats'] = prfs_d['version'] + prfs_d['input_cats']
    prfs_d['input_ref'] = confmap(cf, "CatsDirs")['input_ref']

    prfs_d['first_star'] = confmap(cf, "CatsOrganization")['first_star']
    prfs_d['first_star'] = int(prfs_d['first_star'])
    prfs_d['first_galaxy'] = confmap(cf, "CatsOrganization")['first_galaxy']
    prfs_d['first_galaxy'] = int(prfs_d['first_galaxy'])
    prfs_d['first_sso'] = confmap(cf, "CatsOrganization")['first_sso']
    prfs_d['first_sso'] = int(prfs_d['first_sso'])

    outputdirs_list = ['plots_dir', 'results_dir', 'images_out', 'fits_out',
                       'report_out', 'dithers_out', 'catalogs_dir', 'tmp_out',
                       'filter_dir']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = confmap(cf, "OutputDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['detections'] = int(confmap(cf, "Misc")['detections'])
    prfs_d['pm_low'] = float(confmap(cf, "Misc")['pm_low'])
    prfs_d['pm_up'] = float(confmap(cf, "Misc")['pm_up'])
    prfs_d['pm_sn'] = float(confmap(cf, "Misc")['pm_sn'])
    pms = confmap(cf, "Misc")['pms']
    pms = pms.replace(",", " ")
    prfs_d['pms'] = [float(x) for x in pms.split()]
    mags = confmap(cf, "Misc")['mags']
    mags = mags.replace(",", " ")
    prfs_d['mags'] = mags.split()
    confidences = confmap(cf, "Misc")['confidences']
    confidences = confidences.replace(",", " ")
    prfs_d['confidences'] = [int(x) for x in confidences.split()]
    cross_ids = confmap(cf, "Misc")['cross_ids']
    cross_ids = cross_ids.replace(",", " ")
    prfs_d['cross_ids'] = [float(x) for x in cross_ids.split()]
    prfs_d['r_fit'] = confmap(cf, "Misc")['r_fit']
    prfs_d['cores_number'] = confmap(cf, "Misc")['cores_number']
    if prfs_d['cores_number'] == '0':
        prfs_d['cores_number'] = int(str(cpu_count()))
        # TODO should leave free at least 20% of processors
    else:
        prfs_d['cores_number'] = int(prfs_d['cores_number'])
    prfs_d['tolerance'] = float(confmap(cf, "Misc")['tolerance'])

    return prfs_d


def pipeline_help(logger):
    """ returns a helpful log information
    To-do Need to be improved

    @param logger: a logger object.

    @return True: if everything goes alright
    """
    logger.info('An invalid option has been chosen')
    logger.info('Select:')
    logger.info('- reports')
    logger.info('    - runs the full pipeline creating a report file')
    logger.info('- stats')
    logger.info('    - runs the full pipeline without reports analysing data')
    logger.info('- analysis')
    logger.info('    - only runs the analysis pipeline')
    logger.info('- scamp')
    logger.info('    - only runs the scamp pipeline')
    logger.info('- catalogue')
    logger.info('    - only runs the catalogue creation script')
    logger.info('- help')
    logger.info('    - shows an helpful message')

    return True


def get_ticks(min_number, max_number, resolution):
    """ given a number returns a list with all its submultiples

    @param min_number:
    @param max_number:
    @param resolution:

    @return ticks:
    """
    number_ticks = int(max_number / resolution + 1)
    ticks = []
    for i in range(int(min_number), number_ticks, 1):
        ticks.append(i * resolution)

    return ticks


def get_limits(first_list, second_list):
    """ given two lists returns the minimum and the maximum value of them

    @param first_list: first list to be analysed
    @param second_list: second list to be analysed

    @return limits: a list of two values the minimun and the maximun one
    """

    if max(first_list) > max(second_list):
        high = max(first_list)
    else:
        high = max(second_list)

    if min(first_list) > min(second_list):
        low = min(second_list)
    else:
        low = min(first_list)

    limits = [low, high]

    return limits


def pm_compute(logger, merged_df, full_df):
    """ given a merged catalogue and a full one return a full catalogue with
        proper motions in arcseconds per hour

    @param logger: a logger object
    @param merged_df:
    @param full_df:

    @return db: a dataframe with all proper motions values
    """
    logger.debug('Computing right ascension proper motion')
    pmalpha = merged_df['PMALPHA_J2000'].divide(8.75e6)
    pmalpha_l = []
    logger.debug('Computing declination proper motion')
    pmdelta = merged_df['PMDELTA_J2000'].divide(8.75e6)
    pmdelta_l = []
    logger.debug('Computing right ascension proper motion error')
    pmealpha = merged_df['PMALPHAERR_J2000'].divide(8.75e6)
    pmealpha_l = []
    logger.debug('Computing declination proper motion error')
    pmedelta = merged_df['PMDELTAERR_J2000'].divide(8.75e6)
    pmedelta_l = []

    logger.debug('Computing proper motion')
    pm = Series(np.sqrt(np.array(pmalpha**2 + pmdelta**2), dtype=float))
    pm_l = []
    logger.debug('Computing proper motion error')
    pme = Series(np.sqrt(np.array(pmealpha**2 + pmedelta**2), dtype=float))
    pme_l = []

    for idx_merged, source in enumerate(merged_df['SOURCE_NUMBER']):
        # print('{} - {}'.format(idx_merged, len(merged_db['SOURCE_NUMBER'])))
        full_p_db = full_df[full_df['SOURCE_NUMBER'].isin([source])]

        for idx in full_p_db['SOURCE_NUMBER']:
            pmalpha_l.append(pmalpha.iloc[idx_merged])
            pmdelta_l.append(pmdelta.iloc[idx_merged])
            pmealpha_l.append(pmealpha.iloc[idx_merged])
            pmedelta_l.append(pmedelta.iloc[idx_merged])
            if pm.iloc[idx_merged] < 0:
                print(pm.iloc[idx_merged])
            pm_l.append(pm.iloc[idx_merged])
            pme_l.append(pme.iloc[idx_merged])

    # Series creation
    pmalpha_s = Series(pmalpha_l, name='PMALPHA', dtype=float)
    pmdelta_s = Series(pmdelta_l, name='PMDELTA', dtype=float)
    pmealpha_s = Series(pmealpha_l, name='PMALPHAERR', dtype=float)
    pmedelta_s = Series(pmedelta_l, name='PMDELTAERR', dtype=float)

    pm_s = Series(pm_l, name='PM', dtype=float)
    pme_s = Series(pme_l, name='PMERR', dtype=float)

    df = concat([full_df['SOURCE_NUMBER'].reset_index(),
                 full_df['CATALOG_NUMBER'].reset_index(),
                 full_df['EXTENSION'].reset_index(),
                 full_df['ASTR_INSTRUM'].reset_index(),
                 full_df['PHOT_INSTRUM'].reset_index(),
                 full_df['X_IMAGE'].reset_index(),
                 full_df['Y_IMAGE'].reset_index(),
                 full_df['ERRA_IMAGE'].reset_index(),
                 full_df['ERRB_IMAGE'].reset_index(),
                 full_df['ERRTHETA_IMAGE'].reset_index(),
                 full_df['ALPHA_J2000'].reset_index(),
                 full_df['DELTA_J2000'].reset_index(),
                 full_df['ERRA_WORLD'].reset_index(),
                 full_df['ERRB_WORLD'].reset_index(),
                 full_df['ERRTHETA_WORLD'].reset_index(),
                 full_df['EPOCH'].reset_index(), full_df['MAG'].reset_index(),
                 full_df['MAGERR'].reset_index(),
                 full_df['FLAGS_EXTRACTION'].reset_index(),
                 full_df['FLAGS_SCAMP'].reset_index(),
                 full_df['FLAGS_IMA'].reset_index(), pmalpha_s,
                 pmdelta_s, pmealpha_s, pmedelta_s, pm_s, pme_s], axis=1)
    del df['index']

    return df


def pm_filter(full_db, pm_low, pm_up):
    """ filters proper motions values according to their values

    :param full_db:
    :param pm_low:
    :param pm_up:
    :return: full_db
    """
    mask = (full_db['PM'] > float(pm_low)) & (full_db['PM'] < float(pm_up))
    full_db = full_db[mask]

    return full_db


def b_filter(full_db, b_low, b_up):
    """ filters proper motions values according to their values

    :param full_db:
    :param b_low:
    :param b_up:
    :return: full_db
    """
    up = full_db['MEAN_B_IMAGE'] > float(b_low)
    down = full_db['MEAN_B_IMAGE'] < float(b_up)
    mask = up & down
    full_db = full_db[mask]

    return full_db


def confidence_filter(db, r):
    """ filter objects with a non-coherence movement

    :param db:
    :param r:
    :return:
    """
    passed = []

    # Get Coordintates, Errors, and Epochs
    for i, source in enumerate(set(db['SOURCE_NUMBER'])):
        ra = db.loc[db['SOURCE_NUMBER'] == source, 'ALPHA_J2000'].tolist()
        dec = db.loc[db['SOURCE_NUMBER'] == source, 'DELTA_J2000'].tolist()
        epoch = db.loc[db['SOURCE_NUMBER'] == source, 'EPOCH'].tolist()
        # cats = cats_coherence(db.loc[db['SOURCE_NUMBER'] == source,
        #                       'CATALOG_NUMBER'].tolist())
        # Calculate WLS fit and confidence interval
        for dimension in [ra, dec]:
            x = np.array(epoch)
            y = np.array(dimension)
            sigma = []
            if dimension == ra:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRA_WORLD'].tolist()
            if dimension == dec:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRB_WORLD'].tolist()

            edim = np.array([1 / var for var in sigma])
            x = sm.add_constant(x)
            # Model: y~x+c
            model = sm.WLS(y, x, weigths=edim)
            fitted = model.fit()
            if fitted.rsquared >= float(r):
                passed.append(source)

    # Passed if both dimension have required rsquared
    passed = [p for p in passed if passed.count(p) >= 1]
    passed = list(set(passed))
    db = db[db['SOURCE_NUMBER'].isin(passed)]

    return db


def setting_logger(prfs_d):
    """ sets-up a logger object ready to be used

    TODO improve logger definition

    @return logger:
    """
    print(prfs_d['logger_config'])
    config.fileConfig(prfs_d['logger_config'])

    # TODO implement logger level setting
    """
    if argv[4] == '-INFO':
        logger = getLogger("main_process").setLevel(INFO)
    elif argv[4] == '-DEBUG':
        logger = getLogger("main_process").setLevel(DEBUG)
    else:
        raise Exception
    """
    logger = getLogger("main_process")
    logger.info("Pipeline started")
    # logger.FileHandler('spam.log')

    return logger


def clear_data(logger, prfs_d, catalogues):
    """ Clear old fits and jpeg images. Removes catalogue files too.

    TODO it is still useful?

    @param logger: a logger object
    @param prfs_d:
    @param catalogues: if True removes catalogues too

    @return True: if everything goes alright
    """

    if catalogues:
        logger.debug('setting variables to be remove')
        pipeline_dir = prfs_d['home'] + '/pipeline/'
        swarp_xml = pipeline_dir + 'swarp.xml'
        merged_fits = pipeline_dir + 'coadd.fits'
        merged_weight_fits = pipeline_dir + 'coadd.weight.fits'

        try:
            logger.debug("removing old xml swarp's file")
            os.remove(swarp_xml)
        except OSError:
            pass

        try:
            logger.debug('removes merged fits {}'.format(merged_fits))
            os.remove(merged_fits)
            logger.debug('removes merged weight {}'.format(merged_weight_fits))
            os.remove(merged_weight_fits)
        except OSError:
            logger.error("merged fits cannot be removed")
            pass

        try:
            cats_files = os.listdir(prfs_d['fits_dir'])
            for cat in cats_files:
                if cat[:3] == 'mag' and cat[-4:] == '.cat':
                    cat_name = prfs_d['output_cats'] + '/' + cat
                    logger.debug('removing...{}'.format(cat_name))
                    os.remove(prfs_d['output_cats'] + '/' + cat)
        except OSError:
            logger.error("old sextracted catalogues cannot be removed")
            pass

        try:
            for region_idx in range(0, len(prfs_d['mags']), 1):
                logger.debug('removes regions {}')
                mags = prfs_d['mags']
                mag = mags[region_idx]
                home = prfs_d['home']
                os.remove(home + '/regions_final_m{}.reg'.format(mag))
        except OSError:
            logger.error('final regions cannot be removed')

    return True


def compare_floats(f1, f2, allowed_error):
    """ given a pair of floats and a tolerance return true if
        the difference between them is littler than tolerance

        although there is a native method at numpy library
        this approach seems to be easier

    @param f1: a float to be compared
    @param f2: a float to be compared
    @param allowed_error: tolerance admissible between floats

    @return True: if floats are close enought
    """
    return abs(f1 - f2) <= allowed_error


def check_distance(x1, x2, y1, y2, allowed_error):
    """

    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param allowed_error:
    :return:
    """
    distance = hypot(x2 - x1, y2 - y1)

    if distance <= allowed_error:
        return True, distance
    else:
        return False, distance


def speeds_range(prfs_d, confidence):
    """ given a confidence value returns a dict with speeds

    @param prfs_d:
    @param confidence:

    @return speeds_dict:
    """
    speeds_dict = {}
    for pm_ in prfs_d['pms']:
        speeds_dict[pm_] = [pm_ - pm_ * confidence / 100,
                            pm_ + pm_ * confidence / 100]

    return speeds_dict


def measure_distance(x2, x1, y2, y1):
    """

    :param x2:
    :param x1:
    :param y2:
    :param y1:
    :return:
    """

    distance = np.sqrt((x2 - x1) * (x2 - x1) - (y2 - y1) * (y2 - y1))

    return distance


def significant_l(number):
    """

    :param number:
    :return:
    """
    len_ = len(str(number))
    a = ('%.' + str(len_) + 'E') % Decimal(number)
    significate_d = a.split(".")[0]
    times = a.split("E")[1]
    result = int(significate_d) * (10 ** int(times))

    return result


def create_folder(logger, folder):
    """ creates a folder if it is necessary

    :param logger:
    :param folder:
    :return:
    """

    try:
        if not os.path.isfile(folder):
            os.makedirs(folder)
        return True
    except OSError:
        logger.debug('Folder {} already created'.format(folder))
        return True
    except Exception as e:
        print('Unexpected error: {}'.format(e))


def get_dither(catalog_n):
    """

    :param catalog_n:
    :return: dither_n
    """
    cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
            ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
            ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
            ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
            ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
            ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
            ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
            ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
            ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
            ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
            ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
            ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

    ccd = ''
    dither = 0
    for cat_ in cats:
        if cat_[2] == catalog_n:
            ccd = cat_[0]
            dither = cat_[1]

    return ccd, int(dither)


def get_norm_speed(o_pm):
    """

    :param o_pm:
    :return: pm_norm
    """
    prfs_d = extract_settings()
    speeds_d = speeds_range(prfs_d, 50)

    pm_norm = 0
    for key_ in speeds_d.keys():
        low = speeds_d[key_][0]
        high = speeds_d[key_][1]
        if low < o_pm < high:
            pm_norm = key_

    return pm_norm


def check_source(o_df, o_alpha, o_delta):
    """

    :param o_df:
    :param o_alpha:
    :param o_delta:
    :return:
    """
    prfs_d = extract_settings()

    o_df = o_df[o_df['ALPHA_J2000'] + prfs_d['tolerance'] > o_alpha]
    o_df = o_df[o_alpha > o_df['ALPHA_J2000'] - prfs_d['tolerance']]
    o_df = o_df[o_df['DELTA_J2000'] + prfs_d['tolerance'] > o_delta]
    o_df = o_df[o_delta > o_df['DELTA_J2000'] - prfs_d['tolerance']]

    return o_df


def check_source_elvis(o_df, o_alpha, o_delta):
    """

    :param o_df:
    :param o_alpha:
    :param o_delta:
    :return:
    """
    prfs_d = extract_settings_elvis()

    o_df = o_df[o_df['ALPHA_J2000'] + prfs_d['tolerance'] > o_alpha]
    o_df = o_df[o_alpha > o_df['ALPHA_J2000'] - prfs_d['tolerance']]
    o_df = o_df[o_df['DELTA_J2000'] + prfs_d['tolerance'] > o_delta]
    o_df = o_df[o_delta > o_df['DELTA_J2000'] - prfs_d['tolerance']]

    return o_df
