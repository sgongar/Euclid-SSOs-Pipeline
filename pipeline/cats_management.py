#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

Todo:
    * Improve log messages

"""

from multiprocessing import Process
from os import listdir, remove

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from images_management import get_fits_limits
import itertools
from misc import check_distance, get_fits, extract_settings
from numpy import genfromtxt, mean, float64
from numpy import asarray, sqrt, isclose, logical_and
from regions import Create_regions
from pandas import concat, DataFrame, Series, read_csv

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def rebase_catalogue(logger, mag, prfs_d, remove_ssos, cat_location):
    """

    @param logger: a logger object
    @param prfs_d:
    @param remove_ssos:
    @param cat_location:

    @return True: if everything goes alright
    """
    # salida del catalogo por sextractor
    logger.debug('opening sextracted catalog from {}'.format(cat_location))
    sex_table = get_output_catalog(prfs_d, cat_location)

    # entrada del catalogo por luca
    cat_number = 10
    complete_list = False
    input_cat = prfs_d['input_cats'] + '/Cat_20-21_d1.dat'
    logger.debug('opening input catalog from {}'.format(input_cat))
    input_table = get_input_cat(logger, prfs_d, cat_number,
                                input_cat, complete_list)
    # salida de los dos
    fits_loc = '{}/m_{}_x0_y0_d1.fits'.format(prfs_d['fits_dir'], mag)
    hdulist = fits.open(fits_loc)
    logger.debug('getting WCS information from {}'.format(fits_loc))
    w = WCS(hdulist[0].header)

    logger.debug('getting list of x/y values')
    regions_list = []
    for source_num in range(input_table['x_values'].as_matrix().size):
        x_value = input_table['x_values'].as_matrix()[source_num]
        y_value = input_table['y_values'].as_matrix()[source_num]
        regions_list.append([x_value, y_value])

    logger.debug('transforming x/y values to ra/dec ones')
    input_regions = w.wcs_pix2world(regions_list, 1)

    logger.debug('creating new list with ra/dec values')
    alpha_list = []
    delta_list = []
    for regions in input_regions:
        alpha_list.append(regions[0])
        delta_list.append(regions[1])

    logger.debug('creating new Pandas series for ra/dec values')
    alpha_j2000 = Series(alpha_list, name='alpha_j2000')
    delta_j2000 = Series(delta_list, name='delta_j2000')

    logger.debug('concatenating series of ra and dec')
    coords_table = concat([alpha_j2000, delta_j2000], axis=1)

    """ TODO implement a improve comparation algorithm
    av_r_motion = (pos_3_ra - pos_1_ra)
    av_d_motion = (pos_3_dec - pos_1_dec)

    av_motion = hypot(abs(av_r_motion), abs(av_d_motion))
    """

    list_sex_alpha_image = []
    list_sex_delta_image = []
    list_sex_number = []
    list_sex_x_image = []
    list_sex_y_image = []
    list_sex_x2_image = []
    list_sex_y2_image = []
    list_sex_xy_image = []
    list_sex_isoarea_image = []
    list_sex_background = []
    list_sex_threshold = []
    list_sex_flux_max = []
    list_sex_a_image = []
    list_sex_b_image = []
    list_sex_theta_image = []
    list_sex_erra_image = []
    list_sex_errb_image = []
    list_sex_flux_iso = []
    list_sex_fluxerr_iso = []
    list_sex_mag_iso = []
    list_sex_magerr_iso = []
    list_sex_flux_aper = []
    list_sex_fluxerr_aper = []
    list_sex_mag_aper = []
    list_sex_magerr_aper = []
    list_sex_alpha_sky = []
    list_sex_delta_sky = []
    list_sex_b_image = []
    list_sex_theta_image = []
    list_sex_errtheta_image = []
    list_sex_mu_max = []
    list_sex_fwhm_image = []
    list_sex_flux_radius = []
    list_sex_elongation = []
    list_sex_ellipticity = []
    list_sex_cxx_image = []
    list_sex_cxy_image = []
    list_sex_cyy_image = []
    list_sex_errcxx_image = []
    list_sex_errcxy_image = []
    list_sex_errcyy_image = []
    list_sex_mag_auto = []
    list_sex_xwin_image = []
    list_sex_ywin_image = []
    list_sex_flux_auto = []
    list_sex_fluxerr_auto = []
    list_sex_magerr_auto = []
    list_sex_snr_win = []
    list_sex_x_world = []
    list_sex_y_world = []
    list_sex_errx2_world = []
    list_sex_erry2_world = []
    list_sex_erryxy_world = []
    list_sex_awin_image = []
    list_sex_bwin_image = []
    list_sex_thetawin_image = []
    list_sex_errawin_image = []
    list_sex_errbwin_image = []
    list_sex_errthetawin_image = []
    list_sex_flags = []
    list_sex_fwhm_world = []
    list_sex_erra_world = []
    list_sex_errb_world = []

    sources_number = alpha_j2000.size  # same as delta_j2000.size
    logger.debug('{} sources to be analysed'.format(sources_number))
    for source_num in range(sources_number):
        # voy mirando una a una las fuentes de sextractor y si esta en el
        # original copio todo a una table nueva
        # getting alpha, delta for each source sextracted
        # FIXME no uso valores absolutos!
        source_alpha = coords_table['alpha_j2000'].iloc[source_num]
        source_delta = coords_table['delta_j2000'].iloc[source_num]

        alpha_mask = isclose(sex_table['ALPHA_J2000'].values[:, None],
                             [source_alpha],
                             atol=prfs_d['alpha_t']).any(axis=1)
        delta_mask = isclose(sex_table['DELTA_J2000'].values[:, None],
                             [source_delta],
                             atol=prfs_d['delta_t']).any(axis=1)

        mask = logical_and(alpha_mask, delta_mask)
        df = sex_table[mask]
        if not df.empty:
            for sex_n in range(df['ALPHA_J2000'].size):
                sex_alpha_image = df['ALPHA_J2000'].iloc[sex_n]
                list_sex_alpha_image.append(sex_alpha_image)
                sex_delta_image = df['DELTA_J2000'].iloc[sex_n]
                list_sex_delta_image.append(sex_delta_image)
                sex_number = df['NUMBER'].iloc[sex_n]
                list_sex_number.append(sex_number)
                sex_x_image = df['X_IMAGE'].iloc[sex_n]
                list_sex_x_image.append(sex_x_image)
                sex_y_image = df['Y_IMAGE'].iloc[sex_n]
                list_sex_y_image.append(sex_y_image)
                sex_x2_image = df['X2_IMAGE'].iloc[sex_n]
                list_sex_x2_image.append(sex_x2_image)
                sex_y2_image = df['Y2_IMAGE'].iloc[sex_n]
                list_sex_y2_image.append(sex_y2_image)
                sex_xy_image = df['XY_IMAGE'].iloc[sex_n]
                list_sex_xy_image.append(sex_xy_image)
                sex_isoarea_image = df['ISOAREA_IMAGE'].iloc[sex_n]
                list_sex_isoarea_image.append(sex_isoarea_image)
                sex_background = df['BACKGROUND'].iloc[sex_n]
                list_sex_background.append(sex_background)
                sex_threshold = df['THRESHOLD'].iloc[sex_n]
                list_sex_threshold.append(sex_threshold)
                sex_flux_max = df['FLUX_MAX'].iloc[sex_n]
                list_sex_flux_max.append(sex_flux_max)
                sex_a_image = df['A_IMAGE'].iloc[sex_n]
                list_sex_a_image.append(sex_a_image)
                sex_b_image = df['B_IMAGE'].iloc[sex_n]
                list_sex_b_image.append(sex_b_image)
                sex_theta_image = df['THETA_IMAGE'].iloc[sex_n]
                list_sex_theta_image.append(sex_theta_image)
                sex_erra_image = df['ERRA_IMAGE'].iloc[sex_n]
                list_sex_erra_image.append(sex_erra_image)
                sex_errb_image = df['ERRB_IMAGE'].iloc[sex_n]
                list_sex_errb_image.append(sex_errb_image)
                sex_flux_iso = df['FLUX_ISO'].iloc[sex_n]
                list_sex_flux_iso.append(sex_flux_iso)
                sex_fluxerr_iso = df['FLUXERR_ISO'].iloc[sex_n]
                list_sex_fluxerr_iso.append(sex_fluxerr_iso)
                sex_mag_iso = df['MAG_ISO'].iloc[sex_n]
                list_sex_mag_iso.append(sex_mag_iso)
                sex_magerr_iso = df['MAGERR_ISO'].iloc[sex_n]
                list_sex_magerr_iso.append(sex_magerr_iso)
                sex_flux_aper = df['FLUX_APER'].iloc[sex_n]
                list_sex_flux_aper.append(sex_flux_aper)
                sex_fluxerr_aper = df['FLUXERR_APER'].iloc[sex_n]
                list_sex_fluxerr_aper.append(sex_fluxerr_aper)
                sex_mag_aper = df['MAG_APER'].iloc[sex_n]
                list_sex_mag_aper.append(sex_mag_aper)
                sex_magerr_aper = df['MAGERR_APER'].iloc[sex_n]
                list_sex_magerr_aper.append(sex_magerr_aper)
                sex_alpha_sky = df['ALPHA_SKY'].iloc[sex_n]
                list_sex_alpha_sky.append(sex_alpha_sky)
                sex_delta_sky = df['DELTA_SKY'].iloc[sex_n]
                list_sex_delta_sky.append(sex_delta_sky)
                sex_errtheta_image = df['ERRTHETA_IMAGE'].iloc[sex_n]
                list_sex_errtheta_image.append(sex_errtheta_image)
                sex_mu_max = df['MU_MAX'].iloc[sex_n]
                list_sex_mu_max.append(sex_mu_max)
                sex_fwhm_image = df['FWHM_IMAGE'].iloc[sex_n]
                list_sex_fwhm_image.append(sex_fwhm_image)
                sex_flux_radius = df['FLUX_RADIUS'].iloc[sex_n]
                list_sex_flux_radius.append(sex_flux_radius)
                sex_elongation = df['ELONGATION'].iloc[sex_n]
                list_sex_elongation.append(sex_elongation)
                sex_ellipticity = df['ELLIPTICITY'].iloc[sex_n]
                list_sex_ellipticity.append(sex_ellipticity)
                sex_cxx_image = df['CXX_IMAGE'].iloc[sex_n]
                list_sex_cxx_image.append(sex_cxx_image)
                sex_cxy_image = df['CXY_IMAGE'].iloc[sex_n]
                list_sex_cxy_image.append(sex_cxy_image)
                sex_cyy_image = df['CYY_IMAGE'].iloc[sex_n]
                list_sex_cyy_image.append(sex_cyy_image)
                sex_errcxx_image = df['ERRCXX_IMAGE'].iloc[sex_n]
                list_sex_errcxx_image.append(sex_errcxx_image)
                sex_errcxy_image = df['ERRCXY_IMAGE'].iloc[sex_n]
                list_sex_errcxy_image.append(sex_errcxy_image)
                sex_errcyy_image = df['ERRCYY_IMAGE'].iloc[sex_n]
                list_sex_errcyy_image.append(sex_errcyy_image)
                sex_mag_auto = df['MAG_AUTO'].iloc[sex_n]
                list_sex_mag_auto.append(sex_mag_auto)
                sex_xwin_image = df['XWIN_IMAGE'].iloc[sex_n]
                list_sex_xwin_image.append(sex_xwin_image)
                sex_ywin_image = df['YWIN_IMAGE'].iloc[sex_n]
                list_sex_ywin_image.append(sex_ywin_image)
                sex_flux_auto = df['FLUX_AUTO'].iloc[sex_n]
                list_sex_flux_auto.append(sex_flux_auto)
                sex_fluxerr_auto = df['FLUXERR_AUTO'].iloc[sex_n]
                list_sex_fluxerr_auto.append(sex_fluxerr_auto)
                sex_magerr_auto = df['MAGERR_AUTO'].iloc[sex_n]
                list_sex_magerr_auto.append(sex_magerr_auto)
                sex_snr_win = df['SNR_WIN'].iloc[sex_n]
                list_sex_snr_win.append(sex_snr_win)
                sex_x_world = df['X_WORLD'].iloc[sex_n]
                list_sex_x_world.append(sex_x_world)
                sex_y_world = df['Y_WORLD'].iloc[sex_n]
                list_sex_y_world.append(sex_y_world)
                sex_errx2_world = df['ERRX2_WORLD'].iloc[sex_n]
                list_sex_errx2_world.append(sex_errx2_world)
                sex_erry2_world = df['ERRY2_WORLD'].iloc[sex_n]
                list_sex_erry2_world.append(sex_erry2_world)
                sex_erryxy_world = df['ERRXY_WORLD'].iloc[sex_n]
                list_sex_erryxy_world.append(sex_erryxy_world)
                sex_awin_image = df['AWIN_IMAGE'].iloc[sex_n]
                list_sex_awin_image.append(sex_awin_image)
                sex_bwin_image = df['BWIN_IMAGE'].iloc[sex_n]
                list_sex_bwin_image.append(sex_bwin_image)
                sex_thetawin_image = df['THETAWIN_IMAGE'].iloc[sex_n]
                list_sex_thetawin_image.append(sex_thetawin_image)
                sex_errawin_image = df['ERRAWIN_IMAGE'].iloc[sex_n]
                list_sex_errawin_image.append(sex_errawin_image)
                sex_errbwin_image = df['ERRBWIN_IMAGE'].iloc[sex_n]
                list_sex_errbwin_image.append(sex_errbwin_image)
                sex_errthetawin_image = df['ERRTHETAWIN_IMAGE'].iloc[sex_n]
                list_sex_errthetawin_image.append(sex_errthetawin_image)
                sex_flags = df['FLAGS'].iloc[sex_n]
                list_sex_flags.append(sex_flags)
                sex_fwhm_world = df['FWHM_WORLD'].iloc[sex_n]
                list_sex_fwhm_world.append(sex_fwhm_world)
                sex_erra_world = df['ERRA_WORLD'].iloc[sex_n]
                list_sex_erra_world.append(sex_erra_world)
                sex_errb_world = df['ERRB_WORLD'].iloc[sex_n]
                list_sex_errb_world.append(sex_errb_world)

    ALPHA_J2000_final = Series(list_sex_alpha_image, name='ALPHA_J2000')
    DELTA_J2000_final = Series(list_sex_delta_image, name='DELTA_J2000')
    NUMBER_final = Series(list_sex_number, name='NUMBER')
    X_IMAGE_final = Series(list_sex_x_image, name='X_IMAGE')
    Y_IMAGE_final = Series(list_sex_y_image, name='Y_IMAGE')
    X2_IMAGE_final = Series(list_sex_x2_image, name='X2_IMAGE')
    Y2_IMAGE_final = Series(list_sex_y2_image, name='Y2_IMAGE')
    XY_IMAGE_final = Series(list_sex_xy_image, name='XY_IMAGE')
    ISOAREA_IMAGE_final = Series(list_sex_isoarea_image, name='ISOAREA_IMAGE')
    BACKGROUND_final = Series(list_sex_background, name='BACKGROUND')
    THRESHOLD_final = Series(list_sex_threshold, name='THRESHOLD')
    FLUX_MAX_final = Series(list_sex_flux_max, name='FLUX_MAX')
    A_IMAGE_final = Series(list_sex_a_image, name='A_IMAGE')
    B_IMAGE_final = Series(list_sex_b_image, name='B_IMAGE')
    THETA_IMAGE_final = Series(list_sex_theta_image, name='THETA_IMAGE')
    ERRA_IMAGE_final = Series(list_sex_erra_image, name='ERRA_IMAGE')
    ERRB_IMAGE_final = Series(list_sex_errb_image, name='ERRB_IMAGE')
    FLUX_ISO_final = Series(list_sex_flux_iso, name='FLUX_ISO')
    FLUXERR_ISO_final = Series(list_sex_fluxerr_iso, name='FLUXERR_ISO')
    MAG_ISO_final = Series(list_sex_mag_iso, name='MAG_ISO')
    MAGERR_ISO_final = Series(list_sex_magerr_iso, name='MAGERR_ISO')
    FLUX_APER_final = Series(list_sex_flux_aper, name='FLUX_APER')
    FLUXERR_APER_final = Series(list_sex_fluxerr_aper, name='FLUXERR_APER')
    MAG_APER_final = Series(list_sex_mag_aper, name='MAG_APER')
    MAGERR_APER_final = Series(list_sex_magerr_aper, name='MAGERR_APER')
    ALPHA_SKY_final = Series(list_sex_alpha_sky, name='ALPHA_SKY')
    DELTA_SKY_final = Series(list_sex_delta_sky, name='DELTA_SKY')
    ERRTHETA_IMAGE_final = Series(list_sex_errtheta_image,
                                  name='ERRTHETA_IMAGE')
    MU_MAX_final = Series(list_sex_mu_max, name='MU_MAX')
    FWHM_IMAGE_final = Series(list_sex_fwhm_image, name='FWHM_IMAGE')
    FLUX_RADIUS_final = Series(list_sex_flux_radius, name='FLUX_RADIUS')
    ELONGATION_final = Series(list_sex_elongation, name='ELONGATION')
    ELLIPTICITY_final = Series(list_sex_ellipticity, name='ELLIPTICITY')
    CXX_IMAGE_final = Series(list_sex_cxx_image, name='CXX_IMAGE')
    CXY_IMAGE_final = Series(list_sex_cxy_image, name='CXY_IMAGE')
    CYY_IMAGE_final = Series(list_sex_cyy_image, name='CYY_IMAGE')
    ERRCXX_IMAGE_final = Series(list_sex_errcxx_image, name='ERRCXX_IMAGE')
    ERRCXY_IMAGE_final = Series(list_sex_errcxy_image, name='ERRCXY_IMAGE')
    ERRCYY_IMAGE_final = Series(list_sex_errcyy_image, name='ERRCYY_IMAGE')
    MAG_AUTO_final = Series(list_sex_mag_auto, name='MAG_AUTO')
    XWIN_IMAGE_final = Series(list_sex_xwin_image, name='XWIN_IMAGE')
    YWIN_IMAGE_final = Series(list_sex_ywin_image, name='YWIN_IMAGE')
    FLUX_AUTO_final = Series(list_sex_flux_auto, name='FLUX_AUTO')
    FLUXERR_AUTO_final = Series(list_sex_fluxerr_auto, name='FLUXERR_AUTO')
    MAGERR_AUTO_final = Series(list_sex_magerr_auto, name='MAGERR_AUTO')
    SNR_WIN_final = Series(list_sex_snr_win, name='SNR_WIN')
    X_WORLD_final = Series(list_sex_x_world, name='X_WORLD')
    Y_WORLD_final = Series(list_sex_y_world, name='Y_WORLD')
    ERRX2_WORLD_final = Series(list_sex_errx2_world, name='ERRX2_WORLD')
    ERRY2_WORLD_final = Series(list_sex_erry2_world, name='ERRY2_WORLD')
    ERRXY_WORLD_final = Series(list_sex_erryxy_world, name='ERRXY_WORLD')
    AWIN_IMAGE_final = Series(list_sex_awin_image, name='AWIN_IMAGE')
    BWIN_IMAGE_final = Series(list_sex_bwin_image, name='BWIN_IMAGE')
    THETAWIN_IMAGE_final = Series(list_sex_thetawin_image,
                                  name='THETAWIN_IMAGE')
    ERRAWIN_IMAGE_final = Series(list_sex_errawin_image, name='ERRAWIN_IMAGE')
    ERRBWIN_IMAGE_final = Series(list_sex_errbwin_image, name='ERRBWIN_IMAGE')
    ERRTHETAWIN_IMAGE_final = Series(list_sex_errthetawin_image,
                                     name='ERRTHETAWIN_IMAGE')
    FLAGS_final = Series(list_sex_flags, name='FLAGS')
    FWHM_WORLD_final = Series(list_sex_fwhm_world, name='FWHM_WORLD')
    ERRA_WORLD_final = Series(list_sex_erra_world, name='ERRA_WORLD')
    ERRB_WORLD_final = Series(list_sex_errb_world, name='ERRB_WORLD')

    final_catalogue = concat([ALPHA_J2000_final,
                              DELTA_J2000_final, NUMBER_final,
                              X_IMAGE_final, Y_IMAGE_final, X2_IMAGE_final,
                              Y2_IMAGE_final, XY_IMAGE_final,
                              ISOAREA_IMAGE_final, BACKGROUND_final,
                              THRESHOLD_final, FLUX_MAX_final, A_IMAGE_final,
                              B_IMAGE_final, THETA_IMAGE_final,
                              ERRA_IMAGE_final, ERRB_IMAGE_final,
                              FLUX_ISO_final, FLUXERR_ISO_final,
                              MAG_ISO_final, MAGERR_ISO_final, FLUX_APER_final,
                              FLUXERR_APER_final, MAG_APER_final,
                              MAGERR_APER_final, ALPHA_SKY_final,
                              DELTA_SKY_final, ERRTHETA_IMAGE_final,
                              MU_MAX_final, FWHM_IMAGE_final,
                              FLUX_RADIUS_final, ELONGATION_final,
                              ELLIPTICITY_final, CXX_IMAGE_final,
                              CXY_IMAGE_final, CYY_IMAGE_final,
                              ERRCXX_IMAGE_final, ERRCXY_IMAGE_final,
                              ERRCYY_IMAGE_final, MAG_AUTO_final,
                              XWIN_IMAGE_final, YWIN_IMAGE_final,
                              FLUX_AUTO_final, FLUXERR_AUTO_final,
                              MAGERR_AUTO_final, SNR_WIN_final, X_WORLD_final,
                              Y_WORLD_final, ERRX2_WORLD_final,
                              ERRY2_WORLD_final, ERRXY_WORLD_final,
                              AWIN_IMAGE_final, BWIN_IMAGE_final,
                              THETAWIN_IMAGE_final, ERRAWIN_IMAGE_final,
                              ERRBWIN_IMAGE_final, ERRTHETAWIN_IMAGE_final,
                              FLAGS_final, FWHM_WORLD_final, ERRA_WORLD_final,
                              ERRB_WORLD_final], axis=1)

    return final_catalogue


def get_max_mag(logger, prfs_d, conf_num, speeds):
    """

    @param logger: a logger object
    @param prfs_d:
    @param conf_num:
    @param speeds:

    @return max_values:
    """
    mags = []
    spds = []

    for speed in range(len(speeds)):
        cat_file = read_csv('{}pipeline/images_{}.csv'.format(prfs_d['fed_home'],
                                                              conf_num))
        spd = speeds[speed]
        cat_file = cat_file[cat_file['original_pm'].isin([float(spd)])]
        # cat_file.sort_values('original_mag')
        magnitudes = cat_file['original_mag'].tolist()
        detections = cat_file['detections'].tolist()

        for i in range(int(max(magnitudes)), int(min(magnitudes)), -1):
            if i != 21:  # workaround gap, mag 21 is not present at dither 1
                idx = magnitudes.index(i)
                if detections[idx] != 0:
                    mags.append(i)
                    spds.append(spd)
                    break

    max_values = list(zip(mags, spds))

    return max_values


def get_std(logger, prfs_d, conf_num):
    """ gets standard deviation of detections number

        to-do catalogue location is "too hardcoded"

    @param logger: a logger object
    @param prfs_d: a dictionary which contains all valuable data
    @param conf_num: number of catalogue to be analysed

    @return std: standard deviation represented as float variable
    """
    cat_file = read_csv('{}pipeline/images_{}.csv'.format(prfs_d['fed_home'],
                                                          conf_num))
    detections = cat_file['detections'].tolist()
    detections = asarray(detections)

    # Standard deviation
    std = sqrt(mean((detections - 1)**2))

    return std


def mag_converter(logger, input_mag, input_pm, x, unique_speeds):
    """

    @param logger: a logger object
    @param input_mag:
    @param input_pm:
    @param x:
    @param unique_speeds:

    @return mag, pm:
    """

    if int(input_pm) != 0:
        mag = input_mag
        pm = input_pm - 1000
    else:
        for idx, val in enumerate(unique_speeds):
            if val[1] == float(x):
                pm = val[0]
        mag = input_mag - 2.5

    return mag, pm


def remove_uncertainties(itered_list, unique_list):
    """

    @param itered_list:
    @param unique_list:

    @return unique_list:
    """
    for i in range(len(itered_list)):
        if itered_list[i][0] - itered_list[i][1] == 1:
            if itered_list[i][0] in unique_list:
                unique_list.remove(itered_list[i][0])
                unique_list.append(mean([itered_list[i][0],
                                         itered_list[i][1]]))
            if itered_list[i][1] in unique_list:
                unique_list.remove(itered_list[i][1])
                unique_list.append(mean([itered_list[i][0],
                                         itered_list[i][1]]))

    unique_list = list(set(unique_list))

    return unique_list


def get_input_cat(log_check, prfs_d, cat_number,
                  cat_file, complete_list):
    """ gets input values from luca's catalogue

    @param logger:
    @param prfs_d:
    @param cat_number:
    @param cat_file:
    @param complete_list: a boolean parameter True returns stars values
                          False returns SSOs

    @return sources_df: a Pandas dataframe with all the valuable data
    """
    # version = prfs_d['version'][-1:]
    version = '10'
    log_check.debug('actual version is {}'.format(version))

    if version is '3':
        catalogue_file = '{}/Cat_d{}.dat'.format(prfs_d['input_cats'],
                                                 cat_number)

        catalogue = genfromtxt(catalogue_file)
        list_x = catalogue[:, 0]
        list_y = catalogue[:, 1]
        list_mag = catalogue[:, 2]
        list_pm = catalogue[:, 3]

        x_values = []
        y_values = []
        pm_values = []
        mag_values = []

        for i in range(len(list_pm)):
            if complete_list:
                if list_pm[i] == 0.:
                    x_values.append(list_x[i])
                    y_values.append(list_y[i])
                s1 = Series(x_values, name='x_values')
                s2 = Series(y_values, name='y_values')
                sources_df = concat([s1, s2], axis=1)
            else:
                if list_pm[i] != 0.:
                    x_values.append(list_x[i])
                    y_values.append(list_y[i])
                    pm_values.append((list_pm[i] - 1000) / 10)
                    mag_values.append(list_mag[i])

                s1 = Series(x_values, name='x_values')
                s2 = Series(y_values, name='y_values')
                s3 = Series(pm_values, name='pm_values')
                s4 = Series(mag_values, name='mag_values')

                sources_df = concat([s1, s2, s3, s4], axis=1)

        return sources_df

    elif version is '5':
        cat_location = '{}/Cat_d{}.dat'.format(prfs_d['input_cats'],
                                               cat_number)

        cat_df = read_csv(cat_location, names=['x_image', 'y_image', 'mag',
                                               'pm', 'pa'], delimiter=" ")

        catalogue_file = genfromtxt(cat_location)
        list_x = catalogue_file[:, 0]
        list_y = catalogue_file[:, 1]

        unique_x = list(set(list_x))
        unique_y = list(set(list_y))

        for x_index in range(len(unique_x)):
            unique_x[x_index] = int(unique_x[x_index])

        for y_index in range(len(unique_y)):
            unique_y[y_index] = int(unique_y[y_index])

        iter_x = list(itertools.permutations(unique_x, 2))
        iter_y = list(itertools.permutations(unique_y, 2))

        unique_x = remove_uncertainties(iter_x, unique_x)
        unique_y = remove_uncertainties(iter_y, unique_y)

        unique_sources = list(itertools.product(unique_x, unique_y))

        x_values = cat_df['x_image'].tolist()
        list_x_int = ['{}'.format(int(x_value)) for x_value in x_values]
        y_values = cat_df['y_image'].tolist()
        list_y_int = ['{}'.format(int(y_value)) for y_value in y_values]
        mag_values = cat_df['mag'].tolist()
        pm_values = cat_df['pm'].tolist()
        # pa_values = cat_df['pa'].tolist()

        v1 = Series(list_x_int, name='x_values')
        v2 = Series(list_y_int, name='y_values')
        v3 = Series(mag_values, name='mag_values')
        v4 = Series(pm_values, name='pm_values')

        round_df = concat([v1, v2, v3, v4], axis=1)

        # gets a speeds list from settings file
        speeds_list = prfs_d['speeds'].split(',')
        speeds_list = [float(i) for i in speeds_list]
        speeds_list = sorted(speeds_list)

        # gets a sorted list of x positions
        unique_x_sorted = [int(j) for j in unique_x]
        unique_x_sorted = sorted(unique_x_sorted)

        # gets a list with all pairs x position and speed
        unique_speeds = []
        for pair in range(0, 18, 1):
            unique_speeds.append([speeds_list[pair],
                                  unique_x_sorted[pair]])

        sizes = []
        list_x = []
        list_y = []
        list_mag = []
        list_pm = []
        for source in range(len(unique_sources)):
            x = str(int(unique_sources[source][0]))
            y = str(int(unique_sources[source][1]))

            df_x = round_df[round_df['x_values'].isin([x])]
            df_y = df_x[df_x['y_values'].isin([y])]
            sizes.append(df_y.size / 4)
            if (df_y.size / 4) != 1 and (df_y.size / 4) != 10:
                list_x.append(x)
                list_y.append(y)
                mag, pm = mag_converter(log_check,
                                        df_y['mag_values'].iloc[0],
                                        df_y['pm_values'].iloc[0],
                                        x, unique_speeds)
                list_mag.append(mag)
                list_pm.append(pm)
            else:
                list_x.append(x)
                list_y.append(y)
                mag, pm = mag_converter(log_check,
                                        df_y['mag_values'].iloc[0],
                                        df_y['pm_values'].iloc[0],
                                        x, unique_speeds)
                list_mag.append(mag)
                list_pm.append(pm)

        s1 = Series(list_x, name='x_values', dtype=float64)
        s2 = Series(list_y, name='y_values', dtype=float64)
        s3 = Series(list_mag, name='mag_values')
        s4 = Series(list_pm, name='pm_values')

        sources_df = concat([s1, s2, s3, s4], axis=1)

        return sources_df

    elif version is '7':
        cats_list = sorted(listdir(prfs_d['input_cats']))

        catalogue_file = prfs_d['input_cats'] + '/' + cats_list[cat_number]
        catalogue = genfromtxt(catalogue_file)

        list_x = catalogue[:, 0]
        list_y = catalogue[:, 1]
        list_mag = catalogue[:, 2]
        list_pm = catalogue[:, 3]

        x_values = []
        y_values = []
        pm_values = []
        mag_values = []

        # generate indexes
        speed_0_001 = range(4, 2304, 23)
        speed_1 = range(14, 2314, 23)
        speed_3 = range(20, 2320, 23)
        speed_5 = range(21, 2321, 23)
        speed_10 = range(22, 2322, 23)

        for index in speed_0_001:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.001
        for index in speed_1:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 1
        for index in speed_3:
            list_pm[index] = list_pm[index] - 1000
        for index in speed_5:
            list_pm[index] = list_pm[index] - 1000
        for index in speed_10:
            list_pm[index] = list_pm[index] - 1000

        indexes = speed_0_001 + speed_1 + speed_3 + speed_5 + speed_10
        indexes = sorted(indexes)

        s1 = Series(list_x, name='x_values', dtype=float64)
        s2 = Series(list_y, name='y_values', dtype=float64)
        s3 = Series(list_mag, name='mag_values', dtype=float64)
        s4 = Series(list_pm, name='pm_values', dtype=float64)

        sources_df = concat([s1, s2, s3, s4], axis=1)
        sources_df = sources_df.iloc[indexes, :]

        return sources_df

    elif version is '8':
        catalogue = genfromtxt(cat_file)

        list_x = catalogue[:, 0]
        list_y = catalogue[:, 1]
        list_mag = catalogue[:, 2]
        list_pm = catalogue[:, 3]

        x_values = []
        y_values = []
        pm_values = []
        mag_values = []

        # generate indexes
        stars = range(0, 15783, 1)
        galaxies = range(15784, 170723, 1)
        speed_0_001 = range(170728, 173278, 75)
        speed_0_003 = range(170738, 173288, 75)
        speed_0_01 = range(170748, 173298, 75)
        speed_0_03 = range(170758, 173308, 75)
        speed_0_1 = range(170768, 173318, 75)
        speed_0_3 = range(170778, 173328, 75)
        speed_1 = range(170788, 173338, 75)

        speed_3 = range(170794, 173339, 75)
        speed_10 = range(170795, 173340, 75)
        speed_30 = range(170796, 173341, 75)
        speed_100 = range(170797, 173342, 75)
        speed_300 = range(170798, 173343, 75)

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

        if complete_list:
            indexes = (stars + galaxies + speed_0_001 + speed_0_003 +
                       speed_0_01 + speed_0_03 + speed_0_1 + speed_0_3 +
                       speed_1 + speed_3 + speed_10 + speed_30 +
                       speed_100 + speed_300)
            indexes = sorted(indexes)

            s1 = Series(list_x, name='x_values', dtype=float64)
            s2 = Series(list_y, name='y_values', dtype=float64)
            s3 = Series(list_mag, name='mag_values', dtype=float64)
            s4 = Series(list_pm, name='pm_values', dtype=float64)

            sources_df = concat([s1, s2, s3, s4], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            return sources_df
        else:
            "no SS0s present"
            indexes = (stars)
            indexes = sorted(indexes)

            s1 = Series(list_x, name='x_values', dtype=float64)
            s2 = Series(list_y, name='y_values', dtype=float64)
            s3 = Series(list_mag, name='mag_values', dtype=float64)
            s4 = Series(list_pm, name='pm_values', dtype=float64)

            sources_df = concat([s1, s2, s3, s4], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            return sources_df

    elif version is '10' or version is '13' or version is '12':
        log_check.debug('catalogue file is {}'.format(cat_file))
        catalogue = genfromtxt(cat_file)

        list_x = catalogue[:, 0]
        list_y = catalogue[:, 1]
        list_mag = catalogue[:, 2]
        list_pm = catalogue[:, 3]

        """
        x_values = []
        y_values = []
        pm_values = []
        mag_values = []
        """

        # generate indexes
        stars = range(prfs_d['first_star'], prfs_d['first_galaxy'] - 1, 1)
        galaxies = range(prfs_d['first_galaxy'], prfs_d['first_sso'] - 1, 1)
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

        if complete_list:
            indexes = (stars + galaxies + speed_0_001 + speed_0_003 +
                       speed_0_01 + speed_0_03 + speed_0_1 + speed_0_3 +
                       speed_1 + speed_3 + speed_10 + speed_30 +
                       speed_100 + speed_300)
            indexes = sorted(indexes)

            s1 = Series(list_x, name='x_values', dtype=float64)
            s2 = Series(list_y, name='y_values', dtype=float64)
            s3 = Series(list_mag, name='mag_values', dtype=float64)
            s4 = Series(list_pm, name='pm_values', dtype=float64)

            sources_df = concat([s1, s2, s3, s4], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            return sources_df
        else:
            # "no SS0s present"
            indexes = (stars)
            indexes = sorted(indexes)

            s1 = Series(list_x, name='x_values', dtype=float64)
            s2 = Series(list_y, name='y_values', dtype=float64)
            s3 = Series(list_mag, name='mag_values', dtype=float64)
            s4 = Series(list_pm, name='pm_values', dtype=float64)

            sources_df = concat([s1, s2, s3, s4], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            return sources_df


def get_output_catalog(logger, cat_location):
    """ gets a Pandas dataframe with all data

    @param logger: a logger object
    @param cat_location: where the study catalogue is

    @return: a Pandas dataframe with all the valuable data
    """
    # Reading full catalogue
    try:
        catalogue = fits.open(cat_location)
        catalogue_data = catalogue[2].data
        catalogue_table = Table(catalogue_data).to_pandas()
    except Exception as e:
        print(e)
        # logger.error('exception {} in file {}'.format(e, cat_location))

    return catalogue_table


def rewriting_catalogue(logger, prfs_d, cat_table, mag):
    """ write a new catalogue only fulfilled with stars

    @param logger: a logger object
    @param prfs_d:
    @param catalogue_table:

    @return True: if everything goes alright
    """
    logger.debug('reading columns from {} magnitude table'.format(mag))

    c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
                     array=cat_table['NUMBER'])
    c2 = fits.Column(name='X_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=cat_table['X_IMAGE'])
    c3 = fits.Column(name='Y_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=cat_table['Y_IMAGE'])
    c4 = fits.Column(name='X2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_table['X2_IMAGE'])
    c5 = fits.Column(name='Y2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_table['Y2_IMAGE'])
    c6 = fits.Column(name='XY_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_table['XY_IMAGE'])
    c7 = fits.Column(name='ISOAREA_IMAGE', format='1J', unit='pixel**2',
                     disp='I9', array=cat_table['ISOAREA_IMAGE'])
    c8 = fits.Column(name='BACKGROUND', format='1E', unit='count',
                     disp='G12.7', array=cat_table['BACKGROUND'])
    c9 = fits.Column(name='THRESHOLD', format='1E', unit='count',
                     disp='G12.7', array=cat_table['THRESHOLD'])
    c10 = fits.Column(name='FLUX_MAX', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUX_MAX'])
    c11 = fits.Column(name='A_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_table['A_IMAGE'])
    c12 = fits.Column(name='B_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_table['B_IMAGE'])
    c13 = fits.Column(name='THETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_table['THETA_IMAGE'])
    c14 = fits.Column(name='ERRA_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_table['ERRA_IMAGE'])
    c15 = fits.Column(name='ERRB_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_table['ERRB_IMAGE'])
    c16 = fits.Column(name='FLUX_ISO', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUX_ISO'])
    c17 = fits.Column(name='FLUXERR_ISO', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUXERR_ISO'])
    c18 = fits.Column(name='MAG_ISO', format='1E', unit='mag',
                      disp='F8.4', array=cat_table['MAG_ISO'])
    c19 = fits.Column(name='MAGERR_ISO', format='1E', unit='mag',
                      disp='F8.4', array=cat_table['MAGERR_ISO'])
    c20 = fits.Column(name='FLUX_APER', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUX_APER'])
    c21 = fits.Column(name='FLUXERR_APER', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUXERR_APER'])
    c22 = fits.Column(name='MAG_APER', format='1E', unit='mag',
                      disp='F8.4', array=cat_table['MAG_APER'])
    c23 = fits.Column(name='MAGERR_APER', format='1E', unit='mag',
                      disp='F8.4', array=cat_table['MAGERR_APER'])
    c24 = fits.Column(name='ALPHA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=cat_table['ALPHA_SKY'])
    c25 = fits.Column(name='DELTA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=cat_table['DELTA_SKY'])
    c26 = fits.Column(name='ERRTHETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_table['ERRTHETA_IMAGE'])
    c27 = fits.Column(name='MU_MAX', format='1E',
                      unit='mag * arcsec**(-2)', disp='F8.4',
                      array=cat_table['MU_MAX'])
    c28 = fits.Column(name='FWHM_IMAGE', format='1E', unit='pixel',
                      disp='F8.2', array=cat_table['FWHM_IMAGE'])
    c29 = fits.Column(name='FLUX_RADIUS', format='1E', unit='pixel',
                      disp='F10.3', array=cat_table['FLUX_RADIUS'])
    c30 = fits.Column(name='ELONGATION', format='1E', disp='F8.3',
                      array=cat_table['ELONGATION'])
    c31 = fits.Column(name='ELLIPTICITY', format='1E', disp='F8.3',
                      array=cat_table['ELLIPTICITY'])
    c32 = fits.Column(name='CXX_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_table['CXX_IMAGE'])
    c33 = fits.Column(name='CXY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_table['CXY_IMAGE'])
    c34 = fits.Column(name='CYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_table['CYY_IMAGE'])
    c35 = fits.Column(name='ERRCXX_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=cat_table['ERRCXX_IMAGE'])
    c36 = fits.Column(name='ERRCXY_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=cat_table['ERRCXY_IMAGE'])
    c37 = fits.Column(name='ERRCYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='G12.7', array=cat_table['ERRCYY_IMAGE'])
    c38 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=cat_table['MAG_AUTO'])
    c39 = fits.Column(name='XWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=cat_table['XWIN_IMAGE'])
    c40 = fits.Column(name='YWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=cat_table['YWIN_IMAGE'])
    c41 = fits.Column(name='FLUX_AUTO', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUX_AUTO'])
    c42 = fits.Column(name='FLUXERR_AUTO', format='1E', unit='count',
                      disp='G12.7', array=cat_table['FLUXERR_AUTO'])
    c43 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=cat_table['MAGERR_AUTO'])
    c44 = fits.Column(name='SNR_WIN', format='1E', disp='G10.4',
                      array=cat_table['SNR_WIN'])
    c45 = fits.Column(name='ALPHA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=cat_table['ALPHA_J2000'])
    c46 = fits.Column(name='DELTA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=cat_table['DELTA_J2000'])
    c47 = fits.Column(name='X_WORLD', format='1D', unit='deg', disp='E18.10',
                      array=cat_table['X_WORLD'])
    c48 = fits.Column(name='Y_WORLD', format='1D', unit='deg', disp='E18.10',
                      array=cat_table['Y_WORLD'])
    c49 = fits.Column(name='ERRX2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_table['ERRX2_WORLD'])
    c50 = fits.Column(name='ERRY2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_table['ERRY2_WORLD'])
    c51 = fits.Column(name='ERRXY_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_table['ERRXY_WORLD'])
    c52 = fits.Column(name='AWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_table['AWIN_IMAGE'])
    c53 = fits.Column(name='BWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_table['BWIN_IMAGE'])
    c54 = fits.Column(name='THETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_table['THETAWIN_IMAGE'])
    c55 = fits.Column(name='ERRAWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_table['ERRAWIN_IMAGE'])
    c56 = fits.Column(name='ERRBWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_table['ERRBWIN_IMAGE'])
    c57 = fits.Column(name='ERRTHETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_table['ERRTHETAWIN_IMAGE'])
    c58 = fits.Column(name='FLAGS', format='1I', disp='I3',
                      array=cat_table['FLAGS'])
    c59 = fits.Column(name='FWHM_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_table['FWHM_WORLD'])
    c60 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_table['ERRA_WORLD'])
    c61 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_table['ERRB_WORLD'])

    coldefs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                            c13, c14, c15, c16, c17, c18, c19, c20, c21, c22,
                            c23, c24, c25, c26, c27, c28, c29, c30, c31, c32,
                            c33, c34, c35, c36, c37, c38, c39, c40, c41, c42,
                            c43, c44, c45, c46, c47, c48, c49, c50, c51, c52,
                            c53, c54, c55, c56, c57, c58, c59, c60, c61])
    tbhdu = fits.BinTableHDU.from_columns(coldefs)

    logger.debug('reading columns from {} magnitude table'.format(mag))
    cat_location = prfs_d['output_cats'] + '/catalogue_{}.cat'.format(mag)

    catalogue = fits.open(cat_location)
    logger.debug('saving new table from {} magnitude'.format(mag))
    catalogue[2] = tbhdu
    catalogue[2].header['EXTNAME'] = 'LDAC_OBJECTS'

    new_cat = prfs_d['output_cats'] + '/catalogue_{}.cat'.format(mag)
    logger.debug('writing cat {}/catalogue_{}.cat'.format(prfs_d['output_cats'],
                                                          mag))
    catalogue.writeto(new_cat, overwrite=True)

    return True


def cut_catalog(catalog, keys, margin):
    """ returns a resized version of input catalog

    @param catalog:
    @param keys:
    @param margin:

    @return cutted_catalog
    """
    for key_ in keys:
        catalog = catalog[(key_[0] + margin) > catalog[key_[1]]]
        catalog = catalog[(key_[0] - margin) < catalog[key_[1]]]

    return catalog


def shrink_catalog(catalog, margins_d):
    """

    @param catalog:
    @param margins_d:

    @return catalog:
    """

    catalog = catalog[catalog['ALPHA_J2000'] < margins_d['above_ra']]
    catalog = catalog[catalog['ALPHA_J2000'] > margins_d['below_ra']]
    catalog = catalog[catalog['DELTA_J2000'] < margins_d['above_dec']]
    catalog = catalog[catalog['DELTA_J2000'] > margins_d['below_dec']]

    return catalog


def look_for_ssos(logger, prfs_d, mag, scmp_d, f_conf):
    """

    tolerance has been adjusted to 1'5 arcseconds -> 0.000416666
    0'75 -> 0.000208333

    @param logger:
    @param prfs_d:
    @param mag:
    @param scmp_d
    """
    sex_keys = [['sex_number', 'NUMBER'], ['sex_x_image', 'X_IMAGE'],
                ['sex_y_image', 'Y_IMAGE'],
                ['sex_alpha_sky', 'ALPHA_SKY'], ['sex_delta_sky', 'DELTA_SKY'],
                ['sex_snr_win', 'SNR_WIN'],
                ['sex_alpha_j2000', 'ALPHA_J2000'],
                ['sex_delta_j2000', 'DELTA_J2000'], ['sex_x_world', 'X_WORLD'],
                ['sex_y_world', 'Y_WORLD']]

    scmp_keys = [['scamp_source_number', 'SOURCE_NUMBER'],
                 ['scmp_x_image', 'X_IMAGE'], ['scmp_y_image', 'Y_IMAGE'],
                 ['scmp_alpha_j2000', 'ALPHA_J2000'],
                 ['scmp_delta_j2000', 'DELTA_J2000'],
                 ['scmp_epoch', 'EPOCH'], ['scmp_mag', 'MAG']]

    out_dict = {}
    out_dict['sources'] = []
    out_dict['dithers'] = []
    out_dict['sextractor'] = []
    out_dict['sextractor_flag'] = []
    out_dict['scamp'] = []
    out_dict['scamp_flag'] = []
    out_dict['scamp_distance'] = []
    out_dict['CCD'] = []

    for key_ in sex_keys:
        out_dict[key_[0]] = []

    for key_ in scmp_keys:
        out_dict[key_[0]] = []

    unique_fits = get_fits(True)  # Gets fits from fits folder

    for fits_idx in range(0, len(unique_fits), 8):
        try:
            ssos_j = []
            for proc in range(0, 9, 1):
                idx = fits_idx + proc  # index

                image = unique_fits[idx]
                ssos_p = Process(target=look_for_ssos_thread,
                                 args=(logger, prfs_d, image, out_dict,
                                       sex_keys, scmp_keys, mag, scmp_d,
                                       f_conf,))
                ssos_j.append(ssos_p)
                ssos_p.start()

            active_ssos = list([job.is_alive() for job in ssos_j])
            while True in active_ssos:
                active_ssos = list([job.is_alive() for job in ssos_j])
                pass
        except IndexError:
            logger.debug('extraction finished')

    return True


def look_for_ssos_thread(logger, prfs_d, image, out_dict, sex_keys,
                         scamp_keys, mag, scmp_d, f_conf):
    """

    @param logger:
    @param prfs_d:
    @param image:
    @param out_dict:
    @param sex_keys:
    @param scamp_keys:
    @param mag:
    @param scmp_d:
    """
    # Booleans variables
    logger.debug('loading image {}'.format(image))

    # Gets sky limits of each image
    d_limits = {}
    for d in range(1, 5, 1):
        i_image = '{}/{}{}{}'.format(prfs_d['fits_dir'], image[:15],
                                     d, image[-5:])
        d_limits['{:d}'.format(d)] = get_fits_limits(i_image)

    save = True  # Not file needed
    complete = True
    # Gets ra/dec coordinates for each dither
    d_ssos = {}
    for d in range(1, 5, 1):
        d_cat = prfs_d['input_cats'] + '/Cat_20-21_d{}.dat'.format(d)
        d_ssos['{}'.format(d)] = Create_regions(d_cat, prfs_d).luca(save,
                                                                    complete)

    # Substract not present SSOs in image selected
    cat_list = []
    for d in range(1, 5, 1):
        idx = str(d)
        limits = d_limits[idx]  # get fits limits
        out_cat = d_ssos['{}'.format(d)]
        out_cat = out_cat[limits['above_ra'] > out_cat['alpha_j2000']]
        out_cat = out_cat[limits['below_ra'] < out_cat['alpha_j2000']]
        out_cat = out_cat[limits['above_dec'] > out_cat['delta_j2000']]
        out_cat = out_cat[limits['below_dec'] < out_cat['delta_j2000']]
        cat_list.append(out_cat)

    # Merge all catalog dither into a single one
    ssos = concat(cat_list, ignore_index=True)
    ssos_n = '{}_{}_{}'.format(image[:13], f_conf, mag)
    ssos.to_csv('{}/{}_ssos.csv'.format(prfs_d['catalogs_dir'], ssos_n))

    # Creates a list populated with all single sources
    unique_sources = list(set(ssos['source'].tolist()))

    # Creates lists for final catalogue in out_dict
    # TODO out_dict should be a shared dict

    counter_source = 0
    # Loops over CCD sources - Custom catalog
    for source_number in unique_sources:

        logger.debug('look for source {}'.format(source_number))
        counter_source = counter_source + 1
        tmp_cat = ssos.loc[ssos['source'] == source_number]
        # Looks for duplicate sources in same dither
        if tmp_cat['source'].size is not tmp_cat['dither_values'].size:
            raise Exception

        logger.debug('checking {} of {}/{}'.format(source_number,
                                                   counter_source,
                                                   len(unique_sources)))
        for d in tmp_cat['dither_values'].tolist():
            idx = str(int(d))
            logger.debug('checking dither {} of {}'.format(d, source_number))
            sex_flag = False
            sex_counter = 0
            scmp_flag = False
            scmp_counter = 0

            # Margin represents how much catalog will be cutted
            margin = 4 * prfs_d['tolerance']

            # Appending values to lists
            out_dict['sources'].append(source_number)
            out_dict['dithers'].append(d)

            # Looks for a particular dither
            tmp_out = tmp_cat[tmp_cat['dither_values'].isin([d])]
            # Looks for correspondence in sextractor catalogue
            # Coordinates to be founded
            if tmp_out['alpha_j2000'].size > 1:
                raise Exception
            if tmp_out['delta_j2000'].size > 1:
                raise Exception

            (alpha, delta) = (tmp_out['alpha_j2000'].iloc[0],
                              tmp_out['delta_j2000'].iloc[0])

            # Gets limits of desired image
            limits = d_limits[idx]

            # Open sextractor catalog
            sex_cat = '{}/{}{}.cat'.format(prfs_d['output_cats'], image[:15],
                                           int(d))
            hdu_list = fits.open(sex_cat)
            sex_table = Table(hdu_list[2].data).to_pandas()

            margins = [[alpha, 'X_WORLD'], [delta, 'Y_WORLD']]
            sex_table = cut_catalog(sex_table, margins, margin)

            # loop over sextractor output
            for sex_source in sex_table['NUMBER'].tolist():
                t = sex_table[sex_table['NUMBER'].isin([sex_source])].iloc[0]
                ra = t['X_WORLD']
                dec = t['Y_WORLD']

                close, distance = check_distance(alpha, ra, delta, dec,
                                                 prfs_d['tolerance'])
                if close:
                    df = sex_table[sex_table['NUMBER'].isin([sex_source])]
                    sex_flag = True
                    sex_counter = sex_counter + 1

            if sex_counter > 1:
                out_dict['sextractor_flag'].append(1)
            else:
                out_dict['sextractor_flag'].append(0)

            if sex_flag:
                out_dict['sextractor'].append(True)
                out_dict['CCD'].append(image[8:13])
                for key_ in sex_keys:
                    out_dict[key_[0]].append(df[key_[1]].iloc[0])
            else:
                out_dict['sextractor'].append(False)
                out_dict['CCD'].append(image[8:13])
                for key_ in sex_keys:
                    out_dict[key_[0]].append(False)

            # open scamp catalog
            full_cat_n = '{}_{}'.format(f_conf, mag)
            scamp_cat = '{}/full_{}_1.cat'.format(prfs_d['results_dir'],
                                                  full_cat_n)
            hdu_list = fits.open(scamp_cat)
            scmp_table = Table(hdu_list[2].data).to_pandas()
            dithers = [range(1, 10, 1), range(10, 19, 1),
                       range(19, 28, 1), range(28, 37, 1)]
            scmp_table = scmp_table[scmp_table['CATALOG_NUMBER'].isin(dithers[int(d)-1])]

            # scamp's data is for the entire fpa, we should reduce it
            margins = [[alpha, 'ALPHA_J2000'], [delta, 'DELTA_J2000']]
            scmp_table = cut_catalog(scmp_table, margins, margin)

            # loop over scamp output
            scmp_dataframes = {}    # Sometimes we get multiple results for the
                                    # same coordinates so we will store them
                                    # and get the closer one
            scmp_distances = {}  # Need also a place for store distances
            for scmp_source in scmp_table['SOURCE_NUMBER'].tolist():
                t = scmp_table[scmp_table['SOURCE_NUMBER'].isin([scmp_source])].iloc[0]
                ra = t['ALPHA_J2000']
                dec = t['DELTA_J2000']

                close, distance = check_distance(alpha, ra, delta, dec,
                                                 prfs_d['tolerance'] * 2)

                if close:
                    df_2 = scmp_table[scmp_table['SOURCE_NUMBER'].isin([scmp_source])]
                    scmp_dataframes[scmp_counter] = df_2
                    scmp_distances[scmp_counter] = distance
                    scmp_flag = True
                    scmp_counter = scmp_counter + 1

            if scmp_counter > 1:
                out_dict['scamp_flag'].append(1)
            else:
                out_dict['scamp_flag'].append(0)

            if scmp_flag:
                i = scmp_distances.values().index(min(scmp_distances.values()))
                df_3 = scmp_dataframes[scmp_dataframes.keys()[i]]
                # Append distance between input catalog and scamp's catalog
                distance_seconds = min(scmp_distances.values()) / 0.00027777778
                out_dict['scamp_distance'].append(distance_seconds)

                # TODO Get scamp proper motion value

                # out_dict['scamp_pm'].append('pm')
                out_dict['scamp'].append(True)
                for key_ in scamp_keys:
                    out_dict[key_[0]].append(df_3[key_[1]].iloc[0])
            else:
                out_dict['scamp_distance'].append(False)
                out_dict['scamp'].append(False)
                for key_ in scamp_keys:
                    out_dict[key_[0]].append(False)

    ssos_Cat = DataFrame(out_dict)
    sso_cat_n = '{}_{}_{}'.format(image[:13], f_conf, mag)
    ssos_Cat.to_csv('{}/{}_sso_cat.csv'.format(prfs_d['catalogs_dir'],
                                               sso_cat_n))


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


def merge_sso_cat(logger, prfs_d, conf, mag):
    """ creates a single sso_cat.csv file for each configuration

    @param logger:
    @param prfs_d:
    @param conf:
    @param mag:

    @return dataframe: a dataframe with all the merged data
    """
    file_name = '_{}_{}_{}_{}_{}_'.format(conf[0], conf[1], conf[2],
                                          conf[3], mag)
    idx = 0
    csv_files = []
    dataframes_dict = {}
    dataframes_list = []

    files = listdir(prfs_d['home'] + '/pipeline')
    for csv_ in files:
        if csv_[-12:] == '_sso_cat.csv' and file_name in csv_:
            idx += 1
            csv_files.append(csv_)
            dataframes_dict[idx] = read_csv(csv_, index_col=0)

    # Remove old files
    for file_ in csv_files:
        remove(file_)

    for key_ in dataframes_dict.keys():
        dataframes_list.append(dataframes_dict[key_])

    output_file_name = file_name[1:] + 'sso_cat.csv'

    dataframe = concat(dataframes_list)
    dataframe.to_csv(output_file_name)

    logger.debug('output catalog {} created'.format(output_file_name))

    return dataframe


def merge_ssos(logger, prfs_d, conf, mag):
    """ creates a single ssos.csv file for each configuration

    """
    file_name = '_{}_{}_{}_{}_{}_'.format(conf[0], conf[1], conf[2],
                                          conf[3], mag)
    print(file_name)
    idx = 0
    csv_files = []
    dataframes_dict = {}
    dataframes_list = []

    files = listdir(prfs_d['home'] + '/pipeline')
    for csv_ in files:
        if csv_[-9:] == '_ssos.csv' and file_name in csv_:
            idx += 1
            csv_files.append(csv_)
            dataframes_dict[idx] = read_csv(csv_, index_col=0)

    # Remove old files
    for file_ in csv_files:
        remove(file_)

    for key_ in dataframes_dict.keys():
        dataframes_list.append(dataframes_dict[key_])

    output_file_name = file_name[1:] + 'ssos.csv'

    dataframe = concat(dataframes_list)
    dataframe.to_csv(output_file_name)

    return dataframe

    logger.debug('output catalog {} created'.format(output_file_name))


def merge_catalog(prfs_d, folder_sex):
    """

    @param prfs_d:
    @param folder_sex:
    @return cat_n: a reference in string format.
    """
    """
    cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                 folder_scmp)
    """
    fits_files = get_fits(unique=True)

    # Using a mutable dictionary instead a fixed lists will allow us
    # to change catalog columns without a problem
    # Needs
    cat_dict = {}
    for idx in range(0, len(fits_files), 1):
        cat_loc = '{}/{}/{}.cat'.format(prfs_d['fits_dir'], folder_sex,
                                        fits_files[idx][:-5])
        # Gets data from CCD catalog
        cat_data = get_output_catalog(prfs_d, cat_loc)
        cat_n = fits_files[idx][:-5]
        cat_dict[cat_n] = {}  # Creates a dict for selected CCD
        for column_ in cat_data.columns:
            cat_dict[cat_n][column_] = cat_data[column_].tolist()

    cat_final = {}
    for ccd_key in cat_dict:
        for param_key in cat_dict[ccd_key].keys():
            for value_ in cat_dict[ccd_key][param_key]:
                try:
                    cat_final[param_key].append(value_)
                except KeyError:
                    cat_final[param_key] = []
                    cat_final[param_key].append(value_)

    cat_df = DataFrame(cat_final)
    rewrite_catalog(cat_df, prfs_d, folder_sex)
        # cat_df.to_csv('test.csv')


def rewrite_catalog(cat_df, prfs_d, folder_sex):
    """
    writes a new catalogue only fulfilled with stars
    @param logger: a logger object
    @param prfs_d:
    @param catalogue_table:
    @return True: if everything goes alright
    """
    c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
                     array=cat_df['NUMBER'])
    c2 = fits.Column(name='X_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=cat_df['X_IMAGE'])
    c3 = fits.Column(name='Y_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=cat_df['Y_IMAGE'])
    c4 = fits.Column(name='X2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_df['X2_IMAGE'])
    c5 = fits.Column(name='Y2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_df['Y2_IMAGE'])
    c6 = fits.Column(name='XY_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_df['XY_IMAGE'])
    c7 = fits.Column(name='ISOAREA_IMAGE', format='1J', unit='pixel**2',
                     disp='I9', array=cat_df['ISOAREA_IMAGE'])
    c8 = fits.Column(name='BACKGROUND', format='1E', unit='count',
                     disp='G12.7', array=cat_df['BACKGROUND'])
    c9 = fits.Column(name='THRESHOLD', format='1E', unit='count',
                     disp='G12.7', array=cat_df['THRESHOLD'])
    c10 = fits.Column(name='FLUX_MAX', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUX_MAX'])
    c11 = fits.Column(name='A_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_df['A_IMAGE'])
    c12 = fits.Column(name='B_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_df['B_IMAGE'])
    c13 = fits.Column(name='THETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_df['THETA_IMAGE'])
    c14 = fits.Column(name='ERRA_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_df['ERRA_IMAGE'])
    c15 = fits.Column(name='ERRB_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_df['ERRB_IMAGE'])
    c16 = fits.Column(name='FLUX_ISO', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUX_ISO'])
    c17 = fits.Column(name='FLUXERR_ISO', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUXERR_ISO'])
    c18 = fits.Column(name='MAG_ISO', format='1E', unit='mag',
                      disp='F8.4', array=cat_df['MAG_ISO'])
    c19 = fits.Column(name='MAGERR_ISO', format='1E', unit='mag',
                      disp='F8.4', array=cat_df['MAGERR_ISO'])
    c20 = fits.Column(name='FLUX_APER', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUX_APER'])
    c21 = fits.Column(name='FLUXERR_APER', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUXERR_APER'])
    c22 = fits.Column(name='MAG_APER', format='1E', unit='mag',
                      disp='F8.4', array=cat_df['MAG_APER'])
    c23 = fits.Column(name='MAGERR_APER', format='1E', unit='mag',
                      disp='F8.4', array=cat_df['MAGERR_APER'])
    c24 = fits.Column(name='ALPHA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=cat_df['ALPHA_SKY'])
    c25 = fits.Column(name='DELTA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=cat_df['DELTA_SKY'])
    c26 = fits.Column(name='ERRTHETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_df['ERRTHETA_IMAGE'])
    c27 = fits.Column(name='MU_MAX', format='1E',
                      unit='mag * arcsec**(-2)', disp='F8.4',
                      array=cat_df['MU_MAX'])
    c28 = fits.Column(name='FWHM_IMAGE', format='1E', unit='pixel',
                      disp='F8.2', array=cat_df['FWHM_IMAGE'])
    c29 = fits.Column(name='FLUX_RADIUS', format='1E', unit='pixel',
                      disp='F10.3', array=cat_df['FLUX_RADIUS'])
    c30 = fits.Column(name='ELONGATION', format='1E', disp='F8.3',
                      array=cat_df['ELONGATION'])
    c31 = fits.Column(name='ELLIPTICITY', format='1E', disp='F8.3',
                      array=cat_df['ELLIPTICITY'])
    c32 = fits.Column(name='CXX_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_df['CXX_IMAGE'])
    c33 = fits.Column(name='CXY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_df['CXY_IMAGE'])
    c34 = fits.Column(name='CYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_df['CYY_IMAGE'])
    c35 = fits.Column(name='ERRCXX_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=cat_df['ERRCXX_IMAGE'])
    c36 = fits.Column(name='ERRCXY_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=cat_df['ERRCXY_IMAGE'])
    c37 = fits.Column(name='ERRCYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='G12.7', array=cat_df['ERRCYY_IMAGE'])
    c38 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=cat_df['MAG_AUTO'])
    c39 = fits.Column(name='XWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=cat_df['XWIN_IMAGE'])
    c40 = fits.Column(name='YWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=cat_df['YWIN_IMAGE'])
    c41 = fits.Column(name='FLUX_AUTO', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUX_AUTO'])
    c42 = fits.Column(name='FLUXERR_AUTO', format='1E', unit='count',
                      disp='G12.7', array=cat_df['FLUXERR_AUTO'])
    c43 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=cat_df['MAGERR_AUTO'])
    c44 = fits.Column(name='SNR_WIN', format='1E', disp='G10.4',
                      array=cat_df['SNR_WIN'])
    c45 = fits.Column(name='ALPHA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=cat_df['ALPHA_J2000'])
    c46 = fits.Column(name='DELTA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=cat_df['DELTA_J2000'])
    c47 = fits.Column(name='X_WORLD', format='1D', unit='deg',
                      disp='E18.10', array=cat_df['X_WORLD'])
    c48 = fits.Column(name='Y_WORLD', format='1D', unit='deg',
                      disp='E18.10', array=cat_df['Y_WORLD'])
    c49 = fits.Column(name='ERRX2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_df['ERRX2_WORLD'])
    c50 = fits.Column(name='ERRY2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_df['ERRY2_WORLD'])
    c51 = fits.Column(name='ERRXY_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_df['ERRXY_WORLD'])
    c52 = fits.Column(name='AWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_df['AWIN_IMAGE'])
    c53 = fits.Column(name='BWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_df['BWIN_IMAGE'])
    c54 = fits.Column(name='THETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_df['THETAWIN_IMAGE'])
    c55 = fits.Column(name='ERRAWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_df['ERRAWIN_IMAGE'])
    c56 = fits.Column(name='ERRBWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_df['ERRBWIN_IMAGE'])
    c57 = fits.Column(name='ERRTHETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_df['ERRTHETAWIN_IMAGE'])
    c58 = fits.Column(name='FLAGS', format='1I', disp='I3',
                      array=cat_df['FLAGS'])
    c59 = fits.Column(name='FWHM_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_df['FWHM_WORLD'])
    c60 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_df['ERRA_WORLD'])
    c61 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_df['ERRB_WORLD'])

    coldefs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,
                            c12, c13, c14, c15, c16, c17, c18, c19, c20,
                            c21, c22, c23, c24, c25, c26, c27, c28, c29,
                            c30, c31, c32, c33, c34, c35, c36, c37, c38,
                            c39, c40, c41, c42, c43, c44, c45, c46, c47,
                            c48, c49, c50, c51, c52, c53, c54, c55, c56,
                            c57, c58, c59, c60, c61])
    tbhdu = fits.BinTableHDU.from_columns(coldefs)

    cat_ref = '{}/{}/mag_20-21_CCD_x1_y1_d1.cat'.format(prfs_d['fits_dir'],
                                                        folder_sex)

    catalog = fits.open(cat_ref)
    catalog[2] = tbhdu
    catalog[2].header['EXTNAME'] = 'LDAC_OBJECTS'
    cat_loc = '{}/{}/catalog.cat'.format(prfs_d['fits_dir'], folder_sex)
    catalog.writeto(cat_loc, overwrite=True)


def check_source(catalog_n, i_df, o_alpha, o_delta):
    """

    :param catalog_n:
    :param i_df:0
    :param o_alpha:
    :param o_delta:
    :return:
    """
    tolerance = 0.0000833  # 0.3 arcsecond

    if catalog_n != 0:
        i_df = i_df[i_df['catalog'].isin([catalog_n])]

    i_df = i_df[i_df['alpha_j2000'] + tolerance > o_alpha]
    i_df = i_df[o_alpha > i_df['alpha_j2000'] - tolerance]
    i_df = i_df[i_df['delta_j2000'] + tolerance > o_delta]
    i_df = i_df[o_delta > i_df['delta_j2000'] - tolerance]

    return i_df


def get_input_dicts(mag):
    """ FIXME bug! why should I set each time catalog name?

    :param mag:
    :return: i_df
    """
    prfs_d = extract_settings()

    # Input sources
    input_cat_d = {}
    for d in range(1, 5, 1):
        cat_loc = '{}/{}/Catalogs'.format(prfs_d['fits_dir'], mag)
        cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
        input_cat_d[d] = '{}.dat'.format(cat_name)
    input_stars_d = Create_regions(input_cat_d).check_stars(mag, True)

    input_cat_d = {}
    for d in range(1, 5, 1):
        cat_loc = '{}/{}/Catalogs'.format(prfs_d['fits_dir'], mag)
        cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
        input_cat_d[d] = '{}.dat'.format(cat_name)
    input_galaxies_d = Create_regions(input_cat_d).check_galaxies(mag, True)

    input_cat_d = {}
    for d in range(1, 5, 1):
        cat_loc = '{}/{}/Catalogs'.format(prfs_d['fits_dir'], mag)
        cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
        input_cat_d[d] = '{}.dat'.format(cat_name)
    input_ssos_d = Create_regions(input_cat_d).check_ssos(mag, True)

    # Creates a DataFrame from an input dictionary
    input_stars_l = []
    for key_ in input_stars_d.keys():
        input_stars_l.append(input_stars_d[key_])

    i_stars_df = concat(input_stars_l, axis=0)
    i_stars_df = i_stars_df.reset_index(drop=True)

    # Creates a DataFrame from an input dictionary
    input_galaxies_l = []
    for key_ in input_galaxies_d.keys():
        input_galaxies_l.append(input_galaxies_d[key_])

    i_galaxies_df = concat(input_galaxies_l, axis=0)
    i_galaxies_df = i_galaxies_df.reset_index(drop=True)

    # Creates a DataFrame from an input dictionary
    input_ssos_l = []
    for key_ in input_ssos_d.keys():
        input_ssos_l.append(input_ssos_d[key_])

    i_ssos_df = concat(input_ssos_l, axis=0)
    i_ssos_df = i_ssos_df.reset_index(drop=True)

    i_df = {'stars': i_stars_df, 'galaxies': i_galaxies_df, 'ssos': i_ssos_df}

    return i_df
