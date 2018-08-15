#!/usr/bin/python
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from misc import setting_logger, extract_settings
from cats_management import get_input_catalogue
from subprocess import Popen
from os import remove, listdir
from pandas import concat, Series
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
# from astropy.nddata.utils import NoOverlapError


# def disable_graphics(logger):
#     cmd_output_1 = 'export DISPLAY=:1'
#     cmd_output_2 = ' Xvfb :1 -screen 0 1024x768x16 &'

#     logger.debug('disabling graphical output 1')
#     disable_process_1 = Popen(cmd_output_1, shell=True)
#     disable_process_1.wait()

#     logger.debug('disabling graphical output 2')
#     disable_process_2 = Popen(cmd_output_2, shell=True)
#     disable_process_2.wait()

#     return True

def database_extract(logger, prefs_dict, cat_number):
    """ returns a pandas dataframe fulfilled with desired catalog.

    @param logger: a logger object.
    @param prefs_dict: a preferences dictionary object.
    @param cat_number: catalog number choosen.

    @return db: dataframe filtered.
    """
    full_catalogue = fits.open(prefs_dict['results_dir'] + '/full_1.cat')
    data_table = Table(full_catalogue[2].data).to_pandas()

    db = data_table.loc[data_table['CATALOG_NUMBER'].isin([int(cat_number)])]

    return db


def clear_data(logger, prefs_dict):
    """

    @param logger:
    @param prefs_dict:
    """
    logger.info('Removing old data')

    for cat_num in range(1, 5, 1):
        files_in_dir = listdir(prefs_dict['images_out'] + '/' + str(cat_num))
        for file_index in range(len(files_in_dir)):
            file_name = '{}/{}/{}'.format(prefs_dict['images_out'], cat_num,
                                          files_in_dir[file_index])
            remove(file_name)
            logger.debug('Removing {}'.format(file_name))

        files_in_dir = listdir(prefs_dict['fits_out'] + '/' + str(cat_num))
        for file_index in range(len(files_in_dir)):
            file_name = '{}/{}/{}'.format(prefs_dict['fits_out'], cat_num,
                                          files_in_dir[file_index])
            remove(file_name)
            logger.debug('Removing {}'.format(file_name))


def images_extraction(logger, prefs_dict, cat_number):
    t1 = get_input_catalogue(prefs_dict, cat_number, False)
    ok = 0
    no = 0

    source_number = []
    has_image = []
    location_image = []
    location_cutout = []

    for source in range(int(t1['x_values'].size)):
        ccd_image = '{}/CCD_x0_y0_d{}.fits'.format(prefs_dict['fits_dir'],
                                                   cat_number)
        logger.debug('looking for image {}'.format(ccd_image))
        f = fits.open(ccd_image)
        w = WCS(f[0].header)

        # Get object's position
        x_position = t1['x_values'].loc[source]
        y_position = t1['y_values'].loc[source]
        position = (x_position, y_position)

        if 0 < x_position < 4096 and 0 < y_position < 4136:
            shape = (140, 140)
            cutout = Cutout2D(f[0].data, position, shape, wcs=w)
            hdu = fits.PrimaryHDU(data=cutout.data,
                                  header=cutout.wcs.to_header())

            output_dir = '{}/{}/'.format(prefs_dict['fits_out'], cat_number)
            output_name = 'source_{}.fits'.format(source)
            output_fits = output_dir + output_name
            logger.debug('writing image {}'.format(output_fits))
            hdu.writeto(output_fits)

            # Output data
            cutout_dir = '{}/{}/'.format(prefs_dict['images_out'], cat_number)
            cutout_name = 'source_{}.jpg'.format(source)
            output_cutout = cutout_dir + cutout_name

            source_number.append(source + 1)
            has_image.append(True)
            location_image.append(output_fits)
            location_cutout.append(output_cutout)
            ok = ok + 1
        else:
            # Output data
            source_number.append(source + 1)
            has_image.append(False)
            location_image.append(' ')
            location_cutout.append(' ')
            no = no + 1

    s1 = Series(source_number, name='source_num')
    s2 = Series(has_image, name='has_image')
    s3 = Series(location_image, name='where_is')
    s4 = Series(location_cutout, name='where_out')

    images_table = concat([s1, s2, s3, s4], axis=1)

    logger.info('{} SSOs in for {} dither'.format(ok, cat_number))
    logger.info('{} SSOs out for {} dither'.format(no, cat_number))

    return images_table


def images_conversion(logger, prefs_dict, table):
    logger.info('Extraction process started')
    # Needs to read configuration from file

    for source in range(table['source_num'].size):
        if table['has_image'].iloc[source]:
            input_image = table['where_is'].iloc[source]
            output_image = table['where_out'].iloc[source]
            cmd_11 = 'ds9 {} -zoom to fit -histequ'.format(input_image)
            cmd_12 = ' -saveimage jpeg {} -exit'.format(output_image)
            cmd_1 = cmd_11 + cmd_12

            ds9_process = Popen(cmd_1, shell=True)
            ds9_process.wait()


def pdf_creation(logger, table_dicts, prefs_dict):
    """

    """

    x_page = [0, 10, 10, 10, 10]
    y_page = [0, 625, 420, 215, 10]

    c = canvas.Canvas(prefs_dict['report_out'], pagesize=A4)

    for source in range(table_dicts[1]['source_num'].size):
        for cat_num in range(1, 5, 1):
            has_image = table_dicts[cat_num]['has_image'].iloc[source]
            database_extract(logger, prefs_dict, cat_num)
            if has_image:
                output_image = table_dicts[cat_num]['where_out'].iloc[source]
                logger.debug('placing image {}'.format(output_image))
                c.drawImage(output_image, x_page[cat_num], y_page[cat_num],
                            width=288, height=175)

                text = output_image
                x, y = x_page[cat_num] + 300, y_page[cat_num]
                logger.debug('writing legend for {}'.format(output_image))
                c.drawString(x, y, text)
                c.setFont("Helvetica", 12)
            else:
                output_image = '/pcdisk/holly/sgongora/Documents/Euclid/empty_sso.png'
                logger.debug('placing image {}'.format(output_image))
                c.drawImage(output_image, x_page[cat_num], y_page[cat_num],
                            width=288, height=175)
                text = "Empty image"
                x, y = x_page[cat_num] + 300, y_page[cat_num]
                logger.debug('writing legend for {}'.format(output_image))
                c.drawString(x, y, text)
                c.setFont("Helvetica", 12)
        c.showPage()
    c.save()


if __name__ == '__main__':
    logger = setting_logger()
    prefs_dict = extract_settings()
    table_dicts = {}

    # disable_graphics(logger)
    clear_data(logger, prefs_dict)
    for cat_number in range(1, 5, 1):
        table_dicts[cat_number] = images_extraction(logger, prefs_dict,
                                                    cat_number)

        database_extract(logger, prefs_dict, cat_number)
        images_conversion(logger, prefs_dict, table_dicts[cat_number])
    pdf_creation(logger, table_dicts, prefs_dict)
