#!/usr/bin/python
# -*- coding: utf-8 -*-

from subprocess import Popen


def ds9_visualization(logger, prfs_d, table_dicts, fits_num):
    """ create a jpeg image from a fits one

    @param logger:
    @param prfs_d:
    @param table_dicts:
    @param fits_num: catalog number

    @return where_out:
    """

    for source_fits in range(table_dicts[str(fits_num)]['where_is'].size):
        input_fits = table_dicts[str(fits_num)]['where_is'].loc[source_fits]
        output_image = table_dicts[str(fits_num)]['where_out'].loc[source_fits]

        total_number = table_dicts[str(fits_num)]['where_is'].size
        logger.debug('saving image {} {}/{}'.format(output_image, source_fits,
                                                    total_number))

        cmd_11 = 'ds9 {} -zoom to fit -histequ '.format(input_fits)
        cmd_12 = '-regions -format xy '
        # load inputs regions
        cmd_13 = '{}pipeline/output_{}.reg '.format(prfs_d['fed_home'],
                                                    fits_num)
        cmd_14 = '-saveimage jpeg {} -exit'.format(output_image)
        cmd_1 = cmd_11 + cmd_12 + cmd_13 + cmd_14

        ds9_process = Popen(cmd_1, shell=True)
        ds9_process.wait()

    return True
