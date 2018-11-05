# !/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1

Todo:
    * Improve log messages
    * Get out check_source

*GNU Terry Pratchett*
"""
from multiprocessing import Process
from os import listdir

from astropy.io import fits
from astropy.time import Time

from misc import extract_settings_elvis, setting_logger

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def change_times():
    """
    hardcoded times:
    2021-06-26T09:00:00.00000
    2021-06-26T09:16:43.00000
    2021-06-26T09:33:26.00000
    2021-06-26T09:50:09.00000

    hardcoded dir:

    """
    prfs_d = extract_settings_elvis()
    logger = setting_logger(prfs_d)

    core_number = int(prfs_d['cores_number'])
    files = listdir(prfs_d['fits_dir'])
    fits_files = []

    for file_ in files:
        if file_[-5:] == '.fits':
            fits_files.append(file_)

    for idx_fits in range(0, len(fits_files), core_number):
        try:
            time_j = []
            for proc in range(0, core_number, 1):
                idx = idx_fits + proc
                time_p = Process(target=change_times_thread,
                                 args=(prfs_d, fits_files[idx],))
                time_j.append(time_p)
                time_p.start()

            active_time = list([job.is_alive() for job in time_j])
            while True in active_time:
                active_time = list([job.is_alive() for job in time_j])
                pass
        except IndexError:
            logger.debug('extraction finished')

    return True


def change_times_thread(prfs_d, fits_image):
    """

    @param prfs_d:
    @param fits_image:

    @return True:
    """
    data, header = fits.getdata('{}/{}'.format(prfs_d['fits_dir'], fits_image),
                                header=True)
    dither = fits_image[-6:-5]

    print('dither {}'.format(dither))

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

    print('fits_image {}'.format(fits_image))
    fits.writeto('{}/{}_t.fits'.format(prfs_d['fits_dir'], fits_image[:-5]),
                 data, header)

    return True


# if __name__=='__main__':
#     prfs_d = extract_settings_elvis()
#     logger = setting_logger(prfs_d)
#
#     change_times()
