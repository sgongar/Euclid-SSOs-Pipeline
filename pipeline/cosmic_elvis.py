# !/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for cosmic rays removal process

Versions:
- 0.1 Initial release
- 0.2 Move method

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from os import listdir

from astropy.io import fits
from astroscrappy import detect_cosmics
from multiprocessing import Process

from misc import extract_settings_elvis


class CosmicELViS:

    def __init__(self, logger):
        """

        :param logger:
        """
        self.prfs_d = extract_settings_elvis()
        self.logger = logger

        active_cosmic = []
        fits_files = listdir(self.prfs_d['fits_dir'])

        for cosmic_idx in range(0, len(fits_files),
                                self.prfs_d['cores_number']):
            try:
                cosmic_j = []
                for proc in range(0, self.prfs_d['cores_number'], 1):
                    idx = cosmic_idx + proc  # index
                    fits_file = fits_files[idx]
                    cosmic_p = Process(target=self.cosmic_thread,
                                       args=(fits_file,))
                    cosmic_j.append(cosmic_p)
                    cosmic_p.start()

                    active_cosmic = list([job.is_alive() for job in cosmic_j])
                while True in active_cosmic:
                    active_cosmic = list([job.is_alive() for job in cosmic_j])
                    pass
            except IndexError:
                print('Extraction finished')

        print('Extraction process of fits images finished')

    def cosmic_thread(self, fits_file):
        """

        :param fits_file:
        :return:
        """
        self.logger.debug('Cleans image {}'.format(fits_file))
        data_, header = fits.getdata('{}/{}'.format(self.prfs_d['fits_dir'],
                                                    fits_file),
                                     header=True)

        (cr_mask,
         cleaned_array) = detect_cosmics(data_, sigclip=4.5, sigfrac=0.3,
                                         objlim=5.0, gain=3.1, readnoise=4,
                                         satlevel=64535.0, pssl=0.0, niter=4,
                                         sepmed=True, cleantype='meanmask',
                                         fsmode='median', psfmodel='gauss',
                                         psffwhm=0.18, psfsize=7, psfk=None,
                                         psfbeta=4.765, verbose=False)

        fits.writeto('{}/{}'.format(self.prfs_d['fits_dir'], fits_file),
                     cleaned_array, header, clobber=True)
