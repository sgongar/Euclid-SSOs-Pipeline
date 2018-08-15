# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 16:44:28 2014

@author: lconversi
"""

import glob
from subprocess import call
from postproc import tileCCD
from postproc import tileFPA
from support import logger as lg

for dither in range(1, 5):
    """
    for xCCD in range(6):
        for yCCD in range(6):
            print
            print '============================================'
            print '  Merging quadrants of xCCD #%d and yCCD #%d' % (xCCD, yCCD)
            print '============================================'
            #
            files = glob.glob('Q*_x%d_y%d_d1.fits' % (xCCD, yCCD))
            #
            if len(files) > 0:
                output = 'CCD_x%d_y%d_d1.fits' % (xCCD, yCCD)
                inputs = dict(files=files, ext=1, output=output)
                #
                tileCCD.tileCCD(inputs,
                                lg.setUpLogger('tileCCDs.log')).runAll()
    """
    print
    print '=============================================='
    print '  Dither n.', dither
    print '  Merging all CCDs from same dither position'
    print '=============================================='
    #
    files = glob.glob('CCD*_d%d.fits' % (dither))
    print len(files)
    print files
    from time import sleep
    sleep(10)
    output = 'FPA_Dith%d.fits' % (dither)
    inputs = dict(files=files, ext=0, output=output)
    tileFPA.tileFPA(inputs, lg.setUpLogger('tileFPA.log')).runAll()

# Deleting all single CCD and log files
call('rm *log', shell=True)
