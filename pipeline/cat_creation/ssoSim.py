# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 12:21:59 2016

@author: lconversi
"""

# CCD Image Simulation
from simulator.simulator import *
# from subprocess import call
from os import remove, getcwd

opts, args = processArgs()

# General options
path = '/pcdisk/holly/sgongora/Documents/Euclid/Euclid-tests/method_creation/'
opts.configfile = path + 'sso.config'
opts.section = 'Default'

for xCCD in range(6):
    for yCCD in range(6):
        for quadrant in range(4):
            for dither in range(1, 2):
                print getcwd()
                print
                print '==============================================='
                print
                print '  Simulating quadrant #%d, xCCD #%d, yCCD #%d and dither #%d' %(quadrant, xCCD, yCCD, dither)
                print
                print '==============================================='
                print
                #
                opts.sourcelist = path + 'Cat_d%d.dat' % (dither)
                print opts.sourcelist
                opts.quadrant = quadrant
                opts.xCCD = xCCD
                opts.yCCD = yCCD
                opts.dither = dither
                #
                # Run the simulator
                simulate = VISsimulator(opts)
                simulate.simulate()

                nobackgroundfile = "nobackgroundQ%d_0%d_0%d.fits" % (quadrant, xCCD, yCCD)
                nonoisefile = "nonoiseQ%d_0%d_0%d.fits" % (quadrant, xCCD, yCCD)
                readoutnoisefile = "readoutnoiseQ%d_0%d_0%d.fits" % (quadrant, xCCD, yCCD)

                print nobackgroundfile, nonoisefile, readoutnoisefile
                remove(getcwd() + '/' + nobackgroundfile)
                remove(getcwd() + '/' + nonoisefile)
                remove(getcwd() + '/' + readoutnoisefile)
