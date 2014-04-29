#!/usr/bin/env python

#
# Copyright (C) 2014, Ralf Kotulla
#                     kotulla@uwm.edu
#
# All rights reserved
#

import os, sys
d,_=os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append("%s/../"%d)

import podi_asteroids
import podi_swarpstack
from podi_commandline import *
from podi_definitions import *
import podi_logging
import logging
import astropy.io.votable
import math
import numpy
import ephem
import time
import scipy.stats

if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    logger = logging.getLogger("AstroidField")

    params = podi_swarpstack.read_swarp_params()

    if (cmdline_arg_isset('-filelist')):
        fl = open(get_cmdline_arg('-filelist'), 'r')
        lines = fl.readlines()
        inputlist = []
        for line in lines:
            filename = line.split()[0]
            if (os.path.isfile(filename)):
                inputlist.append(filename)
        print "filelist read from input file:"
        print "\n".join(inputlist)
    else:
        inputlist = get_clean_cmdline()[2:]

    target_name = get_clean_cmdline()[1]

    smaller_region = 5 # arcmin on a side

    print params
    print options

    # First, stack all frames without any non-sidereal correction
    outputfile = "%s__reference.fits" % (target_name)

    logger.info("Creating the averaged reference stack")
    returned = podi_swarpstack.swarpstack(outputfile, inputlist, params, options, keep_intermediates=True)
    if (returned == None):
        logger.error("something went wrong while creating the reference stack")
        
    else:

        modified_files, single_prepared_files, bgsub_files, unique_singledir = returned
        print single_prepared_files

        # Open the reference frame
        # This we will need for each single frame
        logger.info("opening reference file")
        hdu_ref = astropy.io.fits.open(outputfile)
        ref_frame = hdu_ref[0].data

        for sglfile in single_prepared_files:

            logger.info("Creating difference image from %s" % (sglfile))
            hdu_sgl = astropy.io.fits.open(sglfile)

            hdu_sgl[0].data -= ref_frame
            
            # Strip the fits from the filename and append a ".diff.fits"
            _, base = os.path.split(sglfile)
            sgl_diff_filename = "%s___%s.diff.fits" % (target_name, base)

            logger.debug("Also copying the weight image")
            sgl_weight = sglfile[:-5]+".weight.fits"
            hdu_weight = astropy.io.fits.open(sgl_weight)

            # Mask out the entire region that has 0 weight
            hdu_sgl[0].data[hdu_weight[0].data <= 0] = numpy.NaN

            logger.debug("Writing difference image to %s" % (sgl_diff_filename))
            clobberfile(sgl_diff_filename)
            hdu_sgl.writeto(sgl_diff_filename)
            
            diff_weight = "%s___%s.diff.weight.fits" % (target_name, base)
            clobberfile(diff_weight)
            hdu_weight.writeto(diff_weight, clobber=True)

        logger.info("all done!")


    podi_logging.shutdown_logging(options)
