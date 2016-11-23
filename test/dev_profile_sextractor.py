#!/usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

"""

Convert the multi-file output of the NOAO-supplied, official AuCaP pipeline into
a single  multi-extension FITS file that more closely resembles the output of
podi_collectcells.

To run
------

``podi_aucap2mef.py output.fits aucap1.fits aucap2.fits ...``


"""

import sys
import os
import pyfits
import numpy
import scipy
import time

gain_correct_frames = False
from podi_definitions import *


do_not_copy_headers = ("NAXIS", "NAXIS1", "NAXIS2")

verbose = False


if __name__ == "__main__":

    fitsfile = sys.argv[1]

    dummy_log = open("%s.performance.log" % (fitsfile[:-5]), "a")
    for detsigma in range(1,20):
        for detarea in range(1,10,1):

            if (detsigma * math.sqrt(detarea) < 2):
                continue

            catname = "%s.cat.sigma%02d.detarea%02d.%s" % (
                fitsfile[:-5], detsigma, detarea, time.strftime("%Y%m%dT%H%M%S")
            )

            opts = {
                'catname': catname,
                'filename': fitsfile,
                'detarea': detarea,
                'detsigma': detsigma,
            }

            # sex-command
            sex = "sex -c /work/podi_devel/config/wcsfix.sex \
            -CATALOG_NAME %(catname)s -DETECT_MINAREA %(detarea)02d \
            -DETECT_THRESH %(detsigma)02d -ANALYSIS_THRESH %(detsigma)02d \
            %(filename)s -VERBOSE_TYPE NORMAL " % opts

            # print sex
            print """\
#################
#
#  det-sigma=%(detsigma)d    det-area=%(detarea)d
#
#################""" % opts


            start_clock = time.clock()
            start_time = time.time()
            os.system(sex)
            end_clock = time.clock()
            end_time = time.time()

            # Run Sextractor

            # Load file
            data = numpy.loadtxt(catname)
            number_sources = data.shape[0]

            good_sources = data[:,7] == 0
            n_noflags = numpy.sum(good_sources)

            elapsed_clock = (end_clock - start_clock)
            elapsed_time = (end_time - start_time)

            print detsigma, detarea, number_sources, n_noflags, elapsed_time, elapsed_clock, catname
            print >>dummy_log, detsigma, detarea, number_sources, n_noflags, elapsed_time, elapsed_clock, catname
            dummy_log.flush()

    dummy_log.close()

            
                
