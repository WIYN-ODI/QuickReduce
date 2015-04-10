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

gain_correct_frames = False
from podi_definitions import *


do_not_copy_headers = ("NAXIS", "NAXIS1", "NAXIS2")

verbose = False

if __name__ == "__main__":

    outputfile = sys.argv[1]
    stdout_write("### podi_pipeline2mef\n")

    ota_list = []

    primhdu = pyfits.PrimaryHDU()
    ota_list.append(primhdu)

    for filename in sys.argv[2:]:

        tmpfile = None
        if (filename.endswith(".fz")):
            _, fn = os.path.split(filename)
            tmpfile = "/N/dc2/scratch/odiuser/aucap_tiles/%s" % (fn[:-3])
            os.system("funpack -O %s %s" % (tmpfile, filename))
            hdulist = pyfits.open(tmpfile)
        else:
            hdulist = pyfits.open(filename)
        print "Adding in",filename
        ota_list.append(pyfits.ImageHDU(header=hdulist[0].header,
                                        data=hdulist[0].data))
        # translate FPPOS into OTA
        fppos = 'xy00' if 'FPPOS' not in hdulist[0].header else hdulist[0].header['FPPOS']
        ota = int(fppos[2:4])
        ota_list[-1].name = 'OTA%02d.SCI' % (ota)
        if (not tmpfile == None):
            os.remove(tmpfile)


    stdout_write("\n   writing output file %s ..." % (outputfile))

    hdulist = pyfits.HDUList(ota_list)
    if (os.path.isfile(outputfile)):
	os.remove(outputfile)
    hdulist.writeto(outputfile, clobber=True)

    stdout_write(" done!\n")
    stdout_write("### complete!\n")
