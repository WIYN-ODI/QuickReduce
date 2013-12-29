#! /usr/bin/env python
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

Extract a single extension from the input FITS file.

Usage:
---------

``podi_extractextension input.fits EXT.NAME output.fits``

"""

import sys
import os
import pyfits
import numpy
from podi_definitions import *

if __name__ == "__main__":

    inputfile = sys.argv[1]
    extension = sys.argv[2]
    outputfile = sys.argv[3]

    hdulist = pyfits.open(inputfile)

    primhdu = pyfits.PrimaryHDU()
    if (cmdline_arg_isset("-keepprimary")):
        primhdu = pyfits.PrimaryHDU(header=hdulist[0].header)

    outlist = [primhdu]
    
    try:
        ext = hdulist[extension]
    except:
        print "couldn't find extension"
        sys.exit(0)

    outlist.append(ext)

    print "writing"
    hdu_out = pyfits.HDUList(outlist)
    hdu_out.writeto(outputfile, clobber=True)
