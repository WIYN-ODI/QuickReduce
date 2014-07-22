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

import sys
import os
import pyfits
import numpy
from podi_definitions import *

if __name__ == "__main__":

    for inputfile in sys.argv[1:]:
        
        hdulist = pyfits.open(inputfile)
        outlist = []

        for ext in hdulist:
            #print type(ext)
            if (is_image_extension(ext)):
                outlist.append(ext)

        print "writing",inputfile
        hdu_out = pyfits.HDUList(outlist)
        hdu_out.writeto(inputfile+".notab.fits", clobber=True)
