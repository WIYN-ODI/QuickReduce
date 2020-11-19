#! /usr/bin/env python3
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
import astropy.io.fits as pyfits
import numpy
from podi_definitions import *
from podi_commandline import *

if __name__ == "__main__":

    extension_name = get_clean_cmdline()[1]
    out_prefix = get_clean_cmdline()[2]


    for inputfile in get_clean_cmdline()[3:]:

        _,basename = os.path.split(inputfile)

        print("Extracting from",inputfile,)
        hdulist = pyfits.open(inputfile)

        if (cmdline_arg_isset("-cleanprimary")):
            primhdu = pyfits.PrimaryHDU()
        else:
            primhdu = pyfits.PrimaryHDU(header=hdulist[0].header)

        outlist = [primhdu]

        try:
            outlist.append(hdulist[extension_name])
        except:
            print("--> extension not found")
            continue

        
        if (len(outlist) > 1):
            out_filename = "%s.%s.fits" % (basename[:-5], out_prefix)
            print(" --> writing to",out_filename,)
            hdu_out = pyfits.HDUList(outlist)
            hdu_out.writeto(out_filename, overwrite=True)
            print("done!")

        else:
            print("couldn't find any extension")
