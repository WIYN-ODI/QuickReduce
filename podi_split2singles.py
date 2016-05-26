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
from podi_definitions import is_image_extension

if __name__ == "__main__":


    for filename in sys.argv[1:]:

        #filename = sys.argv[1]
        print filename

        hdulist = pyfits.open(filename)

        for extension in range(1, len(hdulist)):

            if (not is_image_extension(hdulist[extension])):
                continue

            extname = hdulist[extension].header['EXTNAME']
            # if (extname[0:3] != "OTA" or extname[-3:] != "SCI"):
            #     continue

            outfits = filename[:-4]+extname+".fits"

            primhdu = pyfits.PrimaryHDU(header=hdulist[extension].header,
                                        data=hdulist[extension].data)

            # copy all primary header entries into the new primary HDU
            for (key, value) in hdulist[0].header.iteritems():
                # print key, value
                primhdu.header[key] = value

            out_hdulist = pyfits.HDUList([primhdu])
            out_hdulist.writeto(outfits, clobber=True)
