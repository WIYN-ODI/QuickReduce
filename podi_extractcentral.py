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

Extract only the central 3x3 extensions frm the input FITS file and write to
output FITS file.

How to run:
-----------

``podi_extractcentral.py /some/output/directory file1.fits file2.fits ...``

"""


import sys
import os
import astropy.io.fits as pyfits
import numpy
from podi_definitions import *

if __name__ == "__main__":

    output_dir = sys.argv[1]

    for filename in sys.argv[2:]:

        stdout_write("Reading %s ..." % (filename))

        hdulist = pyfits.open(filename)
        out_list = [hdulist[0]]

        outfits = "%s/%s.central3x3.fits" % (output_dir, filename[:-5])
        if (os.path.isfile(outfits)):
            print(" file exists, skipping")
            continue

        for extension in range(1, len(hdulist)):
            
            extname = hdulist[extension].header['EXTNAME']
            if (extname[0:3] != "OTA" or extname[-3:] != "SCI"):
                continue

            ota_x = int(extname[3])
            ota_y = int(extname[4])
            if (ota_x >= 2 and ota_x <= 4 and ota_y >= 2 and ota_y <= 4):
                out_list.append(hdulist[extension])

        out_hdulist = pyfits.HDUList(out_list)
        stdout_write(" writing %s ..." % (outfits))
        out_hdulist.writeto(outfits, overwrite=True)
        #print "-->",outfits
        stdout_write(" done!\n")
