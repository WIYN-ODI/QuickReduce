#! /usr/bin/env python

#
# (c) Ralf Kotulla for WIYN/pODI
#

import sys
import os
import pyfits
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
            print " file exists, skipping"
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
        out_hdulist.writeto(outfits, clobber=True)
        #print "-->",outfits
        stdout_write(" done!\n")


