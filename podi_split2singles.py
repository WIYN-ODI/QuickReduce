#! /usr/bin/env python

#
# (c) Ralf Kotulla for WIYN/pODI
#

import sys
import os
import pyfits
import numpy

if __name__ == "__main__":


    for filename in sys.argv[1:]:

        #filename = sys.argv[1]
        print filename

        hdulist = pyfits.open(filename)

        for extension in range(1, len(hdulist)):
            
            extname = hdulist[extension].header['EXTNAME']
            if (extname[0:3] != "OTA" or extname[-3:] != "SCI"):
                continue

            outfits = filename[:-4]+extname+".fits"

            primhdu = pyfits.PrimaryHDU(header=hdulist[extension].header,
                                        data=hdulist[extension].data)
            out_hdulist = pyfits.HDUList([primhdu])
            out_hdulist.writeto(outfits, clobber=True)
