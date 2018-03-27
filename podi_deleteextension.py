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

Small tool to delete one or more extensions from a FITS file.

Usage:
-----------------

``podi_deleteextension.py input.fits output.fits extension1 extension2 ...``

"""

import sys
import os
import pyfits
import numpy
from podi_definitions import *

if __name__ == "__main__":

    inputfile = sys.argv[1]
    outputfile = sys.argv[2]

    hdulist = pyfits.open(inputfile)

    outlist = [pyfits.PrimaryHDU()]

    for i in range(1, len(hdulist)):

        try:
            extname = hdulist[i].header['EXTNAME']

            delete_this_one = False
            for ext_to_delete in sys.argv[3:]:
                if (extname == ext_to_delete):
                    delete_this_one = True
                    break

            if (delete_this_one):
                print("Deleting extension %s ..." % (extname))
                continue

        except:
            pass

        print("Keeping extension %s ..." % (extname))
        outlist.append(hdulist[i])



    print("writing")
    hdu_out = pyfits.HDUList(outlist)

    hdu_out.writeto(outputfile, clobber=True)
