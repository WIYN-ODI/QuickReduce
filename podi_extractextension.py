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
from podi_commandline import *


def extract_headers(inputfile, outputfile, extensions):


    hdulist = pyfits.open(inputfile)

    if (cmdline_arg_isset("-cleanprimary")):
        primhdu = pyfits.PrimaryHDU()
    else:
        primhdu = pyfits.PrimaryHDU(header=hdulist[0].header)

    outlist = [primhdu]

    for extension in extensions:

        print extension, extension=='PRIMARY'
        try:
            if (extension == 'PRIMARY'):
                ext = pyfits.ImageHDU(header=hdulist[0].header,
                                      data=hdulist[0].data,
                                      name="ORIG_PRIMARY")
            else:
                ext = hdulist[extension]
        except KeyError:
            try:
                iext = int(extension)
                ext = hdulist[iext]
            except:
                print "couldn't find extension",extension
                continue
        except:
            raise

        print "Adding ext", ext.name
        outlist.append(ext)

    if (cmdline_arg_isset("-singleext")):
        # Merge header from 1st and 2nd extension
        cards = outlist[1].header.cards
        for (keyword, value, comment) in cards:
            if (not keyword in outlist[1].header):
                outlist[1].header[keyword] = (value, comment)
        hdu_out = pyfits.HDUList([pyfits.PrimaryHDU(header=outlist[1].header,
                                                    data=outlist[1].data)])
        hdu_out.writeto(outputfile, clobber=True)
        return

    if (len(outlist) > 0):
        print "writing"
        hdu_out = pyfits.HDUList(outlist)
        hdu_out.info()
        hdu_out.writeto(outputfile, clobber=True)
    else:
        print "couldn't find any extension"



if __name__ == "__main__":

    inputfile = get_clean_cmdline()[1]
    outputfile = get_clean_cmdline()[-1]
    extensions = get_clean_cmdline()[2:-1]

    print extensions

    extract_headers(inputfile, outputfile, extensions)
