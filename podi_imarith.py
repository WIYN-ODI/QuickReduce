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
import scipy
import scipy.stats

from podi_definitions import *

if __name__ == "__main__":

    # Read in the input parameters
    input_1 = sys.argv[1]
    op = sys.argv[2]
    input_2 = sys.argv[3]
    output = sys.argv[4]

    stdout_write("\nOpening input files ...")
    # Open both input fits files
    hdu_1 = pyfits.open(input_1)

    numeric_2 = None
    if (not os.path.isfile(input_2)):
        numeric_2 = float(input_2)
    else:
        hdu_2 = pyfits.open(input_2)

    stdout_write(" done!\n")

    rebin_fac = int(cmdline_arg_set_or_default("-bin", 1))
    
    # Now go though each extension and perform the operation
    for idx_img1 in range(1, len(hdu_1)):
        img1 = hdu_1[idx_img1]
        
        fppos1 = img1.header['EXTNAME']
        stdout_write("\rComputing extension %s (%2d of %2d) ..." % (img1.header['EXTNAME'], idx_img1, len(hdu_1)-1))
        if (numeric_2 != None):
            if (op == "+"):
                img1.data += numeric_2
            elif (op == "-"):
                img1.data -= numeric_2
            elif (op == "/"):
                img1.data /= numeric_2
            elif (op == "x"):
                img1.data *= numeric_2
            elif (op == "^"):
                img1.data = numpy.pow(img1.data, numeric_2)
            else:
                stdout_write("Unknown operation %s\n" % (op))

        else:
            for img2 in hdu_2[1:]:
                fppos2 = img2.header['EXTNAME']
                if (fppos2 == fppos1):
                    # This is the one

                    if (op == "+"):
                        img1.data += img2.data

                    elif (op == "-"):
                        img1.data -= img2.data

                    elif (op == "/"):
                        img1.data /= img2.data

                    elif (op == "x"):
                        img1.data *= img2.data

                    else:
                        stdout_write("Unknown operation %s\n" % (op))

        if (rebin_fac > 1):
            img1.data = rebin_image(img1.data, rebin_fac)

    # Now all the work is done, all final data is stored in img2, write results to new file                    
    stdout_write(" writing output ...")
    clobberfile(output)
    hdu_1.writeto(output)
    stdout_write(" done!\n\n")

