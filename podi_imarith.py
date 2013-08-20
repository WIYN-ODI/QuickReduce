#! /usr/bin/env python

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
    hdu_2 = pyfits.open(input_2)
    stdout_write(" done!\n")

    rebin_fac = int(cmdline_arg_set_or_default("-bin", 1))
    
    # Now go though each extension and perform the operation
    for idx_img1 in range(1, len(hdu_1)):
        img1 = hdu_1[idx_img1]
        
        fppos1 = img1.header['EXTNAME']
        stdout_write("\rComputing extension %s (%2d of %2d) ..." % (img1.header['EXTNAME'], idx_img1, len(hdu_1)-1))
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
                    img1 *= img2.data

                else:
                    stdout_write("Unkwnon operation %s\n" % (op))

        if (rebin_fac > 1):
            img1.data = rebin_image(img1.data, rebin_fac)

    # Now all the work is done, all final data is stored in img2, write results to new file                    
    stdout_write(" writing output ...")
    clobberfile(output)
    hdu_1.writeto(output)
    stdout_write(" done!\n\n")

