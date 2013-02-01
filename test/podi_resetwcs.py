#! /usr/bin/env python

import sys
import os
import pyfits
import numpy

from podi_definitions import *

if __name__ == "__main__":

    # Read in the input parameters
    input = sys.argv[1]
    output = sys.argv[2]
    
    # Open both input fits files
    stdout_write("\nOpening input file ...")
    hdu = pyfits.open(input)

    stdout_write(" working ")

    for img1 in hdu[1:]:
        stdout_write(".")
        img1.header['CRVAL1'] = 0
        img1.header['CRVAL2'] = 0

    # Now all the work is done, all final data is stored in img2, write results to new file                    
    stdout_write(" writing results ... ")
    clobberfile(output)
    hdu.writeto(output)
    stdout_write(" done!\n\n")

