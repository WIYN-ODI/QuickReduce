#!/usr/bin/env python


import sys
import numpy
import os
import pyfits
import datetime
import scipy
import scipy.stats
import math
import scipy.spatial
import itertools

import time
sys.path.append("../")
from podi_definitions import *





if __name__ == "__main__":
    verbose=False

    input = sys.argv[1]

    crval1 = float(sys.argv[2])
    crval2 = float(sys.argv[3])

    output = sys.argv[4]

    stdout_write("Loading ...")
    hdulist = pyfits.open(input)

    stdout_write(" fudging ...")
    for i in range(len(hdulist)):
        if (not is_image_extension(hdulist[i])):
            continue

        if ("CRVAL1" in hdulist[i].header):
            hdulist[i].header['CRVAL1'] = crval1

        if ("CRVAL2" in hdulist[i].header):
            hdulist[i].header['CRVAL2'] = crval2


    stdout_write(" writing ...")
    hdulist.writeto(output, clobber=True)
    stdout_write(" done!\n")


    sys.exit(0)

