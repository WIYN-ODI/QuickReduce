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

sys.path.append("../")

from  podi_definitions import *
import podi_search_ipprefcat
import podi_wcs


import time

import ds9
ds9.ds9_xpans()


def fudge_wcs(input_file, ds9, output_file, target_coords):


#    ds9.set("file "+input_file)
    hdulist = pyfits.open(input_file)

  
#    ds9.set_pyfits(hdulist[extension])
#    ds9.set_np2arr(hdulist[extension].data)
    ds9.set("mosaicimage iraf "+input_file)
#    ds9.set("multiframe "+input_file)

    coords = ds9.get("imexam coordinate image")
    print coords
    coords_ref = [float(coords.split()[0]), float(coords.split()[1])]

    # 
    # With the user-selected pixel coordinates, compute the Ra/Dec value of this pixel
    # Since the frames contain distortion, make sure to use the podi_wcs routines.
    #
    wcs_poly = podi_wcs.header_to_polynomial(hdulist[extension].header)
    xy = numpy.array([[coords_ref[0], coords_ref[1]],
                      ])
    print xy.shape
    print xy[:,0:2]
    ra_dec = podi_wcs.wcs_pix2wcs(xy, wcs_poly)
    print ra_dec

    # Now change the CRVAL values in all frames so thie chosen pixel has ra/dec= 0,0
#    d_crval1 = hdulist[extension].header["CRVAL1"] - ra_dec[0,0]
#    d_crval2 = hdulist[extension].header["CRVAL2"] - ra_dec[0,1]
    d_crval1 = hdulist[extension].header["CRVAL1"] - ra_dec[0,0]
    d_crval2 = hdulist[extension].header["CRVAL2"] - ra_dec[0,1]

    print d_crval1, d_crval2

    for ext in range(len(hdulist)):
        if (not is_image_extension(hdulist[ext])):
            continue

#        hdulist[ext].header['CRVAL1'] -= d_crval1
#        hdulist[ext].header['CRVAL2'] -= d_crval2
        hdulist[ext].header['CRVAL1'] = d_crval1 + target_coords[0]
        hdulist[ext].header['CRVAL2'] = d_crval2 + target_coords[1]


    # output_file = sys.argv[3]
    print "writing output file",output_file
    hdulist.writeto(output_file, clobber=True)


if __name__ == "__main__":
    verbose=False

    extension = sys.argv[1]

    ds9 = ds9.ds9(target="whirc", start="-scale zscale -zoom 0.2 -scale scope global", wait=15)

    target_coords = [180., 0.]
    for input_file in sys.argv[2:]:
        output_file = input_file[:-5]+".fudge.fits"
        fudge_wcs(input_file, ds9, output_file, target_coords)
     



    sys.exit(0)

