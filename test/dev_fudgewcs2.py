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


import time
import astLib
import astLib.astWCS

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

    # Load the reference file
    reffile = sys.argv[1]

    hdulist = pyfits.open(reffile)

    # Check if the 3rd parameter is a filename or a number
    if (os.path.isfile(sys.argv[3])):

        # Start ds9
        ds9 = ds9.ds9(target="whirc", start="-scale scope global", wait=15)
        extension = sys.argv[2]

        #
        # This is the mode to determine the offsets.
        #

        ds9.set_np2arr(hdulist[extension].data.astype(numpy.float32))
        ds9.set("cmap grey")
        ds9.set("scale mode minmax")

        coords = ds9.get("imexam coordinate image")
        #print coords
        coords_ref = [float(coords.split()[0]), float(coords.split()[1])]

        # convert the pixel coords to Ra/Dec
        refwcs = astLib.astWCS.WCS(hdulist[extension].header, mode='pyfits')

        ra_dec = numpy.array(refwcs.pix2wcs(coords_ref[0], coords_ref[1]))
        #print ra_dec

        # Now we have the reference coordinates

        # Load the second frame to get the offset
        align_file = sys.argv[3]
        align_ext = sys.argv[4]

        hdulist2 = pyfits.open(align_file)
        ds9.set_np2arr(hdulist2[align_ext].data.astype(numpy.float32))
        ds9.set("cmap grey")
        ds9.set("scale mode minmax")
        coords2 = ds9.get("imexam coordinate image")
        #print coords2
        coords_align = [float(coords2.split()[0]), float(coords2.split()[1])]
        alignwcs = astLib.astWCS.WCS(hdulist2[align_ext].header, mode='pyfits')

        ra_dec_align = numpy.array(alignwcs.pix2wcs(coords_align[0], coords_align[1]))
        #print ra_dec_align

        d_radec = ra_dec - ra_dec_align
        #print d_radec

        print "%s %s %.6f %.6f %s.fudged.fits" % (sys.argv[0], align_file, d_radec[0], d_radec[1], align_file[:-5])

    else:

        # 
        # This is the mode where we apply the offsets
        # 
        d_ra = float(sys.argv[2])
        d_dec = float(sys.argv[3])

        for ext in hdulist:
            if ("CRVAL1" in ext.header):
                ext.header['CRVAL1'] += d_ra
            if ("CRVAL2" in ext.header):
                ext.header['CRVAL2'] += d_dec

        output_file = sys.argv[4]
        hdulist.writeto(output_file, clobber=True)


    sys.exit(0)

