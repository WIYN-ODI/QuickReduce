#!/usr/local/bin/python
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
import time
import podi_persistency

if __name__ == "__main__":

#    for i in range(len(sys.argv)):
#        print i, sys.argv[i]

    reference_file = sys.argv[1]
    extension = sys.argv[2]

    output_fits = sys.argv[3]
    output_txt = sys.argv[4]

    # Create mask from saturated pixels in reference file.

    hdu = pyfits.open(reference_file)
    data = hdu[extension].data

    #filetime = time.strptime(reference_file[0:15],"%Y%m%dT%H%M%S")

    ref_mjd = podi_persistency.get_mjd_from_timestamp(reference_file[0:15])
    #hdu[0].header['MJD-OBS']

    mask = (data >= 64500) #65535)

    hdu.close()
    del  hdu

    data[:] = 0
    data[mask] = 1

    hdulist_out = [pyfits.PrimaryHDU(data=data)]

    txtout = open(output_txt, "w")

    for file in sys.argv[5:]:
        print >>txtout, "#", file, extension

    print >>txtout, ref_mjd, 0, 65535, 65535

    # Now go through each of the files and get the median level of the formerly saturated pixels
    for file in sys.argv[5:]:
        print "Measuring file",file
        this_hdu = pyfits.open(file)
        this_data = this_hdu[extension].data
        this_mjd = mjd = podi_persistency.get_mjd_from_timestamp(file[0:15])
#this_hdu[0].header['MJD-OBS']

        saturated = this_data[mask]
        #print saturated[0]

        mean = numpy.mean(saturated)
        median = numpy.median(saturated)

        hdudump = pyfits.ImageHDU(data=saturated)
        hdudump.header.update("MJD", this_mjd)
        #hdudump.header.update("MEAN", mean)
        #hdudump.header.update("MEDIAN", median)
        
        hdulist_out.append(hdudump)

        print >>txtout, this_mjd, this_mjd-ref_mjd, mean, median

        del this_data
        this_hdu.close()
        del this_hdu

    hduout = pyfits.HDUList(hdulist_out)
    hduout.writeto(output_fits, clobber=True)
