#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy
import pywcs
from astLib import astWCS

from podi_definitions import *



def map_persistency_effects(hdulist):

    mask_thisframe_list = {}
    mask_timeseries_list = {}

    stdout_write("Creating persistency masks ...")
    saturated_pixels_total = 0
    extensions_with_saturated_pixels = 0
    pixels_masked_out_thisframe = 0
    pixels_masked_out_timeseries = 0

    for ext in range(len(hdulist)):
        if (str(type(hdulist[ext])) != "<class 'pyfits.hdu.image.ImageHDU'>"):
            continue

        extname = hdulist[ext].header['EXTNAME']
        # stdout_write("Working on extension %s (%d)\n" % (extname, ext))

        data = hdulist[ext].data
        saturated = (data >= 65535)

        number_saturated_pixels = numpy.sum(saturated)
        if (number_saturated_pixels <= 0):
            continue

        saturated_pixels_total += number_saturated_pixels
        extensions_with_saturated_pixels += 1

        rows, cols = numpy.indices(data.shape)

        mask_thisframe = numpy.zeros(shape=data.shape)
        mask_thisframe = mask_thisframe > 1
        mask_time      = numpy.zeros(shape=data.shape)
        mask_time      = mask_time > 1
        #mask_time.fill(False)

        saturated_rows = rows[saturated]
        saturated_cols = cols[saturated]

        # print "# saturated_pixels =",saturated_rows.shape[0], saturated_cols.shape[0]

        #print "overall:",mask_time.shape, mask_thisframe.shape

        for i in range(saturated_rows.shape[0]):
            mask_up = (cols == saturated_cols[i]) & (rows >= saturated_rows[i])
            mask_down = (cols == saturated_cols[i]) & (rows <= saturated_rows[i])
            #print "this:",mask_up.shape, mask_down.shape

            mask_thisframe = (mask_thisframe) | (mask_up)
            mask_time      = (mask_time)      | (mask_down)

        mask_thisframe_list[extname] = mask_thisframe
        mask_timeseries_list[extname] = mask_time

        pixels_masked_out_thisframe += numpy.sum(mask_thisframe)
        pixels_masked_out_timeseries += numpy.sum(mask_time)
        #data[mask_thisframe] = 100
        #data[mask_time] = mjd
        #data[saturated] = 0

    stdout_write("\n   masked %d/%d pixels caused by %d saturated pixels in %d extensions\n" % (
            pixels_masked_out_thisframe, pixels_masked_out_timeseries, 
            saturated_pixels_total, extensions_with_saturated_pixels))

    return mask_thisframe_list, mask_timeseries_list


if __name__ == "__main__":

    
    inputfile = sys.argv[1]
    persistency_map_file = sys.argv[2]

    hdulist = pyfits.open(inputfile)
    mjd = hdulist[0].header['MJD-OBS']

    if (os.path.isfile(persistency_map_file)):
        persistency_hdu = pyfits.open(persistency_map_file)

    mask_thisframe, mask_timeseries = map_persistency_effects(hdulist)

    for ext in range(0, len(hdulist)):
        if (str(type(hdulist[ext])) != "<class 'pyfits.hdu.image.ImageHDU'>"):
            continue

        extname = hdulist[ext].header['EXTNAME']

        if (extname in mask_thisframe):
            hdulist[ext].data[mask_thisframe[extname]] = 100
            # data[mask_time] = mjd
            # data[saturated] = 0

    hdulist.writeto("persistency.fits", clobber=True)
