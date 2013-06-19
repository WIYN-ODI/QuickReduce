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


def create_new_persistency_map(shape, write_fits=None):

    sy, sx = shape
    px, py = 8*sx, 8*sy

    # Create a primary header.
    # This only contains the MJD of this exposure
    primary_hdu = pyfits.PrimaryHDU()
    primary_hdu.header.update("MJD", 0.0, "MJD of exposure")
    
    primary_hdu.header.update("CELL_X", sx, "x-width of each cell")
    primary_hdu.header.update("CELL_Y", sy, "y-width of each cell")
    
    # Add primary header to HDU list
    hdulist = [primary_hdu]

    # Define some sizes to be used for displaying the frame as "Mosaic IRAF" in ds9
    iraf_gap = 100
    iraf_size_x, iraf_size_y = px+iraf_gap, py+iraf_gap

    stdout_write("Creating mask for OTA")
    for ota_x,ota_y in available_ota_coords:
        ota = ota_x * 10 + ota_y
        stdout_write(" %02d" % (ota))
        
        # Create new array with the full dimensions of the 8x8 cell array, 
        # with overscan still attached
        data = numpy.zeros(shape=(py,px), dtype=numpy.float32)

        # Create extension name
        ext_name = "OTA%02d.PERS" % (ota)

        # Create the ImageHDU
        imghdu = pyfits.ImageHDU(data=data)
        imghdu.update_ext_name(ext_name)

        # Add some additional info so we can display it in ds9:
        detsec = '[%d:%d,%d:%d]' % (ota_x*iraf_size_x, ota_x*iraf_size_x+px, ota_y*iraf_size_y, ota_y*iraf_size_y+py)
        imghdu.header.update("DETSEC", detsec, "Iraf mosaic area of the detector")

        detsize = '[1:%d,1:%d]' % (px, py)
        imghdu.header.update("DETSIZE", detsize, "Iraf total image pixels in full mosaic")

        # Add this OTA to the list of all OTAs in this map
        hdulist.append(imghdu)

    stdout_write(" done!\n")

    if (write_fits != None):
        stdout_write("Writing persistency map (%s) ..." % write_fits)
        fits_hdu = pyfits.HDUList(hdulist)
        clobberfile(write_fits)
        fits_hdu.writeto(write_fits, clobber=True)
        stdout_write(" done!\n")
        return
    else:
        stdout_write("Handing on results ...\n")

    return fits_hdu

if __name__ == "__main__":

    
    if (cmdline_arg_isset('-newmap')):
        inputfile = get_clean_cmdline()[1]
        outputfile = get_clean_cmdline()[2]
        hdulist = pyfits.open(inputfile)
        
        # If this flag is set, simply create a new persistency map
        shape = hdulist[1].data.shape
        create_new_persistency_map(shape, write_fits=outputfile)

        # Quit the program right here
        sys.exit(0)
        

    inputfile = sys.argv[1]
    hdulist = pyfits.open(inputfile)


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
