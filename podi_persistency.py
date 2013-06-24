#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy
import pywcs
from astLib import astWCS
import jdcal

from podi_definitions import *



def map_persistency_effects(hdulist, verbose=False):

    mask_thisframe_list = {}
    mask_timeseries_list = {}

    if (verbose): stdout_write("Creating persistency masks ...")
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

    if (verbose):
        stdout_write("\n   masked %d/%d pixels caused by %d saturated pixels in %d extensions\n" % (
                pixels_masked_out_thisframe, pixels_masked_out_timeseries, 
                saturated_pixels_total, extensions_with_saturated_pixels))

    return mask_thisframe_list, mask_timeseries_list




def mjd_to_time(mjd):

    year, month, day, time = jdcal.jd2gcal(2400000.5, mjd)

    hour = int(math.floor(time * 24.))
    x = time*24 - hour

    minute = int(math.floor(x * 60))
    x = x * 60 - minute

    second = x * 60

    return year, month, day, hour, minute, second


def get_timestamp_from_mjd(mjd):
    year, month, day, hour, minute, second = mjd_to_time(mjd)
    return "%04d%02d%02dT%02d%02d%02d" % (year, month, day, hour, minute, int(math.floor(second)))


def get_mjd_from_timestamp(timestamp):
    
    year, month, day = int(timestamp[0:4]), int(timestamp[4:6]), int(timestamp[6:8])
    hour, minute, second = int(timestamp[9:11]), int(timestamp[11:13]), int(timestamp[13:15])

    off, mjd1 = jdcal.gcal2jd(year, month, day)
    mjd2 = hour/24. + minute/1440. + second/86400.

    return mjd1 + mjd2


mjd_seconds = 1. / 86400.
def find_latest_persistency_map(directory, mjd, verbose=False):

    # Get a list of all files in the specified directory
    filelist = os.listdir(directory)

    min_delta_mjd = 1e9
    latest_map = None

    #
    # Now go over the files and find matching ones, and amonsgst the matching 
    # ones the one file with the closest time stamp
    # filename structure is: persistency_map_20121220T162345.fits
    #                                        |             |
    #                                      the usual timestamp
    #
    for file in filelist:
        #print file[:16]
        if (file[:16] != "persistency_map_" or file[31:] != ".fits"):
            continue

        # Extract timestamp and convert to MJD.
        timestamp = file[16:31]
        file_mjd = get_mjd_from_timestamp(timestamp)

        if (verbose): print file, file_mjd, mjd, "smaller:",(file_mjd<mjd)

        if (file_mjd >= mjd):
            # That's weird, this file shouldn't exist, but maybe it's just a 
            # re-run of the pipeline. Let's ignore it
            continue
        
        # Check if this is a closer match than the one we had before
        # Set 5 seconds as a minimum time requirement to prevent us from potentially 
        # finding the persistency map of this file from an earlier run.
        d_mjd = mjd - file_mjd
        if (d_mjd < min_delta_mjd and d_mjd > 5*mjd_seconds): 
            latest_map = file
            min_delta_mjd = d_mjd
            if (verbose): print "Found better match: %s (MJD=%.6f, or %d secs)" % (latest_map, file_mjd, min_delta_mjd*86400)

    if (latest_map == None):
        return None

    # Create full filename and return
    fullpath = "%s/%s" % (directory, latest_map)
    print "Using",fullpath,"as persistency map"
    return fullpath

def persistency_map_filename(directory, mjd):
    
    # First convert MJD to timestamp
    timestamp = get_timestamp_from_mjd(mjd)

    # And then create and return filename
    filename = "%s/persistency_map_%s.fits" % (directory, timestamp)
    
    return filename


def add_mask_to_map(mask, mjd, map_in):

    # Make a copy of the input frame
    map_out = map_in.copy()

    # Compute the width and height of one cell
    dx, dy = map_in.shape[1]/8, map_in.shape[0]/8

    this_map = numpy.zeros(map_in.shape)

    #print mask
    #print map_in.shape

    for cell in mask:
        # print "in add_mask_to_map",cell, mask[cell].shape

        x,y = int(cell[2]), int(cell[3])
        
        # Compute the bottom left and top right pixel coordinates of this cell in the existing map
        bx, tx, by, ty = cell2ota__get_target_region(x, y)

        # In this small cell region, apply mask and update the MJD with the given timestamp
        mask_datasec = cell2ota__extract_data_from_cell(mask[cell])

        map_out[by:ty,bx:tx][mask_datasec] = mjd

    return map_out
        
def apply_mask_to_data(mask, data):

    out = data
    out[mask] = numpy.NaN

    return out

def get_correction(persistency_map, cell_position, mjd):

    # Compute how big each subframe is
    dx, dy = persistency_map.shape[1]/8, persistency_map.shape[0]/8

    cell_x, cell_y = cell_position

    # From this and the cell-coordinates, determine the 
    # bottom left and top right pixel coordinates
    bx, tx, by, ty = cell2ota__get_target_region(cell_x, cell_y)

    # Now extract frame and compute time-difference 
    # between then and now, and convert delta_MJDs into seconds
    d_mjd = (persistency_map[by:ty,bx:tx] - mjd) * 86400.

    # Add some clever function here...
    # Only use correction if they are within 10 minutes of the frame
    invalid_range_dmjd = (d_mjd < -600) | (d_mjd >= 0)
    correction = 20. * numpy.exp(d_mjd / 125.)
    correction[invalid_range_dmjd] = 0

    return correction

    
def subtract_persistency(persistency_map, image_hdu):

    return


def create_new_persistency_map(shape=None, write_fits=None):

    if (shape == None):
        sx, sy = 480, 494
        px, py = 4096, 4096
    else:
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
        outputfile = "persistency_map_00000000T000000.fits" #get_clean_cmdline()[2]
        
        # If this flag is set, simply create a new persistency map
        create_new_persistency_map(None, write_fits=outputfile)

        # Quit the program right here
        sys.exit(0)
        
    if (cmdline_arg_isset('-findmap')):
        directory = get_clean_cmdline()[1]
        mjd = float(get_clean_cmdline()[2])
        find_latest_persistency_map(directory, mjd, verbose=True)
        sys.exit(0)

    inputfile = sys.argv[1]
    hdulist = pyfits.open(inputfile)
    persistency_map_in = sys.argv[2]

    outputfile = sys.argv[3]
    persistency_map_out = sys.argv[4]

    hdulist = pyfits.open(inputfile)
    mjd = hdulist[0].header['MJD-OBS']


    mask_thisframe, mask_timeseries = map_persistency_effects(hdulist)

    persistency_hdu = pyfits.open(persistency_map_in)
    map_in = persistency_hdu[1].data

    map_out = add_mask_to_map(mask_timeseries, mjd, map_in)
    persistency_hdu[1].data = map_out
    persistency_hdu.writeto(persistency_map_out, clobber=True)

    for ext in range(0, len(hdulist)):
        if (str(type(hdulist[ext])) != "<class 'pyfits.hdu.image.ImageHDU'>"):
            continue

        extname = hdulist[ext].header['EXTNAME']

        if (extname in mask_thisframe):
            hdulist[ext].data[mask_thisframe[extname]] = 100
            # data[mask_time] = mjd
            # data[saturated] = 0

    #hdulist.writeto("persistency.fits", clobber=True)
