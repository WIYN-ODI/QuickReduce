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
import pywcs
from astLib import astWCS
import jdcal

import time
import multiprocessing
import Queue

from podi_definitions import *
import podi_sitesetup as sitesetup



if (sitesetup.number_cpus == "auto"):
    try:
        number_cpus = multiprocessing.cpu_count()
        print "Yippie, found %d CPUs to use in parallel!" % (number_cpus)
        if (number_cpus > sitesetup.max_cpu_count and sitesetup.max_cpu_count > 1):
            number_cpus = sitesetup.max_cpu_count
            print "... but using only %d of them!" % (number_cpus)
    except:
        pass
else:
    number_cpus = sitesetup.number_cpus


def mp_create_saturation_catalog(queue_in, queue_ret, verbose=False):
    """
    This is a small helper routine for the process of creating the saturation catalogs.
    It reads filenames from job queue, creates the arrays of pixel coordinates, and 
    posts the results to a return queue. Actually creating the fits tables is then 
    handled by the main process.
    """

    while (True):
        filename = queue_in.get()

        if (filename == None):
            queue_in.task_done()
            return

        cat_name = create_saturation_catalog_ota(filename, None, verbose=verbose, return_numpy_catalog=True)

        queue_ret.put( cat_name )

        queue_in.task_done()

    return



def create_saturation_catalog(filename, output_dir, verbose=True, mp=False):

    if (os.path.isfile(filename)):
        # This is one of the OTA fits files
        # extract the necessary information to generate the 
        # names of all the other filenames
        hdulist = pyfits.open(filename)
        basename = hdulist[0].header['FILENAME'][:18]
        hdulist.close()

        # Split the input filename to extract the directory part
        directory, dummy = os.path.split(filename)

    elif (os.path.isdir(filename)):
        # As a safety precaution, if the first parameter is the directory containing 
        # the files, extract just the ID string to be used for this script
        if (filename[-1] == "/"):
            filename = filename[:-1]

        basedir, basename = os.path.split(filename)
        directory = filename

    output_filename = "%s/%s.saturated.fits" % (output_dir, basename)
    stdout_write("%s --> %s ...\n" % (filename, output_filename))

    # Setup parallel processing
    queue        = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.JoinableQueue()
    #return_queue = multiprocessing.Queue()
        
    number_jobs_queued = 0
    first_fits_file = None
    ota_list = []

    for (ota_x, ota_y) in available_ota_coords:
        ota = ota_x * 10 + ota_y

        filename = "%s/%s.%02d.fits" % (directory, basename, ota)
        if (not os.path.isfile(filename)):
            filename = "%s/%s.%02d.fits.fz" % (directory, basename, ota)
            if (not os.path.isfile(filename)):
                continue


        queue.put( (filename) )
        number_jobs_queued += 1

        # Remember the very first fits file we find. This will serve as the primary HDU
        if (first_fits_file == None): first_fits_file = filename
       
    # Now start all the workers
    processes = []
    for i in range(number_cpus):
        p = multiprocessing.Process(target=mp_create_saturation_catalog, args=(queue, return_queue, False))
        p.start()
        processes.append(p)
        time.sleep(0.01)

    # Tell all workers to shut down when no more data is left to work on
    for i in range(len(processes)):
        if (verbose): stdout_write("Sending quit command!\n")
        queue.put( (None) )

    # Create a primary HDU from the first found fits-file
    firsthdu = pyfits.open(first_fits_file)
    ota_list.append(pyfits.PrimaryHDU(header=firsthdu[0].header))
    firsthdu.close()
    firsthdu = None

    for i in range(number_jobs_queued):
        if (verbose): print "reading return ",i

        cat_name = return_queue.get()
        if (cat_name != None):
            final_cat, extension_name = cat_name

            columns = [\
                pyfits.Column(name='CELL_X', format='I', array=final_cat[:, 0]),
                pyfits.Column(name='CELL_Y', format='I', array=final_cat[:, 1]),
                pyfits.Column(name='X',      format='I', array=final_cat[:, 2]),
                pyfits.Column(name='Y',      format='I', array=final_cat[:, 3])
                ]
            # Create the table extension
            coldefs = pyfits.ColDefs(columns)
            tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')
            tbhdu.update_ext_name(extension_name, comment="catalog of saturated pixels")

            ota_list.append(tbhdu)
            
        return_queue.task_done()
            
    hdulist = pyfits.HDUList(ota_list)
    output_filename = "%s/%s.saturated.fits" % (output_dir, basename)
    clobberfile(output_filename)
    hdulist.writeto(output_filename, clobber=True)

    return



def create_saturation_catalog_ota(filename, output_dir, verbose=True, return_numpy_catalog=False):

    # Open filename
    if (verbose):
        stdout_write("Creating catalog of saturated pixels\n")
        stdout_write("Input filename: %s\n" % (filename))

    hdulist = pyfits.open(filename)

    mjd = hdulist[0].header['MJD-OBS']
    obsid = hdulist[0].header['OBSID']
    ota = int(hdulist[0].header['FPPOS'][2:4])
    datatype = hdulist[0].header['FILENAME'][0]

    full_coords = numpy.zeros(shape=(0,4)) #, dtype=numpy.int16)
    saturated_pixels_total = 0

    for ext in range(1, len(hdulist)):
        if (type(hdulist[ext]) != pyfits.hdu.image.ImageHDU):
            continue

        # Find all saturated pixels (values >= 65K)
        data = hdulist[ext].data
        saturated = (data >= 65535)

        # Skip this cell if no pixels are saturated
        number_saturated_pixels = numpy.sum(saturated)
        if (number_saturated_pixels <= 0):
            continue

        saturated_pixels_total += number_saturated_pixels
        
        wn_cellx = hdulist[ext].header['WN_CELLX']
        wn_celly = hdulist[ext].header['WN_CELLY']

        if (verbose): print "number of saturated pixels in cell %d,%d: %d" % (wn_cellx, wn_celly, number_saturated_pixels)

        # Do some book-keeping preparing for the masking
        rows, cols = numpy.indices(data.shape)

        saturated_rows = rows[saturated]
        saturated_cols = cols[saturated]

        #print saturated_rows.shape, saturated_cols.shape

        coordinates = numpy.zeros(shape=(number_saturated_pixels,4))
        coordinates[:,0] = wn_cellx
        coordinates[:,1] = wn_celly
        coordinates[:,2] = saturated_cols[:]
        coordinates[:,3] = saturated_rows[:]

        full_coords = numpy.append(full_coords, coordinates, axis=0) #coordinates if full_coords == None else 

    final_cat = numpy.array(full_coords, dtype=numpy.dtype('int16'))

    if (saturated_pixels_total <= 0):
        return None

    # Now define the columns for the table
    columns = [\
        pyfits.Column(name='CELL_X', format='I', array=final_cat[:, 0]),
        pyfits.Column(name='CELL_Y', format='I', array=final_cat[:, 1]),
        pyfits.Column(name='X',      format='I', array=final_cat[:, 2]),
        pyfits.Column(name='Y',      format='I', array=final_cat[:, 3])
        ]
    # Create the table extension
    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')
    extension_name = "OTA%02d.SATPIX" % (ota)
    tbhdu.update_ext_name(extension_name, comment="catalog of saturated pixels")

    if (return_numpy_catalog):
        return final_cat, extension_name

    # Also copy the primary header into the new catalog
    primhdu = pyfits.PrimaryHDU(header=hdulist[0].header)

    # Create a HDUList for output
    out_hdulist = pyfits.HDUList([primhdu, tbhdu])
    
    # And create the output file
    output_filename = "%s/%s%s.%02d.saturated.fits" % (output_dir, datatype, obsid, ota)
    stdout_write("Writing output: %s\n" % (output_filename))

    clobberfile(output_filename)
    out_hdulist.writeto(output_filename, clobber=True)

    if (verbose):
        print "some of the saturated pixels:\n",final_cat[0:10,:]

    #numpy.savetxt("test", final_cat)
    #print full_coords.shape
        
    return final_cat
    


def map_persistency_effects(hdulist, verbose=False):

    mask_thisframe_list = {}
    mask_timeseries_list = {}

    if (verbose): stdout_write("Creating persistency masks ...")
    saturated_pixels_total = 0
    extensions_with_saturated_pixels = 0
    pixels_masked_out_thisframe = 0
    pixels_masked_out_timeseries = 0

    #
    # Check all cells in this file (for on OTA)
    #
    for ext in range(len(hdulist)):
        # Skip extensions that are no Image HDUs
        if (str(type(hdulist[ext])) != "<class 'pyfits.hdu.image.ImageHDU'>"):
            continue

        extname = hdulist[ext].header['EXTNAME']
        if (verbose): stdout_write("Working on extension %s (%d)\n" % (extname, ext))

        # Find all saturated pixels (values >= 65K)
        data = hdulist[ext].data
        saturated = (data >= 65535)

        # Skip this cell if no pixels are saturated
        number_saturated_pixels = numpy.sum(saturated)
        if (number_saturated_pixels <= 0):
            continue

        if (verbose): print "number of saturated pixels:", number_saturated_pixels

        saturated_pixels_total += number_saturated_pixels
        extensions_with_saturated_pixels += 1

        # Do some book-keeping preparing for the masking
        rows, cols = numpy.indices(data.shape)

        mask_thisframe = numpy.zeros(shape=data.shape)
        mask_thisframe = mask_thisframe > 1
        mask_time      = numpy.zeros(shape=data.shape)
        mask_time      = mask_time > 1
        #mask_time.fill(False)

        saturated_rows = rows[saturated]
        saturated_cols = cols[saturated]

        unique_cols = set(saturated_cols)

        #
        # Now convert the list of saturated pixels into a map
        #

        # Old, slow method
        if (False):
            for i in range(saturated_rows.shape[0]):
                mask_up = (cols == saturated_cols[i]) & (rows >= saturated_rows[i])
                mask_down = (cols == saturated_cols[i]) & (rows <= saturated_rows[i])
                #print "this:",mask_up.shape, mask_down.shape

                mask_thisframe = (mask_thisframe) | (mask_up)
                mask_time      = (mask_time)      | (mask_down)

        # New, optimized and way faster method
        row_ids, col_ids = numpy.indices((data.shape[0],1))
        for col in unique_cols:
            #print "working on col",col #saturated[col,:]
            this_col_saturated = row_ids[saturated[:,col]]
            #print "saturated in this col",this_col_saturated
            min_y = numpy.min(this_col_saturated)
            max_y = numpy.max(this_col_saturated)

            mask_thisframe[min_y:, col] = True
            mask_time[:max_y, col] = True

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

    # return two maps:
    # mask_thisframe:  Masks where all saturated pixels and 
    #                  columns above the saturated pixels are masked
    # mask_timeseries: Mask with all persistent pixels (saturated pixels 
    #                  and the rows below) are masked
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
            if (verbose): print "Found better match: %s (MJD=%.6f, or %d secs)" % (
                latest_map, file_mjd, min_delta_mjd*86400)

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


#
# This routine converts the masks into the persistency map that will
# a) be stored on disk to keep track of the persistency and b) be used
# to compute the correction for the science frame
#
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
        
#
# Mask out all saturated or saturation-effected pixels in the current
# frame (i.e. the one where pixels are saturated)
#
def apply_mask_to_data(mask, data):

    out = data
    out[mask] = numpy.NaN

    return out

#
# Compute the actual persistency correction from the persistency map.
#
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
        pers_dir = cmdline_arg_set_or_default("-persistency", "./")
        outputfile = "%s/persistency_map_00000000T000000.fits" % (pers_dir)
        
        # If this flag is set, simply create a new persistency map
        create_new_persistency_map(None, write_fits=outputfile)

        # Quit the program right here
        sys.exit(0)
        
    if (cmdline_arg_isset('-findmap')):
        directory = get_clean_cmdline()[1]
        mjd = float(get_clean_cmdline()[2])
        find_latest_persistency_map(directory, mjd, verbose=True)
        sys.exit(0)

    if (cmdline_arg_isset('-makecat')):
        output_dir = cmdline_arg_set_or_default('-persistency', '.')
        verbose = cmdline_arg_isset("-verbose")
        for filename in get_clean_cmdline()[1:]:
            create_saturation_catalog(filename, output_dir=output_dir, verbose=verbose)
        sys.exit(0)

    inputfile = sys.argv[1]
    hdulist = pyfits.open(inputfile)
    persistency_map_in = sys.argv[2]

    outputfile = sys.argv[3]
    persistency_map_out = sys.argv[4]

    hdulist = pyfits.open(inputfile)
    mjd = hdulist[0].header['MJD-OBS']


    mask_thisframe, mask_timeseries = map_persistency_effects(hdulist, verbose=True)

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
