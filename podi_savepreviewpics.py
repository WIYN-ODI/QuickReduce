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

"""

This module contains a small tool to create preview JPGs to be used and 
displayed in the ODI observing GUI.

It's using multiple processes to speed up execution.

How to use
----------
``./podi_savepreviewpics.py /some/output/dir input1.fits input2.fits \
input....fits -z1=0 -z2=65000 -ncpus=4``

"""

import sys
import os
import pyfits
import numpy
# import scipy
import math

import Image
import ImageDraw

import multiprocessing
import Queue
import time

numpy.seterr(divide='ignore', invalid='ignore')


limit_overexposed = 55000
overexposed = [1.0, 0.0, 0.0]

crossout_missing_otas = True
from podi_definitions import *
from podi_commandline import *

def create_quickview(filename, output_directory, z1=None, z2=None, verbose=False, clobber=True):

    if (not os.path.isfile(filename)):
        return

    dirname, filebase = os.path.split(filename)
    pos = filebase.find(".fits")
    filename_without_fits = filebase[:pos]
    print filename_without_fits
    
    hdulist = pyfits.open(filename)

    binning = get_binning(hdulist[1].header)
    # print "Using binng factor",binning, get_binning(hdulist[1].header)
    # print hdulist[0].header
    
    sizex, sizey = get_collected_image_dimensions(binning)
    sizex /= 8
    sizey /= 8

    dataframe = numpy.zeros(shape=(sizey,sizex))
    dataframe[:,:] = numpy.NaN
    
    #
    # For all cells, perform an overscan subtraction and insert the binned
    # frame into the full OTA view
    #
    for extension in range(1, len(hdulist)):
        if (not type(hdulist[extension]) == pyfits.hdu.image.ImageHDU):
            continue

        cell = hdulist[extension]
        # print cell.header['
        cellx = cell.header['WN_CELLX']
        celly = cell.header['WN_CELLY']
        # print "working on extension",extension, cellx, celly,

        # Extract the science section
        datasec = extract_datasec_from_cell(cell.data, binning)

        # Determine the overscan level
        biassec = extract_biassec_from_cell(cell.data, binning)
        overscan = numpy.mean(biassec[numpy.isfinite(biassec)])

        # print "overscan level",overscan
        datasec -= overscan

        # Bin the cell image x8 
        binned = rebin_image(datasec, 8)

        x1, x2, y1, y2 = cell2ota__get_target_region(cellx, celly, binning=binning)
        x1 = int(math.floor(x1/8))
        x2 = x1 + binned.shape[1]
        y1 = int(math.floor(y1/8))
        y2 = y1 + binned.shape[0]

        # print x1, x2, y1, y2
        dataframe[y1:y2, x1:x2] = binned
        
    #
    # Now we are through all OTA/extensions, compute the median value and stds 
    # so we can scale the frames accordingly
    #
    if (verbose):
        stdout_write("   Finding best intensity levels ...")

    if (z1 == None or z2 == None):
        median = 0
        std = 1e8
        for looper in range(3):
             valid = (dataframe > (median - std)) & (dataframe < (median + 3*std))
             median = numpy.median(dataframe[valid])
             std = numpy.std(dataframe[valid])

    # Compute the intensity levels, making sure not to exceed the maximum possible range
    min_level = float(z1) if z1 != None else numpy.max([median-3*std,0])
    max_level = float(z2) if z2 != None else numpy.min([median+10*std,60000])
    # stdout_write(" using %d ... %d\n" % (min_level, max_level))

    #
    # Now that we have all the scaling factors, go ahead and create the preview images
    #

    if (verbose):
        stdout_write("   Creating jpeg for OTA")

    #
    # Create and save the larger binning x8 image
    #
    greyscale = (dataframe - min_level) / (max_level - min_level)
    greyscale[greyscale < 0] = 0.
    greyscale[greyscale > 1] = 1.
    image_filename = "%s/%s.bin8.jpeg" % (output_directory, filename_without_fits)
    image = Image.fromarray(numpy.uint8(greyscale*255))
    image.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")

    #
    # Create the binned x64 thumbnail
    #
    bin64 = rebin_image(dataframe, 8)
    greyscale = (bin64 - min_level) / (max_level - min_level)
    greyscale[greyscale > 1] = 1.
    greyscale[greyscale < 0] = 0.

    image_filename = "%s/%s.bin64.jpeg" % (output_directory, filename_without_fits)
    image = Image.fromarray(numpy.uint8(greyscale*255))
    image.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")

    return


def multi_convert(queue, xxx):

    while (True):
        filename, output_directory, z1, z2 = queue.get()
        if (filename == None):
            queue.task_done()
            return

        create_quickview(filename, output_directory, z1, z2, verbose=False, clobber=clobber)
    return

if __name__ == "__main__":

#    filename = sys.argv[1]
#    print filename

    output_directory = "."
    output_directory = sys.argv[1]

    clobber = not cmdline_arg_isset("-noclobber")
    if (not clobber):
        stdout_write("Activating no-clobber mode!\n")

    z1 = cmdline_arg_set_or_default("-z1", None)
    z2 = cmdline_arg_set_or_default("-z2", None)
    
    queue = multiprocessing.JoinableQueue()
    processes = []

    number_cpus = int(cmdline_arg_set_or_default("-ncpus", 6))
    
    for filename in sys.argv[2:]:
        queue.put( (filename, output_directory, z1, z2) )
    
    # Create all processes to handle the actual reduction and combination
    for i in range(number_cpus):
        p = multiprocessing.Process(target=multi_convert, args=(queue, None))
        p.start()
        processes.append(p)
        time.sleep(0.01)
        queue.put( (None, None, None, None) )
        
    sys.exit(0)
    
