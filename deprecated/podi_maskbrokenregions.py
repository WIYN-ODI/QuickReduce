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
import astropy.io.fits as pyfits
import numpy
import scipy

from podi_definitions import *

def mask_broken_regions(datablock, regionfile):

    counter = 0
    file = open(regionfile)
    for line in file:
        if (line[0:3] == "box"):
            coords = line[4:-2]
            coord_list = coords.split(",")
                        
            if (datablock is not None):
                x, y = int(float(coord_list[0])), int(float(coord_list[1]))
                dx, dy = int(0.5*float(coord_list[2])), int(0.5*float(coord_list[3]))
                #mask[y-dy:y+dy,x-dx:x+dx] = 1

                x1 = numpy.max([0, x-dx])
                x2 = numpy.min([datablock.shape[1], x+dx])
                y1 = numpy.max([0, y-dy])
                y2 = numpy.min([datablock.shape[0], y+dy])
                datablock[y1:y2, x1:x2] = numpy.nan

                # print x,x+dx,y,y+dy
            counter += 1

    file.close()
    print("Marked",counter,"bad pixel regions")
    return datablock


if __name__ == "__main__":

    cmd_name = sys.argv[0]
    print(cmd_name)
    print("files are in ",os.path.dirname(cmd_name))

    starting_dir = os.path.dirname(cmd_name)
    for filename in sys.argv[1:]:
        
        stdout_write("Next file %s ..." % filename)

        if (filename[-3:] == ".fz"):
            stdout_write(" this is a fz-packed file! ")

            funpacked = filename[:-3]
            if (os.path.isfile(funpacked)):
                stdout_write(" found already fz-unpacked file, using this one instead")
            else:
                os.system("funpack %s" % filename)
                filename = funpacked

        # Open filename and determine which OTA chip it is
        stdout_write(" opening file ")
        hdulist = pyfits.open(filename)
        fppos = hdulist[0].header['FPPOS'].strip()

        # Determine which region file we need
        region_file = "%s/bpm_%s.reg" % (starting_dir, fppos)

        # Make sure the region file exists, otherwise warn the user
        # and continue with the next frame
        if (os.path.isfile(region_file)):

            # Read data block
            image_data = hdulist[0].data

            # Apply the bad pixel regions to file, marking
            # all bad pixels as NaNs
            mask_broken_regions(image_data, region_file)

            # Delete two headers that don't seem to agree with the FITS format
            del hdulist[0].header['RSPRA']
            del hdulist[0].header['RSPDEC']

            # Form the output filename and write the masked frame
            output_filename = filename[0:-5]+".masked.fits"
            print(output_filename)
            hdulist[0].data = image_data
            hdulist.writeto(output_filename, overwrite=True)

        # Close file and free the memory
        hdulist.close()
        del hdulist
        stdout_write("done!")

    sys.exit(0)
    
    

