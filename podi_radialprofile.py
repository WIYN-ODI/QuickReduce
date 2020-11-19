#! /usr/bin/env python3
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
import ephem

from podi_definitions import *


def add_circle(buffer, center_x, center_y, radius, amplitude):

    x, y = numpy.indices(buffer.shape)

    dx = x - center_x
    dy = y - center_y
    d2 = dx*dx + dy*dy

    print(dx[0:10,0:10])
    print(dy[0:10,0:10])
    print(d2[0:10,0:10])

    #print d2[995:1005, 995:1005]
    
    tmp_buffer = numpy.zeros(shape=buffer.shape)
    tmp_buffer[d2 < radius*radius] = amplitude

    buffer += tmp_buffer

    return

def add_annulus(buffer, center_x, center_y, radius_i, radius_o, amplitude):

    x, y = numpy.indices(buffer.shape)

    dx = x - center_x
    dy = y - center_y
    d2 = dx*dx + dy*dy

    #print dx[0:10,0:10]
    #print dy[0:10,0:10]
    #print d2[0:10,0:10]

    #print d2[995:1005, 995:1005]
    
    tmp_buffer = numpy.zeros(shape=buffer.shape)

    selected_pixels = (d2 < radius_o*radius_o) & (d2 > radius_i*radius_i)
    tmp_buffer[selected_pixels] = amplitude

    buffer += tmp_buffer

    return

def create_from_template(command_file, buffer):
    
    # Load command file 
    cmdfile = open(command_file, "r")
    cmds = cmdfile.readlines()
    print(cmds)
    
    for i in range(len(cmds)):
        line = cmds[i]
        if (line[0] == "#"):
            continue
        
        items = line.strip().split()
        print(items)

        shape = items[0]
        if (shape == "fillcircle"):
            center_x = float(items[1])
            center_y = float(items[2])
            radius = float(items[3])
            amplitude = float(items[4])
            
            add_circle(buffer, center_x, center_y, radius, amplitude)

        if (shape == "annulus"):
            center_x = float(items[1])
            center_y = float(items[2])
            radius_i = float(items[3])
            radius_o = float(items[4])
            amplitude = float(items[5])
            
            add_annulus(buffer, center_x, center_y, radius_i, radius_o, amplitude)

            
    return buffer


def optimize_center(data, center_x, center_y):

    dx = center_x - x
    dy = center_y - y
    d = numpy.sqrt(dx*dx + dy*dy)

    radius_1d = d.ravel()
    data_1d = data.ravel()

    # Compute radial scatter
    # to do so, first sort both arrays by distance

    slice = 10 #pixels, approx 1''
    max_radius = math.sqrt(data.shape[0]*data.shape[0] + data.shape[1]*data.shape[1])
    max_radius = cmdline_arg_set_or_default('-maxrad', max_radius)
    n_slices = int(math.ceil(max_radius/slice))

    values = numpy.zeros(shape=(n_slices,4))

    #profile = open("mean_profile.out", "w")
    
    for n_out in range(n_slices):

        inner_radius = n_out * slice
        outer_radius = (n_out+1)    * slice
        
        in_this_bin = (d >= inner_radius) & (d < outer_radius) & (data > -1e9)
        
        values[n_out,0] = numpy.median(data[in_this_bin])
        values[n_out,1] = numpy.std(data[in_this_bin])

        values[n_out,2] = inner_radius
        values[n_out,3] = outer_radius

        #p = matplotlib.pyplot.figure()
        #p.plot(values[:,0], values[2])
        #p.show
        
        #print >> profile, inner_radius, outer_radius, values[n_out,0], values[n_out,1]
        
    return center_x, center_y, radius_1d, data_1d, values





#################################
#
# Important note:
# Many x/y values are swapped, because fits data is arranged in y/x coordinates, not x/y
#
#################################
if __name__ == "__main__":

    # Read in the input parameters
    fitsfile = sys.argv[1]

    extension = sys.argv[2]

    data_output_file = sys.argv[5]
    
    hdulist = pyfits.open(fitsfile)
    
    if (sys.argv[3] == "from" and sys.argv[4] == "file"):
        filter = hdulist[0].header['FILTER']
        if (filter in pupilghost_centers):
            if (extension in pupilghost_centers[filter]):
                center_y, center_x = pupilghost_centers[filter][extension]
            else:
                print("Couldn't find center for this extension")
                sys.exit(0)
        else:
            print("Could find this filter in list")
            sys.exit(0)
    else:
        center_y = float(sys.argv[3])
        center_x = float(sys.argv[4])

    print ("Using center position %d, %d" % (center_y, center_x))

    # Loop over all extensions
    # For now only use the first one, hard enough

    for ota_id in range(0, len(hdulist)):

        ota = hdulist[ota_id]

        extname = "---"
        if ('EXTNAME' in ota.header):
            extname = ota.header["EXTNAME"]
        
        if (extname == extension):
            # Now create the radial profile
            print("found extension %s\n" % (extname))
            
            ota.data = ota_data

            prebinned = rebin_image(ota.data[smaller_box:, smaller_box:], 4)
            x, y = numpy.indices(prebinned.shape)
            cx, cy, radius_1d, data_1d, values = optimize_center(prebinned, (center_x-smaller_box)/4., (center_y-smaller_box)/4.)

            #x, y = numpy.indices(ota.data.shape)
            #cx, cy, radius_1d, data_1d = optimize_center(ota.data, center_x, center_y)

            output = open(data_output_file, "w")
            for i in range(values.shape[0]):
                print (values[i,2], values[i,3], values[i,0], values[i,1], file=output)
            output.close()
            
            output = open(data_output_file+".samples", "w")
            for i in range(0, radius_1d.shape[0]):
                print (i, radius_1d[i], data_1d[i], file=output)
            output.close()
    
