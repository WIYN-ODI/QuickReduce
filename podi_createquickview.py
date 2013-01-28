#! /usr/bin/env python

#
# (c) Ralf Kotulla for WIYN/pODI
#

import sys
import os
import pyfits
import numpy
import scipy
import math

import Image
import ImageDraw

"""
How to use: 
./podi_createquickview file.fits output_dir/

What it does:
The script opens the input file and first bins each OTA extension by 8x8. 
For some OTAs, depending on filter (also see the makeflatfield routine) the 
image data is used to derive some best-fit intensity level with which best to 
display the frame. Once these levels are known, each OTA is exported as jpeg-file,
and all OTAs are also exported as one combined jpeg. All jpegs are written into 
the specified output folder and labeled by their OTA ID and the unique observation
ID (essentially the date of observation).

If specified by the 'crossout_missing_otas' variable, all non-existing OTAs can
be crossed out in the full focalplane jpeg.

"""

limit_overexposed = 55000
overexposed = [1.0, 0.0, 0.0]

crossout_missing_otas = True
from podi_definitions import *

if __name__ == "__main__":

    filename = sys.argv[1]
    print filename

    output_directory = "."
    if (len(sys.argv) > 2):
        output_directory = sys.argv[2]

    hdulist = pyfits.open(filename)
    filter = hdulist[0].header['FILTER']
    obsid  = hdulist[0].header['OBSID'] 

    print "This is filter",filter

    list_of_otas_to_normalize = which_otas_to_use[filter]

    # Allocate some memory to hold the data we need to determine the
    # most suitable intensity levels
    binned_data = numpy.zeros(shape=(13*512*512), dtype=numpy.float32)
    binned_data[:] = numpy.NaN

    datapos = 0
    for extension in range(1, len(hdulist)):
        fppos = int(hdulist[extension].header['FPPOS'][2:4])

        try:
            index = list_of_otas_to_normalize.index(fppos)
        except:
            # We didn't find this OTA in the list, so skip it
            hdulist[extension].header.update("FF_NORM", 0, "Used in normalization")
            extension += 1
            continue

        stdout_write("\rReading OTA %02d" % (fppos))

        # Rebin the frame 8x8 to make it more manageable
        binned = numpy.reshape(hdulist[extension].data, (512,8,512,8)).mean(axis=-1).mean(axis=1)
        one_d = binned.flatten()
        binned_data[datapos:datapos+one_d.shape[0]] = one_d
        datapos += one_d.shape[0]
        del one_d
        del binned

    stdout_write(" done!\n")

    #
    # Now we are through all OTA/extensions, compute the median value and stds 
    # so we can scale the frames accordingly
    #
    stdout_write("Finding best intensity levels ...")
    median = 0
    std = 1e8
    binned_data = binned_data[0:datapos]
    for looper in range(3):
        valid = (binned_data > (median - std)) & (binned_data < (median + 3*std))
        #print numpy.sum(valid)
        median = numpy.median(binned_data[valid])
        std = numpy.std(binned_data[valid])
        #print median, std, median-std, median+3*std

    # Compute the intensity levels, making sure not to exceed the maximum possible range
    min_level = numpy.max([median-1*std,0])
    max_level = numpy.min([median+8*std,60000])
    stdout_write(" using %d ... %d\n" % (min_level, max_level))

    #
    # Now that we have all the scaling factors, go ahead and create the preview images
    #

    # Create space to hold the full 8x8 OTA focal plane
    full_focalplane = numpy.zeros(shape=(4096,4096))

    for extension in range(1, len(hdulist)):
        fppos = int(hdulist[extension].header['FPPOS'][2:4])
        stdout_write("\rCreating jpegs (%02d) ..." % fppos)
        fp_x = fppos % 10
        fp_y = math.floor(fppos / 10)
        
        hdulist[extension].data[numpy.isnan(hdulist[extension].data)] = 0
        binned = numpy.reshape(hdulist[extension].data, (512,8,512,8)).mean(axis=-1).mean(axis=1)
        
        greyscale = (binned - min_level) / (max_level - min_level)
        greyscale[greyscale<0] = 0
        greyscale[greyscale>=1] = 1

        greyscale = numpy.sqrt(greyscale)
        #greyscale = numpy.log10(greyscale+1)/numpy.log10(2.0)

        ffp_x = fp_x * 512
        ffp_y = fp_y * 512
        full_focalplane[ffp_x:ffp_x+512, ffp_y:ffp_y+512] = greyscale[:,:]

        image = Image.fromarray(numpy.uint8(greyscale*255))
        image_filename = "%s/%s.%02d.jpg" % (output_directory, obsid, fppos)
        image.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")
        
        #
        # Mark all overexposed pixels in a different color
        #
        channel_r = greyscale + 0
        channel_g = greyscale + 0
        channel_b = greyscale + 0

        channel_r[binned > limit_overexposed] = overexposed[0]
        channel_g[binned > limit_overexposed] = overexposed[1]
        channel_b[binned > limit_overexposed] = overexposed[2]

        im_r = Image.fromarray(numpy.uint8(channel_r*255))
        im_g = Image.fromarray(numpy.uint8(channel_g*255))
        im_b = Image.fromarray(numpy.uint8(channel_b*255))
        im_rgb = Image.merge('RGB', (im_r, im_g, im_b))
        image_filename = "%s/%s.%02d.rgb.jpg" % (output_directory, obsid, fppos)
        im_rgb.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")

        # Delete all temporary images to keep memory demands low
        del image
        del im_r
        del im_g
        del im_b
        del im_rgb

    #
    # Prepare the preview for the full focal plane
    #
    stdout_write("\rCreating jpegs (full-frame) ...")
    image = Image.fromarray(numpy.uint8(full_focalplane*255))

    # Add lines to indicate detector borders. Make sure to make them wider than 
    # just one pixel, otherwise the lines completely disappear if displaying the 
    # image zoomed out
    draw = ImageDraw.Draw(image)
    for i in range(1,8):
        for linewidth in range(-5,0):
            draw.line((0,i*512+linewidth,image.size[0],i*512+linewidth), fill=128)
            draw.line((i*512+linewidth,0,i*512+linewidth,image.size[1]), fill=128)

    if (crossout_missing_otas):
        # Now loop over all OTAs and mark the ones that do not exist
        for y in range(8):
            for x in range(8):
                tuple = (x,y)
                try:
                    index = available_ota_coords.index(tuple)
                except:
                    # We only get here if the OTA is not listed as available
                    # cross it out in the full focal plane image
                    #print "Strokign out ",tuple
                    for lw in range(-5,0):
                        draw.line((x*512+lw,y*512,(x+1)*512+lw,(y+1)*512), fill=128)
                        draw.line((x*512+lw,(y+1)*512,(x+1)*512+lw,y*512), fill=128)

    image_filename = "%s/%s.jpg" % (output_directory, obsid)
    image.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")
    del image

    stdout_write(" done!\n")
