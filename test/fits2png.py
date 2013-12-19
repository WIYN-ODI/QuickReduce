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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

def rebin_image(data, binfac):

    if (binfac < 1):
        stdout_write("Rebinning at the moment only supports binning to larger pixels with binfac>1\n")
        return None
    elif (binfac == 1):
        return data

    out_size_x, out_size_y = int(math.ceil(data.shape[0]*1.0/binfac)), int(math.ceil(data.shape[1]*1.0/binfac))

    if (out_size_x*binfac != data.shape[0] or out_size_y*binfac != data.shape[1]):
        # The input array size is not a multiple of the new binning
        # Create a slightly larger array to hold the data to be rebinned
        container = numpy.zeros(shape=(out_size_x*binfac, out_size_y*binfac))

        # And insert the original data
        container[0:data.shape[0], 0:data.shape[1]] = data[:,:]
    else:
        container = data 
        
    rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).mean(axis=-1).mean(axis=1)

    return rebinned

def create_quickview(filename, verbose=False, clobber=True):

    hdulist = pyfits.open(filename)

    if (verbose):
        print "Working on file %s ..." % (filename),

    if (verbose):
        print "   Finding contrast: ",
    for extension in range(0, len(hdulist)):

        if (type(hdulist[extension]) not in (pyfits.hdu.image.ImageHDU, pyfits.hdu.image.PrimaryHDU)):
            continue

        try:
            extname = hdulist[extension].header['EXTNAME']
        except:
            extname = "%d" % (extension)
        print extname,

        median = 0
        std = 1e8
        imgdata = rebin_image(hdulist[extension].data, 4)
        for looper in range(3):
            valid = (imgdata > (median - std)) & (imgdata < (median + 3*std))

            median = numpy.median(imgdata[valid])
            std = numpy.std(imgdata[valid])
            #print median, std, median-std, median+3*std

        # Compute the intensity levels, making sure not to exceed the maximum possible range
        min_level = numpy.max([median-1*std,0])
        max_level = numpy.min([median+8*std,60000])
        #stdout_write(" using %d ... %d\n" % (min_level, max_level))

        imgdata[numpy.isnan(imgdata)] = 0
        
        greyscale = (imgdata - min_level) / (max_level - min_level)
        greyscale[greyscale<0] = 0
        greyscale[greyscale>=1] = 1

        greyscale = numpy.sqrt(greyscale)
        #greyscale = numpy.log10(greyscale+1)/numpy.log10(2.0)

        #
        # Mark all overexposed pixels in a different color
        #
        channel_r = greyscale + 0
        channel_g = greyscale + 0
        channel_b = greyscale + 0

        channel_r[imgdata > limit_overexposed] = overexposed[0]
        channel_g[imgdata > limit_overexposed] = overexposed[1]
        channel_b[imgdata > limit_overexposed] = overexposed[2]

        im_r = Image.fromarray(numpy.uint8(channel_r*255))
        im_g = Image.fromarray(numpy.uint8(channel_g*255))
        im_b = Image.fromarray(numpy.uint8(channel_b*255))
        im_rgb = Image.merge('RGB', (im_r, im_g, im_b))

        image_filename = "%s.%s.jpg" % (filename[:-5], extname)
        im_rgb.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")

        # Delete all temporary images to keep memory demands low
        del im_r
        del im_g
        del im_b
        del im_rgb

    print " - done!"



if __name__ == "__main__":

#    filename = sys.argv[1]
#    print filename

    for filename in sys.argv[1:]:
        create_quickview(filename, verbose=True)
