#!/usr/bin/env python

#
# Copyright (C) 2014, Ralf Kotulla
#                     kotulla@uwm.edu
#
# All rights reserved
#

import os, sys
d,_=os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append("%s/../"%d)

import podi_asteroids
import podi_swarpstack
from podi_commandline import *
from podi_definitions import *
import podi_logging
import logging
import astropy.io.votable
import astropy.io.fits
import astropy.wcs
import math
import numpy
import ephem
import time
import scipy.stats
import Image
import ImageDraw
import astLib
import astLib.astWCS
import Image
import ImageDraw

try:
    import cPickle as pickle
except ImportError:
    import pickle


if __name__ == "__main__":

    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    logger = logging.getLogger("LightCurveMovie")

    #############################################################################
    #
    # Get all parameters we need to describe the motion of the object
    #
    #############################################################################

    source_coords = get_clean_cmdline()[1]
    # interpret source_coords
    sc_items = source_coords.split(",")
    src_ra = float(sc_items[0])
    src_dec = float(sc_items[1])

    src_dra, src_ddec, src_mjd = 0., 0., 0.

    if (len(sc_items) == 5):
        src_dra = float(sc_items[2])
        src_ddec = float(sc_items[3])
        if (os.path.isfile(sc_items[4])):
            hdu = astropy.io.fits.open(sc_items[4])
            src_mjd = hdu[0].header['MJD-OBS']
            # src_mjd = hdu[0].header['MJD-STRT']
        else:
            src_mjd = float(sc_items[4])

    # print src_ra, src_dec, src_dra, src_ddec, src_mjd
    logger.info("Reference-MJD: %.6f" % (src_mjd))

    # Internally, all coordinates are corrected for cos(dec)
    cos_declination = math.cos(math.radians(src_dec))
    #src_ra *= cos_declination

    inv_cosdec = numpy.array([1./cos_declination, 1.0])
    basename = get_clean_cmdline()[2]

    filelist = get_clean_cmdline()[3:]
    print filelist

    intensity_range = [-0.01, 0.25]


    dec_width = 25. / 3600. # 10 arcsec
    ra_width = dec_width / cos_declination
    
    #############################################################################
    #
    # Now go through each file and create a small cutout with the 
    # source and its near surrounding
    #
    #############################################################################
    extension_list = [pyfits.PrimaryHDU()]

    for filename in filelist:

        # logger.info("Creating cutout for file %s ..." % (filename))

        # Open file and derive MJD
        hdulist = pyfits.open(filename)
        mjd = hdulist[0].header['MJD-OBS']


        src_radec = numpy.array([src_ra, src_dec]) + \
                    (mjd - src_mjd)*24 * numpy.array([src_dra, src_ddec]) / 3600. * inv_cosdec
        
        corners = numpy.array(
            [[src_radec[0] - ra_width, src_radec[1] - dec_width],
             [src_radec[0] - ra_width, src_radec[1] + dec_width],
             [src_radec[0] + ra_width, src_radec[1] - dec_width],
             [src_radec[0] + ra_width, src_radec[1] + dec_width],
         ])
  
        logger.info("src-coords: %.6f %.6f" % (src_radec[0], src_radec[1]))


        #
        # Now find out in what extension we need to look for the source
        #
        source_extension = None
        for extension in hdulist:
            # print extension
            extname = extension.header['EXTNAME'] if 'EXTNAME' in extension.header else "???"
            if (not is_image_extension(extension)):
                # print extname,"is not a valid image extension"
                continue
            wcs = astLib.astWCS.WCS(extension.header, mode='pyfits')
            coord = numpy.array(wcs.wcs2pix(src_radec[0], src_radec[1]))
            img_size = numpy.array([extension.header['NAXIS2'], extension.header['NAXIS1']])

            print extname,"--->", coord
            if ( ((coord>img_size) | (coord<0)).any() ):
                # This is not the extension we are looking for
                continue

            # print src_radec,"found at position",coord,"in extension",extname
            source_extension = extension
            break

        if (source_extension == None):
            # This means we couldn't locate the extension with the source
            logger.warning("We couldn't locate the extension with the source")
            continue

        # print src_radec,"found at position",coord,"in extension",extname, source_extension.header['EXTNAME']

        #
        # Now determine the pixel-coordinates of the four corners 
        #
        corners_xy = numpy.array(wcs.wcs2pix(corners[:,0], corners[:,1]))
        # print corners_xy

        #########################################################################
        #
        # Cutout the area and insert into a small thumbnail image
        #
        #########################################################################
        data = source_extension.data.T

        corner_min = numpy.floor(numpy.min(corners_xy, axis=0)).astype(numpy.int)
        corner_max = numpy.ceil(numpy.max(corners_xy, axis=0)).astype(numpy.int)

        dimension = corner_max - corner_min
        cutout = numpy.zeros((dimension[0], dimension[1]))
        # print "image dimension:", dimension
        
        # print "stack dimension:", data.shape
            
        trunc_min = numpy.max([ corner_min, [0,0]], axis=0)
        trunc_max = numpy.min([ corner_max, data.shape], axis=0)
        # print "limited to valid area:",trunc_min, trunc_max
        
        # Now extract the data and insert it into the cutout
        insert_min = trunc_min - corner_min
        insert_max = insert_min + (trunc_max - trunc_min)
        
        # print "truncated area:", trunc_max-trunc_min

        # print "insert region:",insert_min, insert_max
        cutout[insert_min[0]:insert_max[0], insert_min[1]:insert_max[1]] = \
            data[trunc_min[0]:trunc_max[0], trunc_min[1]:trunc_max[1]]

        # Subtract the background level using the sky-levels from the header
        cutout -= source_extension.header['SKY_MEDI']

        # Set all NaN pixels to 0.
        cutout[numpy.isnan(cutout)] = 0.
        
        # And calibrate the frames to a zero-point of 25 to compensate for 
        # exposure times, atmospheric variations, etc...
        fluxscale = hdulist[0].header['FLXSCALE'] if 'FLXSCALE' in hdulist[0].header else 1.0
        cutout *= fluxscale

        # Now apply the intensity cuts and scale the image to the 
        # intensity range 0...255 in preparation for the export to jpg/png.
        img_normalized = (cutout - intensity_range[0]) / (intensity_range[1]-intensity_range[0])
        img_8bit = (img_normalized * 255.).T.astype(numpy.uint8)
        
        # Now create the image and write it to disk
        obsid = hdulist[0].header['OBSID']
        img_filename = "%s%s--%12.6f.png" % (basename, obsid, mjd)
        image = Image.fromarray(img_8bit)
        image.transpose(Image.FLIP_TOP_BOTTOM).save(img_filename, "PNG")
        logger.debug("writing png %s", img_filename)

        #
        # Also prepare a small thumbnail fits-file that contains the cutouts 
        # for all frames
        #
        thumb_hdu = pyfits.ImageHDU(data=cutout.T, header=source_extension.header)
        # Fix the WCS information
        extension_list.append(thumb_hdu)

        # break

    datacube = pyfits.HDUList(extension_list)
    datacube.writeto("%s.fits" % (basename), clobber=True)

    podi_logging.shutdown_logging(options)
