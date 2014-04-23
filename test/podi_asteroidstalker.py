#!/usr/bin/env python

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

try:
    import cPickle as pickle
except ImportError:
    import pickle


if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    logger = logging.getLogger("AsteroidStalker")

    target_name = get_clean_cmdline()[1]
    ref_file = get_clean_cmdline()[2]
    inputlist = get_clean_cmdline()[3:]


    hdulist = astropy.io.fits.open(ref_file)
    rel_ra = rel_dec = rel_mjd = 0
    rel_ra = hdulist[0].header['_NSID_RA'] if ('_NSID_RA' in hdulist[0].header) else 0.0
    rel_dec = hdulist[0].header['_NSID_DE'] if ('_NSID_DE' in hdulist[0].header) else 0.0
    rel_mjd = hdulist[0].header['_NSID_RT'] if ('_NSID_RT' in hdulist[0].header) else \
              hdulist[0].header['MJD-OBS'] if ('MJD-OBS' in hdulist[0].header) else 0.0

    mjd_start = hdulist[0].header['MJD-STRT']
    mjd_end   = hdulist[0].header['MJD-END']

    wcs = astropy.wcs.WCS(hdulist[0].header)

    margin = 1 # arcmin

    candidates = None
    with open("asteroidmetry.pickle", "rb") as pf:
        candidates = pickle.load(pf)

    all_bglevels = {}
    all_bglevel_std = {}

    min_level = -0.2
    max_level = 2.5

    add_reference = cmdline_arg_isset("-addref")

    for cand in candidates:
        print cand['rate']

        if (not 'name' in cand):
            continue

        name = cand['name']
        name = name.replace(" ", "_").replace(":", "")

        min_ra, max_ra = numpy.min(cand['sex'][:,0]), numpy.max(cand['sex'][:,0])
        min_dec, max_dec = numpy.min(cand['sex'][:,1]), numpy.max(cand['sex'][:,1])
        cos_dec = numpy.cos(numpy.radians(numpy.max(numpy.fabs(cand['sex'][:,1]))))

        ra_width = margin / 60. / cos_dec
        dec_width = margin / 60.

        corners = [[min_ra - ra_width, min_dec - dec_width],
                   [min_ra - ra_width, max_dec + dec_width],
                   [max_ra + ra_width, min_dec - dec_width],
                   [max_ra + ra_width, max_dec + dec_width],
               ]
        print corners

        corners_xy = wcs.wcs_world2pix(corners, 0)
        corner_min = numpy.floor(numpy.min(corners_xy, axis=0)).astype(numpy.int)
        corner_max = numpy.ceil(numpy.max(corners_xy, axis=0)).astype(numpy.int)

        print wcs
        print corners
        print corners_xy
        print "corner min/max", corner_min, corner_max

        print "\n"*10,name,"\n"*10

        # Now open each of the input files, and extract the sub-image
        for infile in inputlist:
            inhdu = astropy.io.fits.open(infile)
            data = inhdu[0].data.T

            # compute the background level
            logger.info("Computing bg-levels")
            #all_bglevels[infile] = 4.0
            #all_bglevel_std[infile] = 1.5
            if (not infile in all_bglevels):
                valid = data[numpy.isfinite(data)]
                filtered = valid #, mask = three_sigma_clip(valid, return_mask=True)
                std = numpy.std(filtered)
                median = numpy.median(filtered)
                print infile, median, std
                all_bglevels[infile] = median
                all_bglevel_std[infile] = std

            logger.info("median = %f +/- %f" % (all_bglevels[infile], all_bglevel_std[infile]))

            data -= all_bglevels[infile]
            #print "data min/max=",numpy.min(data[numpy.isfinite(data)]), numpy.max(data[numpy.isfinite(data)])

            # min_level = all_bglevels[infile] - 2 * all_bglevel_std[infile]
            # max_level = all_bglevels[infile] + 10 * all_bglevel_std[infile]
            #data[data < min_level] = min_level
            #data[data > max_level] = max_level
            #data -= min_level
            #print "data min/max #2=",numpy.min(data[numpy.isfinite(data)]), numpy.max(data[numpy.isfinite(data)])
            # Now the data is in the range 0 ... max_level-min_level

            dimension = corner_max - corner_min
            cutout = numpy.zeros((dimension[0], dimension[1]))
            # cutout[:,:] = numpy.NaN
            print "image dimension:", dimension
        
            print "stack dimension:", hdulist[0].data.shape
            
            trunc_min = numpy.max([ corner_min, [0,0]], axis=0)
            trunc_max = numpy.min([ corner_max, hdulist[0].data.shape], axis=0)
            print "limited to valid area:",trunc_min, trunc_max
        
            # Now extract the data and insert it into the cutout
            insert_min = trunc_min - corner_min
            insert_max = insert_min + (trunc_max - trunc_min)
        
            print "truncated area:", trunc_max-trunc_min

            print "insert region:",insert_min, insert_max
            cutout[insert_min[0]:insert_max[0], insert_min[1]:insert_max[1]] = \
                data[trunc_min[0]:trunc_max[0], trunc_min[1]:trunc_max[1]]

            cutout[numpy.isnan(cutout)] = 0.
            if (add_reference):
                print "adding reference image"
                cutout[insert_min[0]:insert_max[0], insert_min[1]:insert_max[1]] += \
                    hdulist[0].data[trunc_min[0]:trunc_max[0], trunc_min[1]:trunc_max[1]]
            cutout[cutout < min_level] = min_level
            cutout[cutout > max_level] = max_level
            cutout -= min_level

            print "normalizing:",numpy.arcsinh(max_level-min_level)
            cutout_asinh = numpy.arcsinh(cutout) / (numpy.arcsinh(max_level-min_level)) * 255.
            print numpy.min(cutout_asinh), numpy.max(cutout_asinh)

            image_filename = "%s%s_%s.jpg" % (target_name, name, infile)
            image = Image.fromarray(cutout_asinh.T.astype(numpy.uint8))

            #
            # Find the position of the source in this image
            #
            # get timestamp for frame
            img_mjd = inhdu[0].header['MJD-OBS']+0.5*inhdu[0].header['EXPTIME']/86400.
            # find which is the closest source data
            d_mjd = cand['mjd'] - img_mjd
            # closest
            closest_mjd = numpy.argmin(numpy.fabs(d_mjd))
            # if within 10 seconds of this frame
            if (numpy.fabs(d_mjd[closest_mjd]) < 150./86400.):
                # This source is detected
                src_detected = True
                src_radec = cand['sex'][closest_mjd,0:2]
            else:
                delta_mjd = d_mjd[closest_mjd]
                cos_dec_src = math.cos(math.radians(cand['sex'][closest_mjd,1]))
                d_dec = delta_mjd * 24. / 3600. * cand['rate'][1]
                d_ra = delta_mjd * 24. / 3600. * cand['rate'][0] / cos_dec_src
                src_radec = cand['sex'][closest_mjd,0:2] - [d_ra, d_dec]
                print "undetected src_radec=",src_radec
                src_detected = False

            # Figure out what the pixel coordinates of the source are
            src_xy = wcs.wcs_world2pix([src_radec], 0)[0]
            # convert this into coordinates in the cut-out image
            print "undetected src_xy=",src_xy
            cutout_src_xy = src_xy - corner_min

            print cutout_src_xy
            draw = ImageDraw.Draw(image) 
            # draw.line((100,200, 150,300), fill=128)
            if (src_detected):
                # Now draw two vertical lines pointing at the position of the source
                draw.line((cutout_src_xy[0]-20,cutout_src_xy[1], cutout_src_xy[0]-10, cutout_src_xy[1]), fill=255)
                draw.line((cutout_src_xy[0],cutout_src_xy[1]-20, cutout_src_xy[0], cutout_src_xy[1]-10), fill=255)
            else:
                bw=5
                draw.rectangle([cutout_src_xy[0]-bw, cutout_src_xy[1]-bw, 
                                cutout_src_xy[0]+bw, cutout_src_xy[1]+bw], outline=255)

            image.transpose(Image.FLIP_TOP_BOTTOM).save(image_filename, "JPEG")
            logger.info("writing jpeg %s", image_filename)

            inhdu.close()


    podi_logging.shutdown_logging(options)
