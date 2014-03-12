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
import math
import numpy
import ephem

if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    logger = logging.getLogger("AstroidField")

    params = podi_swarpstack.read_swarp_params()
    inputlist = get_clean_cmdline()[2:]

    target_name = get_clean_cmdline()[1]

    smaller_region = 5 # arcmin on a side

    print params
    print options

    # First, stack all frames without any non-sidereal correction
    outputfile = "%s__sidereal.fits" % (target_name)
    if (not os.path.isfile(outputfile)):
        logger.info("Creating the sidereal stack")
        podi_swarpstack.swarpstack(outputfile, inputlist, params, options)

    # Get all asteroids in the vicinity
    logger.info("Looking for asteroids in frame %s" % (outputfile))
    votab = podi_asteroids.get_asteroid_list_from_fitsfile(outputfile)


    obj_ra = votab.array['RA']
    obj_dec = votab.array['DEC']
    names = votab.array['Name']
    pm_ra = votab.array['dracosdec']
    pm_dec = votab.array['ddec']

    middle_ref = int(math.floor(len(inputlist)/2))

    # Use the first frame in the input list as non-sidereal reference frame
    hdu1 = pyfits.open(inputlist[middle_ref])
    
    options['nonsidereal'] = {}
    options['nonsidereal']['ref'] = inputlist[middle_ref]
    options['nonsidereal']['ref_mjd'] = hdu1[0].header['MJD-OBS']

    for obj in range(names.shape[0]):

        inputlist = get_clean_cmdline()[2:]

        print "\n\n\n"
        logger.info("Starting work on object %s" % (names[obj]))
        print "\n\n\n"

        options['nonsidereal']['dra'] = pm_ra[obj]
        options['nonsidereal']['ddec'] = pm_dec[obj]
        params['use_nonsidereal'] = True

        logger.debug("Creating stack for %s" % (names[obj]))
        logger.info("Non-sidereal rate: Ra=%f ''/hr  Dec=%f ''/hr" % (
            options['nonsidereal']['dra'], options['nonsidereal']['ddec']))

        obj_name = names[obj]
        #print "___%s___" % obj_name,"-->",
        obj_name = obj_name.replace(" ", "_")
        #print obj_name
                
        outputfile = "%s__%s.fits" % (target_name, obj_name)
        outputfile_weight = "%s__%s.weight.fits" % (target_name, obj_name)
        logger.debug("Setting outputfie: %s" % (outputfile))

        # Create the stack-file so we can re-run just this one again later 
        stack_conf = "%s__%s.stack" % (target_name, obj_name)
        logger.debug("Creating the stack-config file: %s" % (stack_conf))
        sf = open(stack_conf, "w")
        print >>sf, outputfile
        for kw in sys.argv[1:]:
            if (kw.startswith("-")):
                print >>sf, kw
        print >>sf, "-nonsidereal=%f,%f,%s" % (pm_ra[obj], pm_dec[obj], inputlist[middle_ref])
        for fn in inputlist:
            print >>sf, fn
        sf.close()
        
        if (not os.path.isfile(outputfile)):
            logger.debug("Calling swarpstack...")
            podi_swarpstack.swarpstack(outputfile, inputlist, params, options)
        else:
            logger.info("%s already exists, skipping" % (outputfile))

        if (not os.path.isfile(outputfile)):
            # Something went wrong
            logger.error("Problem: the stacked output file doesn't exist")
            continue

        logger.info("Creating cut-out centered on object ...")
        cutout_file = "%s__%s.cutout.fits" % (target_name, obj_name)
        cutout_weight = "%s__%s.cutout.weight.fits" % (target_name, obj_name)
        if (smaller_region > 0):
            # Create a thumbnail centered on the asteroids position
            hdu = astropy.io.fits.open(outputfile)
            wcs = astropy.wcs.WCS(hdu[0].header)
            data = hdu[0].data.T

            coord_j2000 = ephem.Equatorial(obj_ra[obj], obj_dec[obj], epoch=ephem.J2000)
            ra = numpy.degrees(coord_j2000.ra)
            dec = numpy.degrees(coord_j2000.dec)

            dec_width = smaller_region / 60. / 2.
            ra_width = dec_width / math.cos(math.radians(dec))
            print ra, dec, ra_width, dec_width

            corners = [[ra - ra_width, dec - dec_width],
                       [ra - ra_width, dec + dec_width],
                       [ra + ra_width, dec - dec_width],
                       [ra + ra_width, dec + dec_width],
            ]
            corners_xy = wcs.wcs_world2pix(corners, 0)
            corner_min = numpy.floor(numpy.min(corners_xy, axis=0)).astype(numpy.int)
            corner_max = numpy.ceil(numpy.max(corners_xy, axis=0)).astype(numpy.int)

            print wcs
            print corners
            print corners_xy
            print "corner min/max", corner_min, corner_max
            
            dimension = corner_max - corner_min
            cutout = numpy.zeros((dimension[0], dimension[1]))
            cutout[:,:] = numpy.NaN
            print "image dimension:", dimension

            print "stack dimension:", data.shape

            trunc_min = numpy.max([ corner_min, [0,0]], axis=0)
            trunc_max = numpy.min([ corner_max, data.shape], axis=0)
            print "limited to valid area:",trunc_min, trunc_max
            
            # Now extract the data and insert it into the cutout
            insert_min = trunc_min - corner_min
            insert_max = insert_min + (trunc_max - trunc_min)

            print "truncated area:", trunc_max-trunc_min

            print "insert region:",insert_min, insert_max
            cutout[insert_min[0]:insert_max[0], insert_min[1]:insert_max[1]] = \
                data[trunc_min[0]:trunc_max[0], trunc_min[1]:trunc_max[1]]

            hdu[0].data = cutout.T
            hdu[0].header['CRPIX1'] -= trunc_min[0]
            hdu[0].header['CRPIX2'] -= trunc_min[1]
            
            clobberfile(cutout_file)
            hdu.writeto(cutout_file, clobber=True)
            hdu.close()

            # Repeat the same procedure with the weight map
            hdu_w = astropy.io.fits.open(outputfile_weight)
            data_w = hdu_w[0].data.T
            cutout_w = numpy.zeros((dimension[0], dimension[1]))
            cutout_w[insert_min[0]:insert_max[0], insert_min[1]:insert_max[1]] = \
                data_w[trunc_min[0]:trunc_max[0], trunc_min[1]:trunc_max[1]]
            hdu_w[0].data = cutout_w.T
            hdu_w[0].header['CRPIX1'] -= trunc_min[0]
            hdu_w[0].header['CRPIX2'] -= trunc_min[1]

            clobberfile(cutout_weight)
            hdu_w.writeto(cutout_weight, clobber=True)
            hdu_w.close()

        # Next run source extractor on the frame
        compute_photometry = True
        if (compute_photometry):
            pass

        #break

    podi_logging.shutdown_logging(options)
