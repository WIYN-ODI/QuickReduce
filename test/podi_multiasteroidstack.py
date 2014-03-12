#!/usr/bin/env python

import os, sys
d,_=os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append("%s/../"%d)

import podi_asteroids
import podi_swarpstack
from podi_commandline import *
import podi_logging
import logging
import astropy.io.votable
import math

if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    params = podi_swarpstack.read_swarp_params()
    inputlist = get_clean_cmdline()[2:]

    target_name = get_clean_cmdline()[1]

    print params
    print options

    # First, stack all frames without any non-sidereal correction
    outputfile = "%s__sidereal.fits" % (target_name)
    if (not os.path.isfile(outputfile)):
        podi_swarpstack.swarpstack(outputfile, inputlist, params, options)

    # Get all asteroids in the vicinity
    print "Looking for asteroids in frame",outputfile
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

        print "\n\n\n\nStarting work on object %s\n\n\n\n" % (names[obj])
        options['nonsidereal']['dra'] = pm_ra[obj]
        options['nonsidereal']['ddec'] = pm_dec[obj]
        params['use_nonsidereal'] = True

        print "Creating stack for",names[obj]
        print "Non-sidereal rate: Ra=%f ''/hr  Dec=%f ''/hr" % (
            options['nonsidereal']['dra'], options['nonsidereal']['ddec'])

        obj_name = names[obj]
        print "___%s___" % obj_name,"-->",
        obj_name = obj_name.replace(" ", "_")
        print obj_name
                
        outputfile = "%s__%s.fits" % (target_name, obj_name)
        print outputfile

        # Create the stack-file so we can re-run just this one again later 
        stack_conf = "%s__%s.stack" % (target_name, obj_name)
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
            podi_swarpstack.swarpstack(outputfile, inputlist, params, options)

        #break

    podi_logging.shutdown_logging(options)
