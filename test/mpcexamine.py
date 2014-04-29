#!/usr/bin/env python

import numpy
import os
import sys
import time
import scipy, scipy.spatial

sys.path.append("../")

import astropy.io.fits

# Import ds9 and start XPA name-server
import ds9 as pyds9
pyds9.ds9_xpans()

import podi_asteroidmetry
from podi_definitions import *

if __name__ == "__main__":

    filename = sys.argv[1]
    # pyds9.ds9.openlist()


    ds9 = pyds9.ds9(target="mpcexamine", 
                    start="-scale zscale -scale scope global -file %s" % (filename), 
                    wait=25)

    hdulist = None
    target_name = "??????"

    # which aperture do we need to use
    magname = '4.0'

    while (True):

        if (hdulist == None):
            print "Loading file and associated catalog"
            ds9.set("file %s" % (filename))

            #
            # Open the fits file to get the timing information from the header
            #
            hdulist = astropy.io.fits.open(filename)

            #
            # Now load the catalog and find the source
            #
            catfile = filename[:-5]+".cat"
            catdata = numpy.loadtxt(catfile)
            # Correct for magzero
            magzero = hdulist[0].header['MAGZERO'] if 'MAGZERO' in hdulist[0].header else 25.
            catdata[:, SXcolumn['mag_aper_2.0']:SXcolumn['mag_aper_12.0']+1] += magzero
            catdata[:, SXcolumn['mag_auto']] += magzero

            cat_tree = scipy.spatial.cKDTree(catdata[:, 0:2])

        print "Load and prepare the file for inspection"
        answer = raw_input("Press enter when ready or 'q' for abort: ")
        print "answer was __%s__" % (answer)

        if (answer == 'q'):
            break

        if (len(answer) > 5 and answer.split()[0] == "file"):
            filename = answer.split()[1]
            hdulist = None
            continue

        if (len(answer) > 5 and answer.split()[0] == "target"):
            target_name = answer.split()[1]
            continue

        
        # hdu = ds9.get_pyfits()
        # print hdu
        # hdu.info()

        # ret = ds9.get(answer)
        # ret = ds9.set(answer)
        # print "\nds9 returned:",ret

        # Select object for which we are to look up the coordinates etc.
        ret = ds9.get('imexam coordinate wcs fk5 degrees')
        
        ra, dec = float(ret.split()[0]), float(ret.split()[1])
        print ra, dec

        # get frame MJD
        mjd = hdulist[0].header['MJD-OBS'] + hdulist[0].header['EXPTIME']/2/86400.

        # Search for the source within 2 arcsec
        coords = [ra, dec]
        counterparts = cat_tree.query_ball_point(coords, r=2./3600., p=2)

        matches = catdata[counterparts]
        print matches[:, 0:2]

        # Draw a blue cross-hair at the location of the source
        # ds9.set("regions fk5; circle %f %f 0.0001" % (
        #     matches[0, SXcolumn['ra']], matches[0, SXcolumn['dec']]))

        ret = podi_asteroidmetry.format_source(src_entry=matches[0], 
                                               mjd=mjd,
                                               magname=magname,
                                               designation=target_name,
                                           )

        print
        print "# Selected from file %s" % (filename)
        print "# Source magnitude & error: %.3f +/- %.3f" % (
            matches[0, SXcolumn['mag_aper_4.0']], matches[0, SXcolumn['mag_err_4.0']])
        print ret
        print

    # time.sleep(10)

    
