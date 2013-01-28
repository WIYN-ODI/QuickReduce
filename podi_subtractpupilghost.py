#! /usr/bin/env python

import sys
import os
import pyfits
import numpy

from podi_definitions import *

import scipy.interpolate


def subtract_pupilghost(filename_hdu, output_file, profile_file, operation="subtract", verbose=True):

    # Make the routine flexible in handling both a filename or HDUList as input parameter
    if (type(filename_hdu) is str):
        if (os.path.isfile(fitsfile)):
            hdulist = pyfits.open(fitsfile)
        else:
            return
    else:
        hdulist = filename_hdu

    filter = hdulist[0].header['FILTER']

    if (os.path.isfile(profile_file) and not operation=="ignore"):
        profile_data = numpy.loadtxt(profile_file)

        # Loop over all extensions
        # For now only use the first one, hard enough
        for ota_id in range(1, len(hdulist)):

            ota = hdulist[ota_id]
            extname = ota.header["EXTNAME"]

            if (filter in pupilghost_centers):
                if (extname in pupilghost_centers[filter]):
                    center_x, center_y = pupilghost_centers[filter][extname]

                    if (verbose): stdout_write("\rCorrecting pupil-ghost in extension %s ..." % (extname))
                    # Create the pixel array and create all distances
                    #        buffer = numpy.zeros(shape=(ota.data.shape[1], ota.data.shape[0]))
                    x, y = numpy.indices(ota.data.shape)                                                                                                           
                    dx = x - center_x                                                                                                                           
                    dy = y - center_y                                                                                                                            
                    d = numpy.sqrt(dx*dx + dy*dy)

                    # Next load the profile, and fit a spline that we can use to evaluate it
                    profile = scipy.interpolate.interp1d(profile_data[:,0]*4, profile_data[:,1], bounds_error=False, fill_value=0)
                    fit = profile(d.ravel())
                    
                    if (operation == "replace"):
                        hdulist[ota_id].data = fit.reshape(ota.data.shape).transpose()
                    else:
                        hdulist[ota_id].data -= fit.reshape(ota.data.shape).transpose()

                    ota.header.add_history("Pupilghost: profile in %s" % profile_file)
                    ota.header.add_history("Pupilghost: center @ %d,%d" % (center_x, center_y))

    if (not output_file == None):
        # Now all the work is done, all final data is stored in img2, write results to new file                    
        if (verbose): stdout_write(" writing output ...")
        clobberfile(outputfile)
        hdulist.writeto(outputfile)
        if (verbose): stdout_write(" done!\n\n")
    else:
        if (verbose): stdout_write(" done, handing on data!\n\n")
        return hdulist
    

if __name__ == "__main__":

    # Read in the input parameters
    fitsfile = sys.argv[1]
    profile_file = sys.argv[2]
    outputfile = sys.argv[3]

    operation = cmdline_arg_set_or_default("-op", "subtract")

    subtract_pupilghost(fitsfile, outputfile, profile_file, operation=operation, verbose=True)
