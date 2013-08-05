#!/usr/bin/env python

#
# (c) Ralf Kotulla for WIYN/pODI
#

"""
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy
import scipy

gain_correct_frames = False
from podi_definitions import *


scaling_factors = {
    "odi_g": 1.0,
    "odi_i": 1.0,
}

def subtract_pupilghost_extension(input_hdu, rotator_angle, pupil_hdu, scaling, rotate=True, verbose=True):
    """
    Contains all functionality to dexecute the actual pupil ghost removal. 
    The pupil ghost is roatted to match the rotation angle of the science 
    frame, then scaled with the specified scaling factor and lastly the 
    scaled template is removed from the science frame. The results stay in 
    the input_hdu variable that is also returned.
    """

    #
    # Check if this extension is included in the pupilghost template
    #
    extname = input_hdu.header['EXTNAME']
    right_ext = -1
    for this_ext in range(1, len(pupil_hdu)):
        if (pupil_hdu[this_ext].header['EXTNAME'] == extname):
            right_ext = this_ext
            break

    if (right_ext < 0):
        # This extension is not listed in the pupil-ghost template
        # We can't do anythin in this case, so let's simply return
        return

    #
    # Now we know we have to do something
    #
    if (verbose): print "Found matching pupil ghost template in extension",right_ext

    if (rotate):
        # print "___%s___" % rotator_angle
        if (rotator_angle == "unknown"):
            rotator_angle = 0.0

        # Rotate to the right angle
        # Make sure to rotate opposite to the rotator angle since we rotate 
        # the correction and not the actual data
        template = pupil_hdu[right_ext].data
        if (math.fabs(rotator_angle) < 0.5):
            if (verbose): print "Rotator angle is small (%.2f), skipping rotation" % (rotator_angle)
            rotated = template
        else:
            if (verbose): print "rotating template by",rotator_angle,"degrees"
            rotated = rotate_around_center(template, -1 * rotator_angle, mask_nans=False)
    else:
        if (verbose): print "No rotation requested, skipping rotation"
        rotated = pupil_hdu[right_ext].data
        
    #
    # Ok, now we have the pupil ghost rotated to the right angle
    #
        
    data_shape = input_hdu.data.shape
    
    # Better: read the center coordinates from the pupil ghost template file
    center_x, center_y = pupilghost_centers[extname]

    # Swap x/y since python does it too
    print "rot.shape=",rotated.shape
    bx = rotated.shape[0] / 2 - center_x
    by = rotated.shape[1] / 2 - center_y
    tx, ty = bx + data_shape[0], by + data_shape[1]

    correction = rotated[by:ty, bx:tx]
    input_hdu.data -= (correction * scaling)

    if (verbose): print "all done, going home!"
    return input_hdu




def subtract_pupilghost(input_hdu, pupil_hdu, scaling, verbose=True):
    """
    This is a wrapper around the subtract_pupilghost_extension routine,
    going through all the extensions that might have a pupil ghost.
    """
    
    # Add some parallelization here

    # Now go through each of the HDUs in the data frame and correct the pupil ghost
    for i in range(1, len(input_hdu)):
        rotator_angle = input_hdu[0].header['ROTSTART']
        subtract_pupilghost_extension(input_hdu[i], rotator_angle, pupil_hdu, scaling=scaling)

    return input_hdu



def create_azimuthal_template(filename, outputfilename):
    """
    This routine will read the full pupilghost template and create 
    a radial profile by averaging out the azimuthal variations. This 
    radial pupil template will then be used in collectcells to remove 
    the pupil ghost from science frames.
    """

    # Open file and extract data
    # Things here are simple, since the file should only have the primary 
    # and one image extension
    hdulist = pyfits.open(filename)
    raw_data = hdulist[1].data

    binning = 1
    import podi_fitpupilghost
    data, radius, angle = podi_fitpupilghost.get_radii_angles(raw_data, (4500,4500), binfac=binning)

    pupilsub, radial_profile, pupil_radial_2d = \
        podi_fitpupilghost.fit_radial_profile(data, radius, angle, data, (500, 4180, 10), binfac=binning, 
                                              verbose=True, show_plots=True,
                                              force_positive=False, zero_edges=False,
                                              save_profile=outputfilename[:-5]+".profile.dat")

    # pyfits.writeto("radial_profile.fits", radial_profile)
    # pyfits.writeto("pupilsub.fits", pupilsub, clobber=True)
    # pyfits.writeto("pupil_radial_2d.fits", pupil_radial_2d, clobber=True)

    hdulist[1].data = pupil_radial_2d
    
    clobberfile(outputfilename)
    hdulist.writeto(outputfilename)

    return


if __name__ == "__main__":

    if (cmdline_arg_isset("-makeradial")):
        inputframe = get_clean_cmdline()[1]
        outputframe = get_clean_cmdline()[2]
        create_azimuthal_template(inputframe, outputframe)
        sys.exit(0)

    # Read filenames from command line
    inputframe = sys.argv[1]
    pupilghost_template = sys.argv[2]

    # and open all fits files
    input_hdu = pyfits.open(inputframe)
    pupil_hdu = pyfits.open(pupilghost_template)

    scaling = 1.0
    if (len(sys.argv) > 4):
        scaling = float(sys.argv[4])
    print "using scaling factor",scaling

    subtract_pupilghost(input_hdu, pupil_hdu, scaling=scaling)

    # Now, using the rotation angle given in the input frame, 
    # rotate and extract the pupil ghost sections for each of the OTAs.
    # And finally apply the correction
    #hdu_matched = subtract_pupilghost(input_hdu, pupil_hdu, scaling)


    output_filename = sys.argv[3]
    input_hdu.writeto(output_filename, clobber=True)

