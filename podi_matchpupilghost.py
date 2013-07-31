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

def subtract_pupilghost(input_hdu, pupil_hdu, scaling, rotate=True):

    if (rotate):
        # Get rotator angle from header
        rotator_angle = input_hdu[0].header['ROTSTART']
        # print "___%s___" % rotator_angle
        if (rotator_angle == "unknown"):
            rotator_angle = 0.0

        # Rotate to the right angle
        # Make sure to rotate opposite to the rotator angle since we rotate 
        # the correction and not the actual data
        template = pupil_hdu[1].data
        if (math.fabs(rotator_angle) < 0.5):
            rotated = template
        else:
            rotated = rotate_around_center(template, -1 * rotator_angle, mask_nans=False)
    else:
        template = pupil_hdu[1].data
        
    # Create a proper HDUList to hold the correction
    primhdu = pyfits.PrimaryHDU()
    hdulist = [primhdu]

    # Now go through each of the extensions for which we know the 
    # centers and cut out the right region of the template.
    for ext_name in pupilghost_centers:
        #print "Lookign for ___%s___" % (ext_name)
        
        # Now search for right OTA
        matched_hdu = None
        for i in range(1,len(input_hdu)):
            try:
                if (input_hdu[i].header['EXTNAME'] == ext_name):
                    matched_hdu = input_hdu[i]
                #print input_hdu[i].header['EXTNAME'],
            except:
                pass
        if (matched_hdu == None): 
            print "Could not find extension",ext_name,", this is a problem that shouldn't happen"
            continue

        # Check if this extension exists in the file to be corrected
        # Generally these OTAs should always exist, but better to be on the safe side
        data_shape = matched_hdu.data.shape
            
        #print ext_name, pupilghost_centers[ext_name]
        
        center_x, center_y = centers = pupilghost_centers[ext_name]
        # print ext_name, center_x, center_y

        # Swap x/y since python does it too
        print "rot.shape=",rotated.shape
        bx = rotated.shape[0] / 2 - center_x
        by = rotated.shape[1] / 2 - center_y
        tx, ty = bx + data_shape[0], by + data_shape[1]

        correction = rotated[by:ty, bx:tx]
        input_hdu[ext_name].data -= (correction * scaling)

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

    # Now, using the rotation angle given in the input frame, 
    # rotate and extract the pupil ghost sections for each of the OTAs.
    # And finally apply the correction
    scaling = 1.0
    if (len(sys.argv) > 4):
        scaling = float(sys.argv[4])
    print "using scaling factor",scaling
    hdu_matched = subtract_pupilghost(input_hdu, pupil_hdu, scaling)

    output_filename = sys.argv[3]
    hdu_matched.writeto(output_filename, clobber=True)

