#!/usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
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
import podi_fitskybackground


scaling_factors = {
    "odi_g": 1.0,
    "odi_i": 0.15,
}

def subtract_pupilghost_extension(input_hdu, rotator_angle, filter, pupil_hdu, scaling, 
                                  rotate=True, verbose=True, non_negative=True):
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

    # Better: read the center coordinates from the pupil ghost template file
    if (filter in pupilghost_centers):
        if (extname in pupilghost_centers[filter]):
            center_x, center_y = pupilghost_centers[filter][extname]
        else:
            if (extname in pupilghost_centers):
                center_x, center_y = pupilghost_centers[extname]
            else:
                # Don't know what center coordinate to use, abort
                return
    elif (extname in pupilghost_centers):
        center_x, center_y = pupilghost_centers[extname]
    else:
        # Don't know what center coordinate to use, abort
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
        # Replace all NaNs with zeros
        template[numpy.isnan(template)] = 0.

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
    
    # Swap x/y since python does it too
    print "rot.shape=",rotated.shape
    bx = rotated.shape[0] / 2 - center_x
    by = rotated.shape[1] / 2 - center_y
    tx, ty = bx + data_shape[0], by + data_shape[1]

    correction = rotated[by:ty, bx:tx]

    if (non_negative):
        correction[correction < 0] = 0

    input_hdu.data -= (correction * scaling)

    input_hdu.header.update("PGCENTER", "%d %d" % (center_x, center_y))

    if (verbose): print "all done, going home!"
    return input_hdu




def subtract_pupilghost(input_hdu, pupil_hdu, scaling, rotate=True, verbose=True):
    """
    This is a wrapper around the subtract_pupilghost_extension routine,
    going through all the extensions that might have a pupil ghost.
    """
    
    # Add some parallelization here

    # Now go through each of the HDUs in the data frame and correct the pupil ghost
    for i in range(1, len(input_hdu)):
        rotator_angle = input_hdu[0].header['ROTSTART']
        filter = input_hdu[0].header['FILTER']
        subtract_pupilghost_extension(input_hdu[i], rotator_angle, filter, pupil_hdu, 
                                      rotate=rotate, scaling=scaling, verbose=verbose)

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


def get_pupilghost_scaling(science_frame, pupilghost_frame, 
                n_samples=750, boxwidth=20, 
                verbose=True):

    sci_hdulist = pyfits.open(science_frame) if (type(science_frame) == str) else science_frame
    pg_hdulist = pyfits.open(pupilghost_frame) if (type(pupilghost_frame) == str) else pupilghost_frame

    merged_all = numpy.zeros(shape=(0,6))

    for ext in range(1,5):
        extname = sci_hdulist[ext].header['EXTNAME']
        filter = sci_hdulist[0].header['FILTER']
        
        center = pupilghost_centers[filter][extname]

        data = sci_hdulist[ext].data
        pg = pg_hdulist[extname].data
        
        #n_samples = 1750
        samples = numpy.array(podi_fitskybackground.sample_background(data, None, None, 
                                                          min_found=n_samples, fit_regions=[], 
                                                          boxwidth=boxwidth))

        box_centers = samples[:,2:4]
        print box_centers.shape

        # Now convert the x/y from the science frame 
        # into x/y in the pupil ghost frame
        dxy = box_centers - center
        pg_xy = dxy + [ 0.5*pg.shape[0], 0.5*pg.shape[1] ]

        #print dxy[0:10,:]
        #print pg_xy[0:10,:]

        #print pg_samples.shape
        print pg_xy.shape
        pg_samples =  numpy.array(podi_fitskybackground.sample_background(pg, None, None, 
                                                                          fit_regions=[],
                                                             box_center=pg_xy))
#        pg_samples = podi_fitskybackground.sample_background(pg, None, None, 
#                                                             box_center=box_centers)

        print pg_samples.shape

        #continue

        ri = 600
        ro = 3000
        radius = numpy.sqrt(numpy.sum( dxy**2, axis=1))
        #print radius[0:10]

        print extname
        print center

        numpy.savetxt("xxx", samples)
        numpy.savetxt("yyy", pg_samples)

        # unaffected data points inside or outside the pupil ghost
        unaffected = (radius < ri) | (radius > ro)
        skylevel = numpy.median(samples[:,4][unaffected])
        print "skylevel=",skylevel
        print numpy.median(samples[:,4])

        # derive maximum intensity of pupil ghost 
        # better: read from file
        pg_max = numpy.max(pg_samples[:,4])

        # Now merge the two arrays:
        merged = numpy.empty(shape=(pg_xy.shape[0],6))
        merged[:,0:2] = samples[:,2:4]
        merged[:,2:4] = pg_samples[:,2:4]
        merged[:,4] = numpy.array(samples)[:,4] - skylevel
        merged[:,5] = numpy.array(pg_samples)[:,4]
        numpy.savetxt("zzz", merged)

        merged_all = numpy.append(merged_all, merged, axis=0)

        ratio = merged[:,4] / merged[:,5]
        
        strong_pg_signal = merged[:,5] > 0.4*pg_max

        valid_ratios = ratio[strong_pg_signal]

        median_ratio = numpy.median(valid_ratios)
        print "median_ratio =",median_ratio
        print "error = ",numpy.std(valid_ratios)
        numpy.savetxt("merged%d" % (ext), merged[strong_pg_signal])


    print "\n\n============================\n\n"
    ratio = merged_all[:,4] / merged_all[:,5]
        
    pg_max = numpy.max(merged_all[:,5])
    strong_pg_signal = merged_all[:,5] > 0.4*pg_max

    valid_ratios = ratio[strong_pg_signal]

    import podi_collectcells
    good_ratios = podi_collectcells.three_sigma_clip(ratio[strong_pg_signal])

    median_ratio = numpy.median(valid_ratios)
    print "median_ratio =",median_ratio
    print "error = ",numpy.std(valid_ratios)
    print "clipped median",numpy.median(good_ratios)
    print "clipped std",numpy.std(good_ratios)
    numpy.savetxt("merged_all", merged_all[strong_pg_signal])
    numpy.savetxt("merged_clipped", good_ratios)

    clipped_median = numpy.median(good_ratios)
    clipped_std = numpy.std(good_ratios)
    return clipped_median, clipped_std

if __name__ == "__main__":

    if (cmdline_arg_isset("-makeradial")):
        inputframe = get_clean_cmdline()[1]
        outputframe = get_clean_cmdline()[2]
        create_azimuthal_template(inputframe, outputframe)
        sys.exit(0)

    if (cmdline_arg_isset("-getscaling")):
        science_frame = get_clean_cmdline()[1]
        pupilghost_frame = get_clean_cmdline()[2]
        #get_scaling(science_frame, pupilghost_frame)

        sci_hdu = pyfits.open(science_frame)
        pg_hdu = pyfits.open(pupilghost_frame)
        get_pupilghost_scaling(sci_hdu, pg_hdu)
        sys.exit(0)

    if (cmdline_arg_isset("-cleanpupilghost")):
        science_frame = get_clean_cmdline()[1]
        pupilghost_frame = get_clean_cmdline()[2]
        output_frame = get_clean_cmdline()[3]
        #get_scaling(science_frame, pupilghost_frame)

        sci_hdu = pyfits.open(science_frame)
        pg_hdu = pyfits.open(pupilghost_frame)
        scaling, scale_std = get_pupilghost_scaling(sci_hdu, pg_hdu)
        subtract_pupilghost(sci_hdu, pg_hdu, scaling=scaling, rotate=False)
        
        sci_hdu.writeto(output_frame, clobber=True)

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

