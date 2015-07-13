#! /usr/bin/env python
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

import sys
import os
import pyfits
import numpy
#import ephem
import matplotlib.pyplot as plot
import scipy.interpolate
#import scipy.signal
import math
import scipy.optimize
import bottleneck
import dev_pgcenter

from podi_definitions import *
from podi_commandline import *
az_knot_limit = [50,600]

write_intermediate = True
use_buffered_files  = True


#
# Get the median intensity level in an annulus between r_inner and r_outer
#
def get_median_level(data, radii, ri, ro):

    selected = (radii > ri) & (radii < ro) #& (numpy.isfinite(data))
    pixelcount = numpy.sum(selected)
    if (pixelcount > 0):
        #cutout = data[selected]
        #median = numpy.median(cutout[0:5001])
        cutout = numpy.array(data[selected], dtype=numpy.float32)
        median = bottleneck.nanmean(cutout)
    else:
        median = numpy.NaN

    return median, pixelcount




def get_radii_angles(data_fullres, center, binfac, verbose=False):

    #
    # Rebin the image 4x to speed up calculations (the pupil ghost 
    # doesn't vary on small scales, so this is ok to do)
    #
    data = rebin_image(data_fullres, binfac)

    center_x, center_y = center
    
    #
    # Convert x/y coordinates to polar coordinates
    #
    if (verbose): stdout_write("   Computing radii ...\n")
    x, y = numpy.indices(data.shape)
    dx = x - center_x/binfac
    dy = y - center_y/binfac

    radius = numpy.sqrt(dx*dx + dy*dy)
    angle = numpy.arctan2(dx, dy)

    return data, radius, angle



def make_pupilghost_slice(filename, binfac, bpmdir, radius_range, clobber=False):

    hdu_ref = pyfits.open(filename)

    hdus = []
    centers = []

    rotator_angle = hdu_ref[0].header['ROTSTART'] 
    stdout_write("\nLoading frame %s ...\n" % (filename))

    #combined_file = "pg_combined_%+04d.fits" % numpy.around(rotator_angle)
    #print "combined-file:",combined_file

    output_filename = "pg_%+04d.fits" % (numpy.round(rotator_angle))
    if (os.path.isfile(output_filename) and not clobber):
        stdout_write("output filename %s already exists, skipping\n" % (output_filename))
        return None

    stdout_write("creating pupilghost slice %s ...\n" % (output_filename))
    
    datas = []
    extnames = []
    rotateds = []

    hdulist = [pyfits.PrimaryHDU()]

    for i in range(1, len(hdu_ref)):

        extname = hdu_ref[i].header['EXTNAME']
        
        pupilghost_centers = ['OTA33.SCI', 'OTA34.SCI', 'OTA43.SCI', 'OTA44.SCI']
        if (extname in pupilghost_centers):
            print "\n\n\n\n\n",extname

            #
            # Determine center position
            #
            # old method: use fixed values
            # center_x, center_y = pupilghost_centers[extname]

            input_hdu = hdu_ref[i]

            fx, fy, fr, vx, vy, vr = dev_pgcenter.find_pupilghost_center(input_hdu, verbose=False)
            center_x = vx
            center_y = vy


            #stdout_write("Using center position %d, %d for OTA %s\n" % (center_y, center_x, extname))
            stdout_write("Adding OTA %s, center @ %d, %d\n" % (extname, center_x, center_y))

            data = hdu_ref[i].data
            if (bpmdir != None):
                bpmfile = "%s/bpm_xy%s.reg" % (bpmdir, extname[3:5])
                mask_broken_regions(data, bpmfile, verbose=False)

            #hdus.append(hdu_ref[i])
            #centers.append((center_y, center_x))

            data = hdu_ref[i].data
            print "data-shape=",data.shape

            # Convert into radii and angles to make sure we can subtract the background
            binned, radius, angle = get_radii_angles(data, (center_y, center_x), binfac)
            print binned.shape, radius.shape, angle.shape, (center_y, center_x)

            # Fit and subtract the background
            bgsub = subtract_background(binned, radius, angle, radius_range, binfac)

            #
            # Insert this rawframe into the larger frame
            #
            combined = numpy.zeros(shape=(9000/binfac,9000/binfac), dtype=numpy.float32)
            combined[:,:] = numpy.NaN

            # Use center position to add the new frame into the combined frame
            # bx, by are the pixel position of the bottom left corner of the frame to be inserted
            bx = combined.shape[1] / 2 - center_x/binfac
            by = combined.shape[0] / 2 - center_y/binfac
            tx, ty = bx + bgsub.shape[0], by + bgsub.shape[1]
            print "insert target: x=", bx, tx, "y=", by, ty
            #combined[bx:tx, by:ty] = binned #bgsub[:,:]
            combined[by:ty, bx:tx] = bgsub[:,:]

            # Rotated around center to match the ROTATOR angle from the fits header
            combined_rotated = rotate_around_center(combined, rotator_angle, mask_nans=True, spline_order=1)

            imghdu = pyfits.ImageHDU(data=combined_rotated)
            imghdu.header.update('EXTNAME', extname)
            imghdu.header.update('ROTANGLE', rotator_angle)

            imghdu.header['OTA'] = int(extname[3:5])
            imghdu.header['PGCNTRFX'] = fx
            imghdu.header['PGCNTRFY'] = fy
            imghdu.header['PGCNTRFR'] = fr
            imghdu.header['PGCNTRVX'] = vx
            imghdu.header['PGCNTRVY'] = vy
            imghdu.header['PGCNTRVR'] = vr

            hdulist.append(imghdu)


    HDUlist = pyfits.HDUList(hdulist)
    HDUlist.writeto(output_filename, clobber=True)

    return rotator_angle


def subtract_background(data, radius, angle, radius_range, binfac):
    """
    This routine takes the input in polar coordinates and fits a straight line
    to the radial profile inside and outside of the allowed range. This is 
    assumed to be the background level (in analogy to the algorithm used in the 
    IRAF task mkpupil).

    Input data:
    - data (the actual intensity values for all pixels)
    - radius (the r in the polar coordianates)
    - angle (the phi in polar coordinates)
    - radius range (r_inner, r_outer, d_radius)
    - binfac (the binning used for the data)
    """

    # Compute the radial bin size in binned pixels
    print "subtracting background - binfac=",binfac
    r_inner, r_outer, dr_full = radius_range
    dr = dr_full/binfac
    r_inner /= binfac
    r_outer /= binfac

    #
    # Compute the number of radial bins
    #
    # Here: Add some correction if the center position is outside the covered area
    max_radius = 1.3 * r_outer #math.sqrt(data.shape[0] * data.shape[1])
    # Splitting up image into a number of rings
    n_radii = int(math.ceil(max_radius / dr))

    #
    # Compute the background level as a linear interpolation of the levels 
    # inside and outside of the pupil ghost
    #
    stdout_write("   Computing background-level ...")
    # Define the background ring levels
    radii = numpy.arange(0, max_radius, dr)
    background_levels = numpy.zeros(shape=(n_radii))
    background_level_errors = numpy.ones(shape=(n_radii)) * 1e9
    background_levels[:] = numpy.NaN
    for i in range(n_radii):

        ri = i * dr
        ro = ri + dr

        if (ri < r_inner):
            ro = numpy.min([ro, r_inner])
        elif (ro > r_outer):
            ri = numpy.max([ri, r_outer])
#        else:
#            # Skip the rings within the pupil ghost range for now
#            continue
        
        #print i, ri, ro
        median, count = get_median_level(data, radius, ri, ro)
        background_levels[i] = median
        background_level_errors[i] = 1. / math.sqrt(count) if count > 0 else 1e9

    # Now fit a straight line to the continuum, assuming it varies 
    # only linearly (if at all) with radius
    # define our (line) fitting function
    
    print "XXXXXXX", radii.shape, background_levels.shape
    numpy.savetxt("radial__%s" % ("x"),
                  numpy.append(radii.reshape((-1,1)),
                               background_levels.reshape((-1,1)), axis=1))
    print "saved"

    # Find average intensity at the largest radii
    avg_level = bottleneck.nanmedian(background_levels[radii>4000])
    print "avg_level=",avg_level

    #
    # Normalize profile
    #
    normalize_region = ((radii < 1100) & (radii > 600)) |  \
                       ((radii > 4000) & (radii < 4600))
    normalize_flux = numpy.mean(background_levels[normalize_region])
    print "normalization flux =", normalize_flux

    #
    # Subtract background and normalize all measurements
    #
    background_levels = (background_levels - normalize_flux) / normalize_flux

    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    bg_for_fit = background_levels
    #bg_for_fit[numpy.isnan(background_levels)] = 0
    bg_for_fit[((radii > ri) & (radii < ro))] = 0
    pinit = [0.0, 0.0] # Assume no slope and constant level of 0
    out = scipy.optimize.leastsq(errfunc, pinit,
                           args=(radii, background_levels, background_level_errors), full_output=1)

    pfinal = out[0]
    covar = out[1]
    stdout_write(" best-fit: %.2e + %.3e * x\n" % (pfinal[0], pfinal[1]))
    #print pfinal
    #print covar

    #
    # Now we have the fit for the background, compute the 2d background 
    # image and subtract it out
    #
    x = numpy.linspace(0, max_radius, 100)
    y_fit = radii * pfinal[1] + pfinal[0]
    background = pfinal[0] + pfinal[1] * radius
    
    bg_sub = ((data - normalize_flux) / normalize_flux) - background

    bg_sub_profile = background_levels - (pfinal[0] + pfinal[1]*radii)
    numpy.savetxt("radial__%s" % ("bgsub"),
                  numpy.append(radii.reshape((-1,1)),
                               bg_sub_profile.reshape((-1,1)), axis=1))
    
    #if (write_intermediate):
    #    bgsub_hdu = pyfits.PrimaryHDU(data=bg_sub)
    #    bgsub_hdu.writeto("bgsub.fits", clobber=True)
    
    return bg_sub



def create_radial_pupilghost(filename, outputfile, radial_opts, verbose=True):
    """
    This function takes a full multi-extension pupil ghost, 
    derives the azimuthal profile and writes it back to a new file.
    This azimuthal profile can then be used to remove the pupilghost
    from science with a range of rotator angles.

    Parameters are:
    - filename of the full 2-d input pupil ghost
    - name for the output file
    """
    
    hdulist = pyfits.open(filename)

    # open a text-file to hold the profile definition.
    profile_txt = open(outputfile+'.dat', 'w')

    # loop over all science exposures, skipping the primary header
    for ext in range(1, len(hdulist)):

        extname = hdulist[ext].header['EXTNAME']
        data_fullres = hdulist[ext].data
        center = (data_fullres.shape[0]/2, data_fullres.shape[1]/2)
        (r_inner, r_outer, dr) = radial_opts
        print radial_opts

        stdout_write("reading extension %d ..." % ext) #(hdulist[ext].header['EXTNAME']))
        # to cut down on computing, bin the frame
        binfac = 4
        data_binned, radius_binned, angle_binned = get_radii_angles(data_fullres, center, binfac)
        
        # For the output we need the polar coordinates in the full resolution
        dummy, radius_fullres, angle_fullres = get_radii_angles(data_fullres, center, 1)

        # Now split the frame into a number of radial rings
        # Use the binned data to keep computing time under control
        r_inner /= binfac
        r_outer /= binfac
        dr /= binfac

        # Allocate an array to hold the data for the pupil ghost, i.e.
        # all radii, average intensity levels, and results from the spline fit.
        n_rings = int(math.ceil((r_outer - r_inner) / dr))
        pupilghost_profile = numpy.zeros(shape=(n_rings,4))

        stdout_write(" computing profile ...")
        for r in range(n_rings):
            ri = r * dr + r_inner
            ro = ri + dr

            in_this_ring = (radius_binned >= ri) & (radius_binned < ro)

            pupilghost_profile[r,0] = (ri+ro)/2 #math.sqrt((ri**2 + ro**2)/2)
            pupilghost_profile[r,1] = numpy.median(data_binned[in_this_ring])

        # Make sure there are no pixels with negative values.
        pupilghost_profile[:,2] = pupilghost_profile[:,1]
        pupilghost_profile[:,2][pupilghost_profile[:,2] < 0] = 0.

        #
        # Now fit the profile with a 1-D spline
        # limit the number of knots to 75, otherwise we'll run into a 
        # weird bug in the spline fitting
        #
        n_knots = (r_outer-r_inner-2*dr)/dr-1
        if (n_knots > 75): 
            n_knots=75
        radial_knots = numpy.linspace(r_inner+0.7*dr, r_outer-0.7*dr, n_knots)
        if (verbose): print "radial knots=",radial_knots[0:5],"...",radial_knots[-5:]

        stdout_write(" fitting ...")
        radial_profile = scipy.interpolate.LSQUnivariateSpline(
            pupilghost_profile[:,0], pupilghost_profile[:,2],
            radial_knots, k=2)
  
        stdout_write(" computing output ...")
        # Now that the fitting is done, compute the spline fit at the 
        # original positions so we can compare things
        pupilghost_profile[:,3] = radial_profile(pupilghost_profile[:,0])

        print >>profile_txt, "#", extname
        numpy.savetxt(profile_txt, pupilghost_profile)
        print >>profile_txt, "\n\n\n\n\n"

        #
        # Compute the 2-d radial profile
        #
        radius_fullres_asbinned = radius_fullres / binfac
        radius_1d = radius_fullres_asbinned.ravel()
        print "rad-1d",radius_1d.shape

        radial_pupilghost = radial_profile(radius_1d).reshape(radius_fullres.shape)
        print "rad pg",radial_pupilghost.shape

        # set all pixels outside the pupil ghost radial range to 0
        radial_pupilghost[(radius_fullres_asbinned > r_outer) | (radius_fullres_asbinned < r_inner)] = 0

        # and save the pupil ghost
        hdulist[ext].data = radial_pupilghost
        stdout_write(" done!\n")

    # Now we are done with all profiles, write the results to the output file
    clobberfile(outputfile)
    stdout_write("writing output file ...")
    hdulist.writeto(outputfile, clobber=True)
    stdout_write(" done!\n")


#################################
#
# Important note:
# Many x/y values are swapped, because fits data is arranged in y/x coordinates, not x/y
#
#################################
if __name__ == "__main__":

    print """\

  fit-pupilghost tool
  part of the pODI QuickReduce pipeline
  (c) 2013, Ralf Kotulla (kotulla@uwm.edu)
  Contact the author with questions and comments.
  Website: members.galev.org/rkotulla
           --> Research --> podi-pipeline
"""

    # Read in the input parameters
    binfac = int(cmdline_arg_set_or_default("-prebin", 1))
    bpmdir = cmdline_arg_set_or_default("-bpm", None)

    if (cmdline_arg_isset("-radial")):
        r_inner = float(cmdline_arg_set_or_default("-ri",  700))
        r_outer = float(cmdline_arg_set_or_default("-ro", 4500))
        dr = float(cmdline_arg_set_or_default("-dr", 20))

        filename = get_clean_cmdline()[1]
        outputfile = get_clean_cmdline()[2]
        radial_opts = (r_inner, r_outer, dr)

        create_radial_pupilghost(filename, outputfile, radial_opts)
    else:
        r_inner = float(cmdline_arg_set_or_default("-ri",  700))
        r_outer = float(cmdline_arg_set_or_default("-ro", 4000))
        dr = float(cmdline_arg_set_or_default("-dr", 20))

        filenames = get_clean_cmdline()[1:]
        print filenames

        radius_range = (r_inner, r_outer, dr)

        for inputfile in filenames:
            make_pupilghost_slice(inputfile, binfac, bpmdir, radius_range, clobber=False)

    sys.exit(0)

