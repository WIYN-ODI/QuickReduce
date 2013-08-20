#! /usr/bin/env python

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

from podi_definitions import *
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




def get_radii_angles(data_fullres, center, binfac):

    #
    # Rebin the image 4x to speed up calculations (the pupil ghost 
    # doesn't vary on small scales, so this is ok to do)
    #
    data = rebin_image(data_fullres, binfac)

    center_x, center_y = center
    
    #
    # Convert x/y coordinates to polar coordinates
    #
    stdout_write("   Computing radii ...\n")
    x, y = numpy.indices(data.shape)
    dx = x - center_x/binfac
    dy = y - center_y/binfac

    radius = numpy.sqrt(dx*dx + dy*dy)
    angle = numpy.arctan2(dx, dy)

    return data, radius, angle


def merge_OTAs(hdus, centers):

    combined = numpy.zeros(shape=(9000,9000), dtype=numpy.float32)
    combined[:,:] = numpy.NaN

    stdout_write("   Adding OTA")
    for i in range(len(hdus)):

        #stdout_write("Adding OTA %s ...\n" % (hdus[i].header['EXTNAME']))
        stdout_write(" %s" % (hdus[i].header['EXTNAME']))
        # Use center position to add the new frame into the combined frame
        # bx, by are the pixel position of the bottom left corner of the frame to be inserted
        bx = combined.shape[0] / 2 - centers[i][0]
        by = combined.shape[1] / 2 - centers[i][1]
        #print bx, by
        tx, ty = bx + hdus[i].data.shape[0], by + hdus[i].data.shape[1]

        combined[by:ty, bx:tx] = hdus[i].data[:,:]
        #combined[bx:tx, by:ty] = hdus[i].data[:,:]
    stdout_write("\n")

    return combined


def make_pupilghost_slice(filename, binfac, bpmdir, radius_range, clobber=False):

    hdu_ref = pyfits.open(filename)

    hdus = []
    centers = []

    rotator_angle = hdu_ref[0].header['ROTSTART'] 
    stdout_write("\nLoading frame %s ...\n" % (filename))

    #combined_file = "pg_combined_%+04d.fits" % numpy.around(rotator_angle)
    #print "combined-file:",combined_file

    output_filename = "pg_%+04d.fits" % (rotator_angle)
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

        if (extname in pupilghost_centers):
            print "\n\n\n\n\n",extname
            center_x, center_y = pupilghost_centers[extname]

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
            
            hdulist.append(imghdu)

    HDUlist = pyfits.HDUList(hdulist)
    HDUlist.writeto(output_filename, clobber=True)

    return rotator_angle


def subtract_background(data, radius, angle, radius_range, binfac):

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
        else:
            # Skip the rings within the pupil ghost range for now
            continue
        
        #print i, ri, ro
        median, count = get_median_level(data, radius, ri, ro)
        background_levels[i] = median
        background_level_errors[i] = 1. / math.sqrt(count) if count > 0 else 1e9

    # Now fit a straight line to the continuum, assuming it varies 
    # only linearly (if at all) with radius
    # define our (line) fitting function
    
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    bg_for_fit = background_levels
    bg_for_fit[numpy.isnan(background_levels)] = 0
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
    
    bg_sub = data - background

    #if (write_intermediate):
    #    bgsub_hdu = pyfits.PrimaryHDU(data=bg_sub)
    #    bgsub_hdu.writeto("bgsub.fits", clobber=True)
    
    return bg_sub



    #
    # Now we have a collection of a bunch of files, possibly each with separate rotator angles
    #

    hdu = pyfits.PrimaryHDU(data = all_data)
    hdu.writeto("all_data.fits", clobber=True)
    hdu = pyfits.PrimaryHDU(data = all_bgsub)
    hdu.writeto("all_bgsub.fits", clobber=True)

    pupil_sub, radial_profile, radial_2d = fit_radial_profile(all_data, all_radius, all_angle, all_bgsub, radius_range)

    pyfits.PrimaryHDU(data=pupil_sub).writeto("pupilsub.fits", clobber=True)
    pyfits.PrimaryHDU(data=radial_2d).writeto("radial2d.fits", clobber=True)
    # create_mapped_coordinates(all_data, all_radius, all_angle, all_bgsub, pupil_sub, radius_range, binfac)

    azimuthal_fits = fit_azimuthal_profiles(all_data, all_radius, all_angle, all_bgsub, pupil_sub, radius_range)

    
    # 
    # Now all the fitting is done, let's compute the output
    #

    # First get a fresh buffer of coordinates
    outbuffer = numpy.zeros(shape=(9000,9000))
    out_data, out_radius, out_angle = get_radii_angles(outbuffer, (outbuffer.shape[0]/2, outbuffer.shape[1]/2), binfac)
    out_angle[out_angle < 0] += 2*numpy.pi

    azimuthal_2d = compute_pupilghost(out_data, out_radius, out_angle, radius_range, binfac,
                                      azimuthal_fits)
    
    pyfits.PrimaryHDU(data=azimuthal_2d).writeto("fit_nonradial.fits", clobber=True)

    # Compute the 2-d radial profile. The extreme values beyond the fitting radius 
    # might be garbage, so set all pixels outside the pupil ghost radial range to 0
    radial_2d = radial_profile(out_radius.ravel()).reshape(out_radius.shape)
    radial_2d[(radius > r_outer/binfac) | (radius < r_inner/binfac)] = 0

    pyfits.PrimaryHDU(data=radial_2d).writeto("fit_radial.fits", clobber=True)
    try:
        full_2d = azimuthal_2d + radial_2d
        full_2d[full_2d<0] = 0

        print "Writing data"
        pyfits.PrimaryHDU(data=full_2d).writeto("fit_rad+nonrad.fits", clobber=True)
    except:
        pass

    #leftover = bg_sub - fullprofile
    #    pyfits.PrimaryHDU(data=leftover).writeto("fit_leftover.fits", clobber=True)


    return

    #------------------------------------------------------------------------------
    #
    # Until now the template is still binned, blow it up to the full resolution
    #
    #------------------------------------------------------------------------------


    print "Interpolating to full resolution"
    xb, yb = numpy.indices(data.shape)
    
    # Prepare the 2-d interpolation spline
    interpol = scipy.interpolate.RectBivariateSpline(xb[:,0], yb[0,:], fullprofile)

    # And use above spline to compute the full-resolution version
    xo, yo = numpy.indices(data_fullres.shape)
    xo = xo * 1.0 / data_fullres.shape[0] * data.shape[0]
    yo = yo * 1.0 / data_fullres.shape[1] * data.shape[1]
    correction = interpol(xo[:,0], yo[0,:]).reshape(data_fullres.shape)

    return correction



    return



def fit_radial_profile(data, radius, angle, bgsub, radius_range, binfac=1, verbose=False, show_plots=False,
                       force_positive=False, zero_edges=False, save_profile=None):
    

    #------------------------------------------------------------------------------
    #
    # Here we assume that all files have their background removed. 
    # Then we continue with the radial profile.
    #
    #------------------------------------------------------------------------------

    r_inner, r_outer, dr_full = radius_range
    dr = dr_full/binfac
    r_inner /= binfac
    r_outer /= binfac

    n_radii = int(math.ceil(r_outer / dr))

    stdout_write("\nComputing radial profile (ri=%d, ro=%d, n=%d, bin=%d)...\n" % (r_inner, r_outer, n_radii, binfac))

    #
    # Next step: Fit the radial profile
    #
    pupil_radii = numpy.zeros(shape=(n_radii))
    pupil_level = numpy.zeros(shape=(n_radii))
    min_r, max_r = 10000, -10000
    for i in range(n_radii):

        ri = i * dr
        ro = ri + dr

        # Only use rings within the pupil gost rings
        if (ri < r_outer and ro > r_inner):
            ri = numpy.max([ri, r_inner])
            ro = numpy.min([ro, r_outer])
            min_r = numpy.min([min_r, i])
            max_r = numpy.max([max_r, i])
        else:
            continue

        if (verbose): stdout_write("radius i=%4d" % (i))
        median, count = get_median_level(bgsub, radius, ri, ro)

        pupil_radii[i] = 0.5 * (ri + ro)
        pupil_level[i] = median
        if (verbose): stdout_write("   ri: %4d    ro: %4d    med: %.4f\n" % (ri, ro, median))

    if (force_positive):
        pupil_level[pupil_level < 0] = 0.
    if (zero_edges):
        pupil_level[min_r] = 0.
        pupil_level[max_r] = 0.

    if (save_profile != None):
        out = open(save_profile, "w")
        print >>out, "# Data: radius, median"
        dummy = numpy.empty(shape=(pupil_radii.shape[0],2))
        dummy[:,0] = pupil_radii[:]
        dummy[:,1] = pupil_level[:]
        #dummy = numpy.append(pupil_radii, pupil_level, axis=1)
        numpy.savetxt(out, dummy)

    #
    # Now fit the profile with a 1-D spline
    #
    n_knots = (r_outer-r_inner-2*dr)/dr-1
    if (n_knots > 75): n_knots=75
    radial_knots = numpy.linspace(r_inner+0.7*dr, r_outer-0.7*dr, n_knots)
    if (verbose): print "radial knots=",radial_knots[0:5],"...",radial_knots[-5:]
    radial_profile = scipy.interpolate.LSQUnivariateSpline(pupil_radii[min_r:max_r+1], pupil_level[min_r:max_r+1], radial_knots, k=2)
    
    # In case of ValueError, this is what causes it (from scipy/intrpolate/interpolate.py): 
    # if not alltrue(t[k+1:n-k]-t[k:n-k-1] > 0,axis=0):
    #     raise ValueError('Interior knots t must satisfy '
    #                      'Schoenberg-Whitney conditions')

    if (save_profile != None):
        print >>out, "\n\n\n\n\n\n\n\n"

        print >>out, "# spline fit: radius, level"
        smooth_radial_1d_x = numpy.linspace(r_inner, r_outer, 1300)
        smooth_radial_1d_y = radial_profile(smooth_radial_1d_x)
        dummy = numpy.empty(shape=(smooth_radial_1d_x.shape[0],2))
        dummy[:,0] = smooth_radial_1d_x[:]
        dummy[:,1] = smooth_radial_1d_y[:]
        numpy.savetxt(out, dummy)

        out.close()


    if (show_plots):
        # create a smooth profile for plotting
        smooth_radial_1d_x = numpy.linspace(r_inner, r_outer, 1300)
        smooth_radial_1d_y = radial_profile(smooth_radial_1d_x)
        knots_y = radial_profile(radial_knots)
        plot.plot(pupil_radii, pupil_level, '.-', radial_knots, knots_y, 'o', smooth_radial_1d_x, smooth_radial_1d_y)
        plot.show()

    #
    # Compute the 2-d radial profile
    #
    radius_1d = radius.ravel()
    pupil_radial_2d = radial_profile(radius.ravel()).reshape(radius.shape)
    # set all pixels outside the pupil ghost radial range to 0
    pupil_radial_2d[(radius > r_outer) | (radius < r_inner)] = 0
    # and subtract the pupil ghost
    
    pupil_sub = bgsub - pupil_radial_2d

    #if (write_intermediate):
    pupil_sub_hdu = pyfits.PrimaryHDU(data = pupil_sub)
    pupil_sub_hdu.writeto("all_pupilsub.fits", clobber=True)

    return pupil_sub, radial_profile, pupil_radial_2d

#    template_radius_1d = template_radius.ravel()
#    template_radial = radial_profile(template_radius.ravel()).reshape(template_radius.shape)
#    template_radial[(template_radius > r_outer) | (template_radius < r_inner)] = 0
#    pupil_sub_hdu = pyfits.PrimaryHDU(data = template_radial)
#    pupil_sub_hdu.writeto("template_radial.fits", clobber=True)


    return




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
    r_inner = float(cmdline_arg_set_or_default("-ri",  700))
    r_outer = float(cmdline_arg_set_or_default("-ri", 4000))
    dr = float(cmdline_arg_set_or_default("-dr", 20))
    binfac = int(cmdline_arg_set_or_default("-prebin", 4))
    bpmdir = cmdline_arg_set_or_default("-bpm", None)

    filenames = get_clean_cmdline()[1:]

    radius_range = (r_inner, r_outer, dr)

    for inputfile in filenames:
        make_pupilghost_slice(inputfile, binfac, bpmdir, radius_range, clobber=False)

    sys.exit(0)

