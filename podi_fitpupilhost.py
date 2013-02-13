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

from podi_definitions import *


def add_circle(buffer, center_x, center_y, radius, amplitude):

    x, y = numpy.indices(buffer.shape)

    dx = x - center_x
    dy = y - center_y
    d2 = dx*dx + dy*dy

    print dx[0:10,0:10]
    print dy[0:10,0:10]
    print d2[0:10,0:10]

    #print d2[995:1005, 995:1005]
    
    tmp_buffer = numpy.zeros(shape=buffer.shape)
    tmp_buffer[d2 < radius*radius] = amplitude

    buffer += tmp_buffer

    return

def add_annulus(buffer, center_x, center_y, radius_i, radius_o, amplitude):

    x, y = numpy.indices(buffer.shape)

    dx = x - center_x
    dy = y - center_y
    d2 = dx*dx + dy*dy

    #print dx[0:10,0:10]
    #print dy[0:10,0:10]
    #print d2[0:10,0:10]

    #print d2[995:1005, 995:1005]
    
    tmp_buffer = numpy.zeros(shape=buffer.shape)

    selected_pixels = (d2 < radius_o*radius_o) & (d2 > radius_i*radius_i)
    tmp_buffer[selected_pixels] = amplitude

    buffer += tmp_buffer

    return

def create_from_template(command_file, buffer):
    
    # Load command file 
    cmdfile = open(command_file, "r")
    cmds = cmdfile.readlines()
    print cmds
    
    for i in range(len(cmds)):
        line = cmds[i]
        if (line[0] == "#"):
            continue
        
        items = line.strip().split()
        print items

        shape = items[0]
        if (shape == "fillcircle"):
            center_x = float(items[1])
            center_y = float(items[2])
            radius = float(items[3])
            amplitude = float(items[4])
            
            add_circle(buffer, center_x, center_y, radius, amplitude)

        if (shape == "annulus"):
            center_x = float(items[1])
            center_y = float(items[2])
            radius_i = float(items[3])
            radius_o = float(items[4])
            amplitude = float(items[5])
            
            add_annulus(buffer, center_x, center_y, radius_i, radius_o, amplitude)

            
    return buffer


#
# Get the median intensity level in an annulus between r_inner and r_outer
#
def get_median_level(data, radii, ri, ro):

    selected = (radii > ri) & (radii < ro) & (numpy.isfinite(data))
    pixelcount = numpy.sum(selected)
    if (pixelcount > 0):
        median = numpy.median(data[selected])
    else:
        median = numpy.NaN

    return median, pixelcount


def fit_pupilghost(data_fullres, center, radius_range, dr_full, 
                   write_intermediate=False, show_plots=False):


    #
    # Rebin the image 4x to speed up calculations (the pupil ghost 
    # doesn't vary on small scales, so this is ok to do)
    #
    binfac = 4
    dr = dr_full/binfac
    data = rebin_image(data_fullres, binfac)

    center_x, center_y = center
    
    r_inner, r_outer = radius_range
    r_inner /= binfac
    r_outer /= binfac

    #
    # Convert x/y coordinates to polar coordinates
    #
    print "Computing radii"
    x, y = numpy.indices(data.shape)
    dx = x - center_x/binfac
    dy = y - center_y/binfac

    radius = numpy.sqrt(dx*dx + dy*dy)
    angle = numpy.arctan2(dx, dy)

    # write some intermediate data products
    if (write_intermediate):
        raw_hdu = pyfits.PrimaryHDU(data=data)
        raw_hdu.writeto("raw.fits", clobber=True)

    #
    # Compute the number of radial bins
    #
    # Here: Add some correction if the center position is outside the covered area
    max_radius = 1.3 * math.sqrt(data.shape[0] * data.shape[1])
    # Splitting up image into a number of rings
    n_radii = int(math.ceil(max_radius / dr))


    #
    # Compute the background level as a linear interpolation of the levels 
    # inside and outside of the pupil ghost
    #
    print "Computing background-level ..."
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
    pinit = [1.0, 0.0] # Assume no slope and constant level of 1
    out = scipy.optimize.leastsq(errfunc, pinit,
                           args=(radii, background_levels, background_level_errors), full_output=1)

    pfinal = out[0]
    covar = out[1]
    print "Best-fit: %.2f + %f * x" % (pfinal[0], pfinal[1])
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

    if (write_intermediate):
        bgsub_hdu = pyfits.PrimaryHDU(data=bg_sub)
        bgsub_hdu.writeto("bgsub.fits", clobber=True)
    

    #------------------------------------------------------------------------------
    #
    # Now we got rid of background gradients, continue with the radial profile
    #
    #------------------------------------------------------------------------------


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

        median, count = get_median_level(bg_sub, radius, ri, ro)

        pupil_radii[i] = 0.5 * (ri + ro)
        pupil_level[i] = median


    #
    # Now fit the profile with a 1-D spline
    #
    radial_knots = numpy.linspace(r_inner+0.7*dr, r_outer-0.7*dr, (r_outer-r_inner-1.4*dr)/dr/1.5)
    radial_profile = scipy.interpolate.LSQUnivariateSpline(pupil_radii[min_r:max_r+1], pupil_level[min_r:max_r+1], radial_knots, k=2)
    
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
    
    pupil_sub = bg_sub - pupil_radial_2d

    if (write_intermediate):
        pupil_sub_hdu = pyfits.PrimaryHDU(data = pupil_sub)
        pupil_sub_hdu.writeto("pupilsub.fits", clobber=True)


    #------------------------------------------------------------------------------
    #
    # Now we have only the non-radial components of the pupil ghost left
    #
    #------------------------------------------------------------------------------


    radial_splines = [None] * n_radii
    for cur_radius in range(n_radii):

        ri = cur_radius * dr
        ro = ri + dr

        # Make sure we ignore points outside the pupil ghost area
        if (ri < r_outer and ro > r_inner):
            ri = numpy.max([ri, r_inner])
            ro = numpy.min([ro, r_outer])
            min_r = numpy.min([min_r, cur_radius])
            max_r = numpy.max([max_r, cur_radius])
        else:
            continue

        # Keep user informed and happy
        sys.stdout.write("\rFitting azimuthal profile in ring %d - %d" % (int(ri), int(ro)))
        sys.stdout.flush()

        # Get all pixels in this ring
        pixels_in_ring = (radius >= ri) & (radius < ro)
        ring_angles = angle[pixels_in_ring]
        ring_data   = pupil_sub[pixels_in_ring]

        # Eliminate all NaN pixels
        valid_pixels_in_ring = numpy.isfinite(ring_data)

        #
        # Get min and max angles
        #
        min_angle = numpy.min(ring_angles[valid_pixels_in_ring])
        max_angle = numpy.max(ring_angles[valid_pixels_in_ring])

        valid_angles = ring_angles[valid_pixels_in_ring]
        valid_data   = ring_data[valid_pixels_in_ring]
        mean_da = (max_angle - min_angle) / valid_angles.shape[0]

        #max_points = 50 #numpy.min([ring_angles.shape[0], 20])
        #angle_knots = numpy.linspace(min_angle+0.01, max_angle-0.01, max_points)


        number_knots = numpy.sum(valid_pixels_in_ring) / 100
        if (number_knots < 30): number_knots = 30
        angle_knots = numpy.linspace(min_angle, max_angle, number_knots) 

        #print "\nmean angle between 2 points:", mean_da

        # Now eliminate knots in regions with no data
        for i in range(angle_knots.shape[0]):
            diff_angle = numpy.fabs(angle_knots[i] - valid_angles)
            min_diffangle = numpy.min(diff_angle)
            if (min_diffangle > 3 * mean_da or angle_knots[i] <= min_angle or angle_knots[i] >= max_angle):
                angle_knots[i] = numpy.NaN

        #print angle_knots.shape[0],"good points out of",number_knots

        # sort all points in this ring, otherwise LSQUnivariateSplie chokes
        si = numpy.argsort(valid_angles)
        sorted_angles = numpy.zeros(valid_angles.shape)
        sorted_data   = numpy.zeros(valid_data.shape)
        for i in range(si.shape[0]):
            sorted_angles[i] = valid_angles[si[i]]
            sorted_data[i]   = valid_data[si[i]]

        # Fit spline
        good_angle_knots = angle_knots[numpy.isfinite(angle_knots)]
        #print good_angle_knots.shape[0],"good points out of",angle_knots.shape[0]
        #print "After sorting", sorted_angles.shape[0], sorted_data.shape[0]
        
        az_profile = scipy.interpolate.LSQUnivariateSpline(sorted_angles, sorted_data, good_angle_knots, k=3)
        

        # Bin the ring in XXX degree boxes

        if (False):
            az_binwidth = math.radians(2)
            binned_angles, binned_data, binned_scatter = bin_azimuthal_profile(ring_angles, ring_data, 
                                                                               az_binwidth, min_angle, max_angle)


            # test 1
            angle_knots = numpy.linspace(binned_angles[1], binned_angles[-2], 20) 

            #
            # Fit a chebychev polynomial to the az data in the given ring
            #
            weight = 1./ (binned_scatter**2)
            cheb_coeff = numpy.polynomial.chebyshev.chebfit(binned_angles, binned_data, deg=9, w=weight)

            az_profile = scipy.interpolate.LSQUnivariateSpline(binned_angles, binned_data, angle_knots, w=weight, k=2)

        fine_profile_x = numpy.linspace(min_angle, max_angle, 1000)
        #fine_profile_y = numpy.polynomial.chebyshev.chebval(fine_profile_x, cheb_coeff)

        radial_splines[cur_radius] = az_profile

        #fine_profile_x = numpy.linspace(min_angle, max_angle, 1000)
        #fine_profile_x = numpy.linspace(binned_angles[0], binned_angles[-1], 1000)
        #az_profile.set_smoothing_factor(0.2)
        fine_profile_y2 = az_profile(fine_profile_x)

        #angle_knots_y = az_profile(angle_knots)

        plot.plot(ring_angles, ring_data, '+', color="#aaaaaa")
        #plot.plot(fine_profile_x, fine_profile_y)
        plot.plot(fine_profile_x, fine_profile_y2)
        #plot.plot(angle_knots, angle_knots_y, 'x')
        #plot.errorbar(x=binned_angles, y=binned_data, yerr=binned_scatter, xerr=0.5*az_binwidth, marker='o', ls='None')

        plot.xlim([-3.2,-1.4])
        plot.ylim([-0.02,0.02])
        #plot.savefig("profile_az_%d-%d.png" % (ri,ro))

        plot.close()
    print " - done!"

    #------------------------------------------------------------------------------
    #
    # Now we have all components, so compute the pupil ghost template
    #
    #------------------------------------------------------------------------------

    #
    # Go through all the pixels in the data block and compute the pupil ghost.
    # Radial profile is simple, for the azimuthal interpolate linearly between the two radii
    #

    nonradial_profile = numpy.zeros(shape=data.shape)
    for x in range(data.shape[0]):
        sys.stdout.write("\rCreating pupil-ghost template, %.1f %% done ..." % ((x+1)*100.0/data.shape[0]))
        sys.stdout.flush()
        for y in range(data.shape[1]):

            radius_here = radius[x,y]
            if (radius_here < r_inner or radius_here > r_outer):
                continue

            angle_here = angle[x,y]

            r_here = float(radius_here) / dr
            ri = int(math.floor(r_here-0.5))
            ro = int(math.ceil(r_here-0.5))

            val_i, val_o = 0,0
            if (radial_splines[ri] != None):
                val_i = radial_splines[ri](angle_here)
            if (radial_splines[ro] != None):
                val_o = radial_splines[ro](angle_here)

            # Now interpolate linearly between the two
            if (ri == ro):
                nonradial_profile[x,y] = val_i
            else:
                slope = (val_o - val_i) / float(ro - ri)
                nonradial_profile[x,y] = (r_here - float(ri)) * slope + val_i

            if (nonradial_profile[x,y] > pupil_radial_2d[x,y]):
                nonradial_profile[x,y] = pupil_radial_2d[x,y]
    print " complete!"

    fullprofile = nonradial_profile + pupil_radial_2d
    fullprofile[fullprofile<0] = 0

    if (write_intermediate):
        print "Writing data"
        # Save non-radial as fits
        pyfits.PrimaryHDU(data=nonradial_profile).writeto("fit_nonradial.fits", clobber=True)
        pyfits.PrimaryHDU(data=pupil_radial_2d).writeto("fit_radial.fits", clobber=True)

        # print "nonradial=",nonradial_profile.shape
        # print "pupil2d=",pupil_radial_2d.shape
        # nonradial_profile = numpy.min([nonradial_profile, pupil_radial_2d])
        # pyfits.PrimaryHDU(data=nonradial_profile).writeto("fit_nonradial2.fits", clobber=True)

        pyfits.PrimaryHDU(data=fullprofile).writeto("fit_rad+nonrad.fits", clobber=True)

        leftover = bg_sub - fullprofile
        pyfits.PrimaryHDU(data=leftover).writeto("fit_leftover.fits", clobber=True)

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


def bin_azimuthal_profile(angles, data, d_angle, min_angle=None, max_angle=None, verbose=False):

    if (min_angle == None): min_angle = numpy.min(angles)
    if (max_angle == None): max_angle = numpy.max(angles)

    n_angles = int(math.floor((max_angle - min_angle) / d_angle))
    
    median = numpy.zeros(shape=(n_angles))
    angle  = numpy.zeros(shape=(n_angles))
    count  = numpy.zeros(shape=(n_angles))
    var    = numpy.zeros(shape=(n_angles))


    for i in range(n_angles):
        a_start = min_angle + i * d_angle
        a_end   = numpy.min([a_start + d_angle, max_angle])
        
        selected = (angles >= a_start) & (angles < a_end)

        this_sector_data = data[selected]
        this_sector_angle = angles[selected]

        count[i] = numpy.sum(selected)

        # Check if any of the pixels are marked as NaN
        nan_count = numpy.sum(numpy.isnan(this_sector_data))

        good_value = numpy.isfinite(this_sector_data)
        if (numpy.sum(good_value) < 1):
            angle[i] = 0.5*(a_start + a_end)
            median[i] = 0
            var[i] = 1e9
        else:
            angle[i] = numpy.median(this_sector_angle[good_value])
            median[i] = numpy.median(this_sector_data[good_value])
            var[i] = numpy.std(this_sector_data[good_value]) / count[i]

        #if (nan_count > 0.1 * count[i]):
        #    # That's too many NaNs
        #    var[i] = 1e9


    median_var = numpy.median(var)
    var[var == 0] = 1e9
    ## compute mean pixels per bin
    #npixels = numpy.median(count)
    #var[count < 0.2*npixels] = 1e9

    if (verbose):
        print "Results from az binning", a_start, a_end
        print angle[:5]
        print median[:5]
        print count[:5]
        print var[:5]
        print


    return angle, median, var

def sort_angles(angles, data):

    si = numpy.argsort(angles)
    out_angle = numpy.zeros(shape=(angles.shape))
    out_data = numpy.zeros(shape=(data.shape))

    for i in range(si.shape[0]):
        out_angle[i] = angles[si[i]]
        out_data[i]  = data[si[i]]

    return out_angle, out_data



#################################
#
# Important note:
# Many x/y values are swapped, because fits data is arranged in y/x coordinates, not x/y
#
#################################
if __name__ == "__main__":


    # Read in the input parameters
    template = sys.argv[1]
    
    output_filename = sys.argv[2]

    
    r_inner = float(cmdline_arg_set_or_default("-ri",  700))
    r_outer = float(cmdline_arg_set_or_default("-ri", 4000))
    dr = float(cmdline_arg_set_or_default("-dr", 20))
    binfac = int(cmdline_arg_set_or_default("-prebin", 4))

    hdu_ref = pyfits.open(template)

    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4115, -190),
                          "OTA44.SCI": (-90, -230),
                          "OTA43.SCI": (-120, 4150),
                          }

    bpmdir = cmdline_arg_set_or_default("-bpm", None)

    hdulist_out = [hdu_ref[0]]

    for i in range(1, len(hdu_ref)):

        extname = hdu_ref[i].header['EXTNAME']
        if (extname in pupilghost_centers):
            center_x, center_y = pupilghost_centers[extname]

            print "Using center position %d, %d" % (center_y, center_x)

            data = hdu_ref[i].data
            if (bpmdir != None):
                bpmfile = "%s/bpm_xy%s.reg" % (bpmdir, extname[3:5])
                mask_broken_regions(data, bpmfile, True)

            correction = fit_pupilghost(data, (center_y, center_x), (r_inner, r_outer), dr)

            hdu_ref[i].data = correction
            hdulist_out.append(hdu_ref[i])

    hdu_out = pyfits.HDUList(hdulist_out)
    hdu_out.writeto(output_filename, clobber=True)

