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
    print "Computing radii"
    x, y = numpy.indices(data.shape)
    dx = x - center_x/binfac
    dy = y - center_y/binfac

    radius = numpy.sqrt(dx*dx + dy*dy)
    angle = numpy.arctan2(dx, dy)

    return data, radius, angle


def merge_OTAs(hdus, centers):

    combined = numpy.zeros(shape=(9000,9000))
    combined[:,:] = numpy.NaN

    for i in range(len(hdus)):

        stdout_write("Adding OTA %s ...\n" % (hdus[i].header['EXTNAME']))

        # Use center position to add the new frame into the combined frame
        # bx, by are the pixel position of the bottom left corner of the frame to be inserted
        bx = combined.shape[0] / 2 - centers[i][0]
        by = combined.shape[1] / 2 - centers[i][1]
        print bx, by
        tx, ty = bx + hdus[i].data.shape[0], by + hdus[i].data.shape[1]

        combined[bx:tx, by:ty] = hdus[i].data[:,:]

    return combined

def fit_pupilghost(hdus, centers, rotator_angles, radius_range, dr_full, 
                   write_intermediate=True, show_plots=False):


    # Choose binning of raw-data to speed up computation
    binfac = 4

    # Compute the radial bin size in binned pixels
    dr = dr_full/binfac
    r_inner, r_outer = radius_range
    r_inner /= binfac
    r_outer /= binfac

    # Allocate some memory to hold the template
    template = numpy.zeros(shape=(9000,9000))
    template_binned, template_radius, template_angle = get_radii_angles(template, (4500,4500), 4)

    combined = merge_OTAs(hdus, centers)
    data, radius, angle = get_radii_angles(combined, (combined.shape[0]/2, combined.shape[1]/2), binfac)
    angle -= rotator_angles[0]

    # write some intermediate data products
    if (write_intermediate):
        raw_hdu = pyfits.PrimaryHDU(data=data)
        raw_hdu.writeto("raw.fits", clobber=True)

        combined_hdu = pyfits.PrimaryHDU(data=combined)
        combined_hdu.writeto("raw_combined.fits", clobber=True)

    sys.exit(0)

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


    template_radius_1d = template_radius.ravel()
    template_radial = radial_profile(template_radius.ravel()).reshape(template_radius.shape)
    template_radial[(template_radius > r_outer) | (template_radius < r_inner)] = 0
    pupil_sub_hdu = pyfits.PrimaryHDU(data = template_radial)
    pupil_sub_hdu.writeto("template_radial.fits", clobber=True)

    sys.exit(0)

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

        # Select valid pixels 
        valid_angles = ring_angles[valid_pixels_in_ring]
        valid_data   = ring_data[valid_pixels_in_ring]

        # compute the mean angle difference between two points
        mean_da = (max_angle - min_angle) / valid_angles.shape[0]

        # Select knots for the spline fitting
        number_knots = numpy.sum(valid_pixels_in_ring) / 100
        if (number_knots < 30): number_knots = 30
        angle_knots = numpy.linspace(min_angle, max_angle, number_knots) 

        # Now eliminate knots in regions with no data
        for i in range(angle_knots.shape[0]):
            diff_angle = numpy.fabs(angle_knots[i] - valid_angles)
            min_diffangle = numpy.min(diff_angle)
            if (min_diffangle > 3 * mean_da or angle_knots[i] <= min_angle or angle_knots[i] >= max_angle):
                angle_knots[i] = numpy.NaN
        good_angle_knots = angle_knots[numpy.isfinite(angle_knots)]

        # sort all points in this ring, otherwise LSQUnivariateSplie chokes
        si = numpy.argsort(valid_angles)
        sorted_angles = numpy.zeros(valid_angles.shape)
        sorted_data   = numpy.zeros(valid_data.shape)
        for i in range(si.shape[0]):
            sorted_angles[i] = valid_angles[si[i]]
            sorted_data[i]   = valid_data[si[i]]
        
        # Fit spline
        az_profile = scipy.interpolate.LSQUnivariateSpline(sorted_angles, sorted_data, good_angle_knots, k=3)
        radial_splines[cur_radius] = az_profile

        if (show_plots):
            fine_profile_x = numpy.linspace(min_angle, max_angle, 1000)
            fine_profile_y = az_profile(fine_profile_x)
            plot.plot(ring_angles, ring_data, '+', color="#aaaaaa")
            plot.plot(fine_profile_x, fine_profile_y)
            plot.plot(angle_knots, angle_knots_y, 'x')
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

    # .1
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4115, -190),
                          "OTA44.SCI": (-90, -230),
                          "OTA43.SCI": (-120, 4150),
                          }

    # .2
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4115, -120),
                          "OTA44.SCI": (-20, -230),
                          "OTA43.SCI": (-120, 4150),
                          }

    # .3
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4115, -120),
                          "OTA44.SCI": (-60, -190),
                          "OTA43.SCI": (-150, 4120),
                          }

    # .4
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4095, -140),
                          "OTA44.SCI": (-30, -210),
                          "OTA43.SCI": (-120, 4150),
                          }

    # .5
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4095, -140),
                          "OTA44.SCI": (-50, -180),
                          "OTA43.SCI": (-120, 4170),
                          }

    # .6
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4095, -140),
                          "OTA44.SCI": (-100, -170),
                          "OTA43.SCI": (-120, 4170),
                          }

    # .7
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4095, -140),
                          "OTA44.SCI": (-115, -150),
                          "OTA43.SCI": (-120, 4170),
                          }

    # .8
    pupilghost_centers = {"OTA33.SCI": (4080, 4180),
                          "OTA34.SCI": (4080, -140),
                          "OTA44.SCI": (-115, -150),
                          "OTA43.SCI": (-120, 4170),
                          }

    # .9
    pupilghost_centers = {"OTA33.SCI": (4205, 4225),
                          "OTA34.SCI": (4205,  -95),
                          "OTA44.SCI": (  10, -105),
                          "OTA43.SCI": (   5, 4215),
                          }

    # .10
    pupilghost_centers = {"OTA33.SCI": (4205, 4135),
                          "OTA34.SCI": (4205, -185),
                          "OTA44.SCI": (  10, -195),
                          "OTA43.SCI": (   5, 4125),
                          }

    # .11
    pupilghost_centers = {"OTA33.SCI": (4190, 4150),
                          "OTA34.SCI": (4205, -185),
                          "OTA44.SCI": (  10, -195),
                          "OTA43.SCI": (   5, 4125),
                          }

    # .12
    pupilghost_centers = {"OTA33.SCI": (4190, 4150),
                          "OTA34.SCI": (4205, -165),
                          "OTA44.SCI": (  10, -185),
                          "OTA43.SCI": (  10, 4120),
                          }

    # .13
    pupilghost_centers = {"OTA33.SCI": (4190, 4150),
                          "OTA34.SCI": (4205, -165),
                          "OTA44.SCI": (  10, -185),
                          "OTA43.SCI": (  10, 4130),
                          }

    bpmdir = cmdline_arg_set_or_default("-bpm", None)

    hdulist_out = [hdu_ref[0]]

    hdus = []
    centers = []

    for i in range(1, len(hdu_ref)):

        extname = hdu_ref[i].header['EXTNAME']
        rotator_angles = [hdu_ref[0].header['ROTOFF']]

        if (extname in pupilghost_centers):
            center_x, center_y = pupilghost_centers[extname]

            print "Using center position %d, %d for OTA %s" % (center_y, center_x, extname)

            data = hdu_ref[i].data
            if (bpmdir != None):
                bpmfile = "%s/bpm_xy%s.reg" % (bpmdir, extname[3:5])
                mask_broken_regions(data, bpmfile, True)

            hdus.append(hdu_ref[i])
            centers.append((center_y, center_x))

    correction = fit_pupilghost(hdus, centers, rotator_angles, (r_inner, r_outer), dr)

            #hdu_ref[i].data = correction
            #hdulist_out.append(hdu_ref[i])

    #hdu_out = pyfits.HDUList(hdulist_out)
    #hdu_out.writeto(output_filename, clobber=True)

