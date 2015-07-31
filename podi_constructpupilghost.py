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
import multiprocessing
from astLib import astWCS
import podi_showassociations as podi_associations

from podi_definitions import *
from podi_commandline import *
import podi_imcombine

az_knot_limit = [50,600]

import logging
import podi_logging

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



def fit_spline_background(radii, flux, logger=None):

    if (logger == None):
        logger = logging.getLogger("FitSpline")

    # data = numpy.loadtxt("/home/work/odi_commissioning/pupilghost/radial__x")
    bad = numpy.isnan(radii) | numpy.isnan(flux)
    exclude = (radii > 1150) & (radii < 3950)
    r = radii.copy()[~bad & ~exclude]
    f = flux.copy()[~bad & ~exclude]
    _min, _max = numpy.min(r), numpy.max(r)
    #print "min/max:", _min, _max

    t = numpy.arange(_min+50,_max-50, 100)
    exclude = (t > 1100) & (t < 4000)
    t = numpy.sort(numpy.append(t[~exclude], [1101,1102, 4001, 4002]))
    #print "t:", t

    # Now fit a spline to the data
    #print "R:",r
    spl = scipy.interpolate.LSQUnivariateSpline(x=r, y=f, t=t)

    # generate a smooth curve for plotting & debugging
    xxx = numpy.arange(_min,_max,10)
    numpy.savetxt("radial.spline",
                  numpy.append(xxx.reshape((-1,1)),
                               spl(xxx).reshape((-1,1)), axis=1))


    numpy.savetxt("radial.t",
                  numpy.append(t.reshape((-1,1)),
                               spl(t).reshape((-1,1)), axis=1))

    # now take the original profile, and normalize/background subtract it
    f_spline = spl(radii)
    f_fixed = (flux-f_spline) / f_spline
    numpy.savetxt("radial.fixed",
                  numpy.append(radii.reshape((-1,1)),
                               f_fixed.reshape((-1,1)), axis=1))

    return spl


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



def mp_pupilghost_slice(job_queue, result_queue, bpmdir, binfac):

    _logger = logging.getLogger("MPSlice")
    _logger.debug("Worker started")

    while (True):

        task = job_queue.get()
        if (task == None):
            job_queue.task_done()
            break

        filename, extname = task

        _, bn = os.path.split(filename)
        logger = logging.getLogger("%s(%s)" % (bn[:-5], extname))
        logger.debug("Starting work")

        hdulist = pyfits.open(filename)
        input_hdu = hdulist[extname]
        rotator_angle = hdulist[0].header['ROTSTART'] 

        logger.info("Searching for center ...")
        centering = dev_pgcenter.find_pupilghost_center(input_hdu, verbose=False)
        fx, fy, fr, vx, vy, vr = centering
        center_x = vx
        center_y = vy


        #stdout_write("Using center position %d, %d for OTA %s\n" % (center_y, center_x, extname))
        logger.info("Adding OTA %s, center @ %d, %d" % (extname, center_x, center_y))

        data = input_hdu.data
        if (bpmdir != None):
            bpmfile = "%s/bpm_xy%s.reg" % (bpmdir, extname[3:5])
            logger.debug("Masking bad pixels from %s" % (bpmfile))
            mask_broken_regions(data, bpmfile, verbose=False)

        # Convert into radii and angles to make sure we can subtract the background
        binned, radius, angle = get_radii_angles(data, (center_y, center_x), binfac)
        #print binned.shape, radius.shape, angle.shape, (center_y, center_x)

        # Fit and subtract the background
        bgsub = subtract_background(binned, radius, angle, radius_range, binfac, logger=logger)

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
        #print "insert target: x=", bx, tx, "y=", by, ty
        #combined[bx:tx, by:ty] = binned #bgsub[:,:]
        combined[by:ty, bx:tx] = bgsub[:,:]

        angle_mismatch = compute_angular_misalignment(input_hdu.header)

        # Rotated around center to match the ROTATOR angle from the fits header
        full_angle = rotator_angle + angle_mismatch
        combined_rotated = rotate_around_center(combined, full_angle, mask_nans=True, spline_order=1)
        #combined_rotated = combined

        imghdu = pyfits.ImageHDU(data=combined_rotated)
        imghdu.header['EXTNAME'] = extname
        imghdu.header['ROTANGLE'] = rotator_angle

        imghdu.header['OTA'] = int(extname[3:5])
        imghdu.header['PGCNTRFX'] = fx
        imghdu.header['PGCNTRFY'] = fy
        imghdu.header['PGCNTRFR'] = fr
        imghdu.header['PGCNTRVX'] = vx
        imghdu.header['PGCNTRVY'] = vy
        imghdu.header['PGCNTRVR'] = vr

        result_queue.put((imghdu,extname, centering, angle_mismatch))
        job_queue.task_done()

    _logger.debug("Worker shutting down")
    return


def make_pupilghost_slice(filename, binfac, bpmdir, radius_range, clobber=False):

    hdu_ref = pyfits.open(filename)

    logger = logging.getLogger("MakePGSlice")

    hdus = []
    centers = []

    rotator_angle = hdu_ref[0].header['ROTSTART'] 
    logger.info("Loading frame %s ..." % (filename))

    #combined_file = "pg_combined_%+04d.fits" % numpy.around(rotator_angle)
    #print "combined-file:",combined_file

    output_filename = "pg_%+04d.fits" % (numpy.round(rotator_angle))
    if (os.path.isfile(output_filename) and not clobber):
        logger.warning("output filename %s already exists, skipping\n" % (output_filename))
        return None

    logger.info("creating pupilghost slice %s ..." % (output_filename))
    
    datas = []
    extnames = []
    rotateds = []

    hdulist = [pyfits.PrimaryHDU()]

    job_queue = multiprocessing.JoinableQueue()
    result_queue = multiprocessing.Queue()

    #
    # Start workers
    #
    processes = []
    for i in range(4):
        p = multiprocessing.Process(
            target=mp_pupilghost_slice,
            kwargs={
                'job_queue': job_queue,
                'result_queue': result_queue,
                'bpmdir': bpmdir,
                'binfac': binfac,
                },
            )
        p.start()
        processes.append(p)

    jobs_ordered = 0
    pupilghost_centers = ['OTA33.SCI', 'OTA34.SCI', 'OTA43.SCI', 'OTA44.SCI']
    for i in range(1, len(hdu_ref)):

        extname = hdu_ref[i].header['EXTNAME']
        
        if (extname in pupilghost_centers):
            #print "\n\n\n\n\n",extname

            #
            # Determine center position
            #
            # old method: use fixed values
            # center_x, center_y = pupilghost_centers[extname]

            job_queue.put((filename, extname))
            jobs_ordered += 1

    #
    # Send quit command
    # 
    for p in processes:
        job_queue.put(None)
        
    #
    # Collect results
    #
    centerings = {}
    all_datas = []
    d_angles = {}
    for i in range(jobs_ordered):
        result = result_queue.get()

        imghdu, extname, centering, angle_mismatch = result
        all_datas.append(imghdu.data)
        hdulist.append(imghdu)
        logger.info("Done with OTA %s" % (extname))
        centerings[extname] = centering
        d_angles[extname] = angle_mismatch

    #
    # Create some keywords to report what was done
    #
    norm_angle = numpy.round(rotator_angle if rotator_angle > 0 else rotator_angle + 360.)
    centerf_str = centerv_str = d_angle_str = ""
    ota_str = ""
    for extname in pupilghost_centers:
        fx, fy, fr, vx, vy, vr = centerings[extname]
        ota_str += "%s;" % (extname)
        centerv_str += "%+05d,%+05d;" % (vx, vy)
        centerf_str += "%+05d,%+05d;" % (fx, fy)
        da = d_angles[extname]
        d_angle_str += "%+6.1f;" % (da * 60.)
    
    #
    # Now combine the slices from each of the contributing OTAs
    # 
    logger.info("Combining all OTA slices")
    combined = podi_imcombine.imcombine_data(all_datas, operation='nanmean.bn')

    comb_hdu = pyfits.ImageHDU(data=combined)
    comb_hdu.name = "COMBINED"
    comb_hdu.header['CNTRF%03d' % (norm_angle)] = (centerf_str[:-1], "PG center, fixed r [px]")
    comb_hdu.header['CNTRV%03d' % (norm_angle)] = (centerv_str[:-1], "PG center, var. r [px]")
    comb_hdu.header['ALPHA%03d' % (norm_angle)] = (d_angle_str[:-1], "OTA angle [arcmin]")
    comb_hdu.header['OTAORDER'] = ota_str[:-1]
    comb_hdu.header['ROTANGLE'] = rotator_angle
    comb_hdu.header['RNDANGLE'] = norm_angle
    hdulist.append(comb_hdu)

    hdulist[0].header['RNDANGLE'] = norm_angle

    #
    # Copy the associations table
    #
    try:
        assoctable = hdu_ref['ASSOCIATIONS']
        logger.debug("Transfering the assocations table")
        hdulist.append(assoctable)
    except:
        logger.warning("No associations table found, unable to transfer table")

    logger.info("All done!")
    HDUlist = pyfits.HDUList(hdulist)
    HDUlist.writeto(output_filename, clobber=True)

    return rotator_angle, norm_angle


def subtract_background(data, radius, angle, radius_range, binfac, logger=None):
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

    if (logger == None):
        logger = logging.getLogger("BGSub")

    # Compute the radial bin size in binned pixels
    logger.debug("subtracting background - binfac=%d" % (binfac))
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
    logger.info("Computing background-level ...")
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
    
    #print "XXXXXXX", radii.shape, background_levels.shape
    numpy.savetxt("radial__%s" % ("x"),
                  numpy.append(radii.reshape((-1,1)),
                               background_levels.reshape((-1,1)), axis=1))
    #print "saved"

    # Find average intensity at the largest radii
    avg_level = bottleneck.nanmedian(background_levels[radii>4000])
    #print "avg_level=",avg_level

    #
    # Normalize profile
    #
    normalize_region = ((radii < 1100) & (radii > 600)) |  \
                       ((radii > 4000) & (radii < 4600))
    normalize_flux = numpy.mean(background_levels[normalize_region])
    #print "normalization flux =", normalize_flux

    #
    # Use the profile and fit a spline to the underlying shape
    #
    spl = fit_spline_background(radii, background_levels, logger=logger)

    # #
    # # Subtract background and normalize all measurements
    # #
    # background_levels = (background_levels - normalize_flux) / normalize_flux

    # fitfunc = lambda p, x: p[0] + p[1] * x
    # errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    # bg_for_fit = background_levels
    # #bg_for_fit[numpy.isnan(background_levels)] = 0
    # bg_for_fit[((radii > ri) & (radii < ro))] = 0
    # pinit = [0.0, 0.0] # Assume no slope and constant level of 0
    # out = scipy.optimize.leastsq(errfunc, pinit,
    #                        args=(radii, background_levels, background_level_errors), full_output=1)

    # pfinal = out[0]
    # covar = out[1]
    # stdout_write(" best-fit: %.2e + %.3e * x\n" % (pfinal[0], pfinal[1]))
    # #print pfinal
    # #print covar

    # #
    # # Now we have the fit for the background, compute the 2d background 
    # # image and subtract it out
    # #
    # x = numpy.linspace(0, max_radius, 100)
    # y_fit = radii * pfinal[1] + pfinal[0]
    # background = pfinal[0] + pfinal[1] * radius
    
    # bg_sub = ((data - normalize_flux) / normalize_flux) - background

    # bg_sub_profile = background_levels - (pfinal[0] + pfinal[1]*radii)
    # numpy.savetxt("radial__%s" % ("bgsub"),
    #               numpy.append(radii.reshape((-1,1)),
    #                            bg_sub_profile.reshape((-1,1)), axis=1))

    background_1d = spl(radius.flatten())
    background_2d = background_1d.reshape(radius.shape)

    bg_sub = (data - background_2d) / background_2d
    
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
        #stdout_write(" done!\n")

    # Now we are done with all profiles, write the results to the output file
    clobberfile(outputfile)
    stdout_write("writing output file ...")
    hdulist.writeto(outputfile, clobber=True)
    #stdout_write(" done!\n")



def compute_angular_misalignment(header, l=256):

    #
    # Make copy so we don't accidently change the input header
    #
    ext = pyfits.ImageHDU(header=header)

    #
    # Get valid WCS system
    #
    ext.header['NAXIS'] = 2
    ext.header['NAXIS1'] = l
    ext.header['NAXIS2'] = l
    ext.header['CRVAL1'] = 0.0
    ext.header['CRVAL2'] = 0.0
    wcs = astWCS.WCS(ext.header, mode='pyfits')

    # compute zero-point
    radec_0_0 = numpy.array(wcs.pix2wcs(0,0))
    wcs.header['CRVAL1'] -= radec_0_0[0]
    if (wcs.header['CRVAL1'] < 0): wcs.header['CRVAL1'] += 360.
    wcs.header['CRVAL2'] -= radec_0_0[1]
    wcs.header['CRVAL1'] += 10. # add some offset to RA to avoid problems around RA=0=360
    wcs.updateFromHeader()

    #
    # compute points at origin and along both axes
    #
    radec_0_0 = numpy.array(wcs.pix2wcs(0,0))
    #print radec_0_0-[10.,0]
    radec_0_100 = numpy.array(wcs.pix2wcs(0,l)) - radec_0_0
    radec_100_0 = numpy.array(wcs.pix2wcs(l,0)) - radec_0_0
    radec_100_100 = numpy.array(wcs.pix2wcs(l,l)) - radec_0_0

    #
    # convert vectors into angles
    #
    #print radec_100_0, radec_0_100
    angle_100_0 = numpy.degrees(numpy.arctan2(radec_100_0[0], radec_100_0[1]))
    angle_0_100 = numpy.degrees(numpy.arctan2(radec_0_100[0], radec_0_100[1]))
    angle_100_100 = numpy.degrees(numpy.arctan2(radec_100_100[0], radec_100_100[1]))

    #print angle_100_100 - 45.

    #
    # Then from the difference between perfect alignment compute misalignment
    #
    d = 90. - angle_100_0
    #print "\n",angle_100_0, angle_0_100, angle_100_0 - angle_0_100, d, 0.5*(d-angle_0_100)
    angle_error = 0.5*(d-angle_0_100)

    angle_error = angle_100_100 - 45.

    rot = wcs.getRotationDeg()
    if (rot > 180): rot -= 360.
    #print rot

    #print "Angle_Misalignment (%s) = %f deg" % (ext.name, angle_error)

    return angle_error


def combine_pupilghost_slices(out_filename, filelist, op='sigclipmean'):

    logger = logging.getLogger("CombinePG")

    print filelist
    
    # 
    # Gather information about rotation angles and center positions
    # Also collect the association data.
    #
    rotangles = numpy.ones(len(filelist)) * -9999
    headers = [None] * len(filelist)

    #
    # Prepare the association data
    #
    assoc_table = {}
    logger.info("Reading data from input files")

    for idx, fn in enumerate(filelist):
        hdulist = pyfits.open(fn)
        
        # Read the center and angle keywords from the "COMBINED" extension
        comb_hdu = hdulist['COMBINED']
        _ra = comb_hdu.header['ROTANGLE']
        rotangles[idx] = _ra if _ra > 0 else _ra + 360.
        headers[idx] = comb_hdu.header
        logger.debug("%s: %.2f" % (fn, _ra))
        
        # comb_hdu.header['CNTRF%03d' % (norm_angle)] = (centerf_str[:-1], "PG center, fixed r [px]")
        # comb_hdu.header['CNTRV%03d' % (norm_angle)] = (centerv_str[:-1], "PG center, var. r [px]")
        # comb_hdu.header['ALPHA%03d' % (norm_angle)] = (d_angle_str[:-1], "OTA angle [arcmin]")
        # comb_hdu.header['OTAORDER'] = ota_str[:-1]

        this_assoc = {'pupilghost-slice': fn }
        assoc_table = collect_reduction_files_used(assoc_table, this_assoc)

        # Read the assocation table of this frame
        in_assoc  = podi_associations.read_associations(hdulist)
        if (not in_assoc == None):
            logger.debug("Found assocations:\n%s" % (str(in_assoc)))
            assoc_table = collect_reduction_files_used(assoc_table, in_assoc)

#    return

    #
    # Now sort the rotator angles from smallest to largest
    #
    angle_sort = numpy.argsort(rotangles)

    #
    # Combine all frames
    #
    logger.info("Stacking all slices into master pupilghost template")
    combined_hdulist = podi_imcombine.imcombine(
        filelist, 
        outputfile=None,
        operation=op, #nanmean.bn',
        return_hdu=True,
        subtract=None, scale=None)

    print combined_hdulist
    combined = combined_hdulist['COMBINED']

    combined.header['STACK_OP'] = op

    #
    # Add the sorted keywords back into the resulting file
    #
    logger.info("Adding metadata")
    primhdu = pyfits.PrimaryHDU()
    
    try:
        prev_hdr = None
        for header, label in [
                ('CNTRF%03d', "center, fixed radius"), 
                ('CNTRV%03d', "center, var. radius"),
                ('ALPHA%03d', "angle mismatch [arcmin]"),
        ]:
            first_hdr = None
            for i in range(rotangles.shape[0]):
                idx = angle_sort[i]
                rotangle = rotangles[idx]
                round_angle = headers[idx]['RNDANGLE']
                keyname = header % round_angle
                #logger.info("Adding key: %s" % (keyname))
                combined.header[keyname] = headers[idx][keyname]
                #primhdu.header[keyname] = headers[idx][keyname]
                first_hdr = keyname if first_hdr == None else first_hdr
                if (prev_hdr == None):
                    print "adding header",keyname," somewhere"
                    primhdu.header.append((keyname, headers[idx][keyname]))
                    prev_hdr = keyname
                else:
                    print "adding header",keyname,"after",prev_hdr
                    primhdu.header.insert(prev_hdr, (keyname, headers[idx][keyname]), after=True)
                prev_hdr = keyname

            add_fits_header_title(primhdu.header, label, first_hdr)

        combined.header['OTAORDER'] = headers[0]['OTAORDER']
    except:
        podi_logging.log_exception()
        pass

    print assoc_table
    assoc_hdu = create_association_table(assoc_table)

    logger.info("Writing output to %s" % (out_filename))
    out_hdulist = pyfits.HDUList([primhdu, combined, assoc_hdu])
    out_hdulist.writeto(out_filename, clobber=True)

    logger.info("Work complete!")

#################################
#
# Important note:
# Many x/y values are swapped, because fits data is arranged in y/x coordinates, not x/y
#
#################################
if __name__ == "__main__":

    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)
    
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

    elif (cmdline_arg_isset("-spline")):

        data = numpy.loadtxt(get_clean_cmdline()[1])

        spl = fit_spline_background(data[:,0], data[:,1])


    elif (cmdline_arg_isset("-angles")):

        filename = get_clean_cmdline()[1]
        hdulist = pyfits.open(filename)
        l = 4096

        for ext in hdulist:
            if (not is_image_extension(ext)):
                continue
     
            angle_error = compute_angular_misalignment(ext.header)
            print "Angle_Misalignment (%s) = %f deg" % (ext.name, angle_error)

    elif (cmdline_arg_isset("-combine")):

        out_filename = get_clean_cmdline()[1]
        filelist = get_clean_cmdline()[2:]

        # combine all images
        combine_pupilghost_slices(out_filename, filelist)

    else:
        r_inner = float(cmdline_arg_set_or_default("-ri",  700))
        r_outer = float(cmdline_arg_set_or_default("-ro", 4000))
        dr = float(cmdline_arg_set_or_default("-dr", 20))

        filenames = get_clean_cmdline()[1:]
        print filenames

        radius_range = (r_inner, r_outer, dr)

        for inputfile in filenames:
            make_pupilghost_slice(inputfile, binfac, bpmdir, radius_range, clobber=False)

    podi_logging.shutdown_logging(options)

    sys.exit(0)

