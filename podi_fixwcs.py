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

Collection of routines used to perform the astrometric correction.

Warning:
--------

    Most functionality is now superseded by ccmatch.

"""

from __future__ import print_function

import sys
import numpy
import os
from query_usno import query_usno
from podi_definitions import *
from podi_commandline import *
import pyfits
#import date
import datetime
#import pywcs  
from astLib import astWCS
import pdb
import scipy
import scipy.stats


import podi_search_ipprefcat
import podi_findstars
import matplotlib.pyplot as pl

arcsec = 1./3600.
number_bright_stars = 100
max_offset = 0.1


N = 200
N_brightest_ota = 70
N_brightest_ref = 150

# Compute matches in smaller blocks to not blow up memory
# Improve: Change execution to parallel !!!
blocksize = 100


use_usno = False
use_sextractor = False
IPP_DIR = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

def compute_wcs_quality(odi_2mass_matched, ota, hdr=None):

    valid_match = odi_2mass_matched[:,2] >= 0
    odi_2mass_matched = odi_2mass_matched[valid_match]

    d_dec = (odi_2mass_matched[:,1] - odi_2mass_matched[:,3])
    d_ra  = ((odi_2mass_matched[:,0] - odi_2mass_matched[:,2]) 
             * numpy.cos(numpy.radians(odi_2mass_matched[:,1])))
    if (ota is not None):
        in_this_ota = odi_2mass_matched[:,10] == ota
        d_ra = d_ra[in_this_ota]
        d_dec = d_dec[in_this_ota]

    numpy.savetxt("wcsquality.test", odi_2mass_matched)
    d_total = numpy.hypot(d_ra, d_dec)
    wcs_scatter = numpy.median(d_total)
    wcs_scatter2 = numpy.std(d_total)
    wcs_mean_dra = numpy.median(d_ra) * 3600.
    wcs_mean_ddec = numpy.median(d_dec) * 3600.
    rms_dra = numpy.sqrt(numpy.mean(d_ra**2)) * 3600.
    rms_ddec = numpy.sqrt(numpy.mean(d_dec**2)) * 3600.
    print("WCS quality:", ota, wcs_mean_dra*3600., wcs_mean_ddec*3600., wcs_scatter*3600., wcs_scatter2*3600., rms_dra, rms_ddec)

    def make_valid(x):
        return x if numpy.isfinite(x) else -9999

    results = {}
    results['RMS-RA'] = rms_dra
    results['RMS-DEC'] = rms_ddec
    results['RMS'] = numpy.hypot(rms_dra, rms_ddec)
    results['MEDIAN-RA'] = wcs_mean_dra
    results['MEDIAN-DEC'] = wcs_mean_ddec

    if (hdr is not None):
        hdr["WCS_RMSA"] = (make_valid(results['RMS-RA']),     "RA r.m.s. of WCS matching [arcsec]")
        hdr["WCS_RMSD"] = (make_valid(results['RMS-DEC']),    "DEC r.m.s. of WCS matching [arcsec]")
        hdr["WCS_RMS"] =  (make_valid(results['RMS']),        "r.m.s. of WCS matching [arcsec]")
        hdr["WCS_ERRA"] = (make_valid(results['MEDIAN-RA']),  "RA median error WCS matching [arcsec]")
        hdr["WCS_ERRD"] = (make_valid(results['MEDIAN-DEC']), "DEC median error of WCS matching [arcsec]")
 
    return results



def count_matches(ota_x, ota_y, ref_x, ref_y, verbose=False, max_offset=0.1, matching_radius_arcsec=3.6):

    ota_count = ota_x.shape[0]
    ref_count = ref_x.shape[0]
    # print "Using %d and %d stars in shift_align_wcs" % (ota_count, ref_count)

    matching_radius = matching_radius_arcsec * arcsec
    max_d2 = matching_radius * matching_radius

    if (verbose): print("Matching radius",matching_radius, max_d2)
    
    index_ota = numpy.arange(ota_count).reshape((ota_count,1)).repeat(ref_count, axis=1)
    index_ref = numpy.arange(ref_count).reshape((1,ref_count)).repeat(ota_count, axis=0)
    #print index_ota[0:10,0:10]
    #print index_ref[0:10,0:10]
    #print index_ota.shape, index_ref.shape

    #print index_ota[280,150], index_ref[280,150]

    match_results = numpy.zeros(shape=(ota_count, ref_count, 3))
    match_results[:,:,:] = numpy.NaN

    #for o in range(ota_count):
    #    print "PO",ota_x[o], ota_y[o]

    #for r in range(ref_count):
    #    print "PR",ref_x[r], ref_y[r]

    matched_pair = 0
    
    median_delta = numpy.median(ref_y)
    cos_delta = math.cos(math.radians(median_delta))

    for o in range(ota_count):
        for r in range(ref_count):
            
            # Take some random shift
            shift_dx = ref_x[r] - ota_x[o]
            shift_dy = ref_y[r] - ota_y[o]
            if (math.fabs(shift_dx*cos_delta) > max_offset or math.fabs(shift_dy) > max_offset):
                continue

            # Apply shift to OTA coordinates
            aligned_x = ota_x + shift_dx
            aligned_y = ota_y + shift_dy

            # Compute distance between every star in the OTA 
            # to every star in the reference frame

            dx1 = aligned_x.reshape((ota_count, 1)).repeat(ref_count, axis=1)
            dx2 = ref_x.reshape((1, ref_count)).repeat(ota_count, axis=0)
            dx = dx1 - dx2
            #print dx.shape

            dy1 = aligned_y.reshape((ota_count, 1)).repeat(ref_count, axis=1)
            dy2 = ref_y.reshape((1, ref_count)).repeat(ota_count, axis=0)
            dy = dy1 - dy2
            #print dy.shape

            d2 = dx*dx*cos_delta**2 + dy*dy
            #print d.shape

            close_match = d2 < max_d2
            matched = numpy.sum(close_match)

            match_results[o,r,0] = matched
            match_results[o,r,1] = shift_dx
            match_results[o,r,2] = shift_dy

            #print "X_DXDY_count",o, r, shift_dx, shift_dy, matched
            if (matched > 1 and False):
                idx_ota, idx_ref = index_ota[close_match], index_ref[close_match]
                #print idx_ota.shape, idx_ref.shape
                for match in range(idx_ota.shape[0]):
                    print("MATCH_%d"%(matched_pair),ota_x[idx_ota[match]], ota_y[idx_ota[match]], \
                        ref_x[idx_ref[match]], ref_y[idx_ref[match]], shift_dx, shift_dy)

                matched_pair += 1
                pass

    #
    # Now that we have the counts for all pairs, find out some 
    # median background level of the number of (random) matches
    #
    good_datapoints = numpy.isfinite(match_results[:,:,0])

    selected_matches = match_results[good_datapoints]
    if (verbose): print("Selected matches=",selected_matches.shape)

    return selected_matches


def shift_align_wcs(ota_x, ota_y, ref_x, ref_y, verbose=False, max_offset=0.1):

    selected_matches = count_matches(ota_x, ota_y, ref_x, ref_y, verbose, max_offset)

    if (selected_matches.shape[0] <= 1):
        # No matches found, return no shift and signal that the solution is invalid (Nmatch<0) 
        return 0, 0, -1, None
        
    median = numpy.median(selected_matches[:,0])
    std = numpy.std(selected_matches[:,0])
    sigma = scipy.stats.scoreatpercentile(selected_matches[:,0], 84) - median

    if (verbose): print("MATCH-RESULTS:")
    if (verbose): print("Random matches:",median,"+/-",std,"   (",sigma,")")
    most_matches_idx = numpy.argmax(selected_matches[:,0])
    if (verbose): print("max position", most_matches_idx)
    most_matches = selected_matches[most_matches_idx,0]
    if (verbose): print("maximum matches:", most_matches)

    #max_index = numpy.unravel_index(max_position, (match_results.shape[0:2]))
    #print max_index
    #print match_results[max_index[0],max_index[1],:]
    if (verbose): print(selected_matches[most_matches_idx])

    #
    # Now to be sure, make sure the peak we found is significant
    #
    if (most_matches >= median + 3*sigma):
        if (verbose): print("Found seemingly significant match:")
    else:
        if (verbose): print("Found match, but not sure how reliable it is:")
    if (verbose): print(selected_matches[most_matches_idx])
    if (verbose): print("Random matches: %d +/- %d (1sigma) (3-sigma: %d)" % (median, sigma, median+3*sigma))

    return selected_matches[most_matches_idx,1], selected_matches[most_matches_idx,2], \
        selected_matches[most_matches_idx,0], \
        selected_matches 


def refine_wcs_shift(ref_x, ref_y, ota_x, ota_y, best_guess, 
                     ota_px_x=None, ota_px_y=None, 
                     alignment_checkfile=None, verbose=False, matching_radius=3):

    #print "INPUT TO REFICE_WCS_SHIFT: =", ref_x.shape, ref_y.shape, ota_x.shape, ota_y.shape, best_guess.shape

    matches = numpy.zeros(shape=(ota_x.shape[0],2))
    # 2 columns: dx, dy
    matched = numpy.zeros(shape=(ota_x.shape[0],8))
    # 8 columns: ota_x, ota_y, ref_x, ref_y, id, matched_id, ota_pixel_x, ota_pixel_y

    if (ota_px_x is None):
        ota_px_x = numpy.ones(shape=(ota_x.shape[0])) * -999
    if (ota_px_y is None):
        ota_px_y = numpy.ones(shape=(ota_y.shape[0])) * -999

        #for i in range(ref_x.shape[0]):
        #    print "REF",ref_x[i], ref_y[i]
        #for i in range(ota_cat.shape[0]):
        #    print "SEX",ota_cat[i,6], ota_cat[i,7]

    #
    # Go through the ODI catalog and search for the closest match from the reference catalog
    # Do this is blocks as a compromise between speed and memory usage
    #

    full_aligned_x = ota_x + best_guess[0]
    full_aligned_y = ota_y + best_guess[1]

    #pl.plot(ref_x, ref_y, "ro", full_aligned_x, full_aligned_y, "b.")
    #pl.title("refine_wcs_shift:   Reference/OTA after shift")
    #pl.show(block=True)

    start = 0
    while(start < ota_x.shape[0]):

        # select a sub-sample o stars to keep memory under control
        aligned_x = full_aligned_x[start:start+blocksize]
        aligned_y = full_aligned_y[start:start+blocksize]

        #
        # Compute distances between each star in the OTA subsample and the full reference catalog
        #

        # x (=RA) needs to be multiplied with cos(declination) 
        # to prevent problems at large declinations
        dx1 = aligned_x.reshape((aligned_x.shape[0], 1)).repeat(ref_x.shape[0], axis=1)
        dx2 = ref_x.reshape((1, ref_x.shape[0])).repeat(aligned_x.shape[0], axis=0)
        dx = (dx2 - dx1) * numpy.cos(numpy.radians(dx1))

        # No special treatment for declinations (=y) necessary
        dy1 = aligned_y.reshape((aligned_y.shape[0], 1)).repeat(ref_y.shape[0], axis=1)
        dy2 = ref_y.reshape((1, ref_y.shape[0])).repeat(aligned_y.shape[0], axis=0)
        dy = dy2 - dy1

        # Compute distances between them
        d2 = dx*dx + dy*dy

        # And sort them by distance
        closest = numpy.argmin(d2, axis=1)

        # Compute the vector from each OTA star to the closest star in the reference catalog
        x = open("dump2", "a")
        for i in range(closest.shape[0]):
            matches[start+i,0] = dx[i, closest[i]]
            matches[start+i,1] = dy[i, closest[i]]
            #print matches[start+i,:]

            # Save the coordinate pairs and which reference star is the closest match
#            matched[start+i,0:6] = [aligned_x[i], aligned_y[i], ref_x[closest[i]], ref_y[closest[i]], start+i, closest[i]]
            matched[start+i,:] = [aligned_x[i], aligned_y[i], ref_x[closest[i]],
                                  ref_y[closest[i]], start+i, closest[i], ota_px_x[start+i], ota_px_y[start+i] ]
            print(aligned_x[i], aligned_y[i], ref_x[closest[i]], ref_y[closest[i]],
                  start+i, closest[i], ota_px_x[start+i], ota_px_y[start+i], file=x)
        x.close()

        start += blocksize
        # Continue with next block

    #
    # Now with all the matches for this OTA, left's determine the best solution
    # Again, use the 3-sigma clipping as above
    #
    # Start with limits excluding stars with separations > matching_radius arcsec

    #return 0, matched, matches

    dmin, dmax = [-matching_radius*arcsec, -matching_radius*arcsec], [matching_radius*arcsec, matching_radius*arcsec]
    clipping = matches.copy()
    for rep in range(2):
        valid_range = (clipping[:,0] < dmax[0]) & (clipping[:,0] > dmin[0]) & (clipping[:,1] > dmin[1]) & (clipping[:,1] < dmax[1])
        if (verbose): print("valid_range.shape=",valid_range.shape)
        if (verbose): print("count=",numpy.sum(valid_range))

        if (numpy.sum(valid_range) <= 0):
            median = [0,0]
            break

        median = numpy.median(clipping[valid_range], axis=0)
        sigma_p = scipy.stats.scoreatpercentile(clipping[valid_range], 84)
        sigma_m = scipy.stats.scoreatpercentile(clipping[valid_range], 16)
        sigma = 0.5*(sigma_p - sigma_m)
        dmin = median - 3*sigma
        dmax = median + 3*sigma
        if (verbose): print("MEDIAN/SIGMA=",median, sigma, dmin, dmax)

    #
    # Dump some info into a file for later analysis and potentially plotting
    #
    if (alignment_checkfile is not None):
        for i in range(matched.shape[0]):
            print (matched[i,0], matched[i,1], matched[i,2], matched[i,3], ext,
                matched[i,0]-median[0],matched[i,1]-median[1],matches[i,0], matches[i,1], file=alignment_checkfile)
        print("\n\n\n\n\n", file=alignment_checkfile)

    if (verbose): print("Ref: %d, ODI: %d, Matched: %d" % (ref_x.shape[0], ota_x.shape[0], matched.shape[0]))
    return median, matched, matches



def pick_brightest(ra, dec, mag, N):

    # print "PICK_BRIGHTEST:", ra.shape, N
    
    if (N >= ra.shape[0]):
        return ra, dec, mag
    
    magsort = numpy.argsort(mag)

    _ra = numpy.zeros(shape=(N,))
    _dec = numpy.zeros(shape=(N,))
    _mag = numpy.zeros(shape=(N,))

    for i in range(N):
        _ra[i] = ra[magsort[i]]
        _dec[i] = dec[magsort[i]]
        _mag[i] = mag[magsort[i]]

    return _ra, _dec, _mag





def get_overall_best_guess(frame_shift):
    
    full_shift = numpy.median(frame_shift, axis=0)
    full_shift_std  = numpy.std(frame_shift, axis=0)
    #print full_shift.shape

    #print frame_shift
    #print full_shift,full_shift_std

    dmin, dmax = [-1, -1], [1,1]
    clipping = frame_shift.copy()
    #print full_shift,full_shift_std,dmin,dmax
    for rep in range(2):
        valid_range \
            = (clipping[:,0] < dmax[0]) & (clipping[:,0] > dmin[0]) \
            & (clipping[:,1] > dmin[1]) & (clipping[:,1] < dmax[1])
               
        median = numpy.median(clipping[valid_range], axis=0)
        sigma_p = scipy.stats.scoreatpercentile(clipping[valid_range], 84)
        sigma_m = scipy.stats.scoreatpercentile(clipping[valid_range], 16)
        sigma = 0.5*(sigma_p - sigma_m)
        dmin = median - 3*sigma
        dmax = median + 3*sigma
        print(median, sigma, dmin, dmax)

    best_guess = median
    print("Best shift solution so far:",best_guess)
    return best_guess

def fixwcs(fitsfile, output_filename, starfinder="findstars", refcatalog="ippref"):

    tmp, dummy = os.path.split(sys.argv[0])
    dot_config_dir = tmp + "/config/"
    print(dot_config_dir)
    
    hdulist = pyfits.open(fitsfile)
    ra, dec = hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']

    if (refcatalog == "usno"):
        #
        # Obtain catalog listing from USNO-B1
        #
        usno_cat = fitsfile[:-5]+".usno.cat"
        download_cat = (not cmdline_arg_isset("-skip_ref")) or (not os.path.isfile(usno_cat))
        if (not download_cat):
            stdout_write("Using local copy of USNO reference catalog for WCS calibration\n")
        usno = query_usno(ra, dec, 45, 100000, usno_cat, download=download_cat)
        stdout_write("Read %s stars from USNO-B1\n" % (usno.shape[0]))

        usno_dumpfile = fitsfile[:-5]+".usno.coord"
        print("USNO dump-file:", usno_dumpfile)
        usno_dump = open(usno_dumpfile, "w")
        for i in range(usno.shape[0]):
            print(usno[i, 0], usno[i, 1], file=usno_dump)
        usno_dump.close()

        ref_ra = usno[:,0]
        ref_dec = usno[:,1]
        ref_mag = usno[:,4]
    else:
        ipp_cat = podi_search_ipprefcat.get_reference_catalog(ra, dec, 0.7, IPP_DIR)
        ref_ra = ipp_cat[:,0]
        ref_dec = ipp_cat[:,1]
        ref_mag = ipp_cat[:,3]

    if (starfinder == "sextractor"):
        #
        # Run Sextractor on the input frame
        #
        sex_catalogfile = fitsfile[:-5]+".source.cat"
        sex_config_file = dot_config_dir+"fixwcs.conf"
        sex_param_file = dot_config_dir+"fixwcs.param"
        sex_logfile = fitsfile[:-5]+".sextractor.log"
        sex_cmd = "sex -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s >& %s" % (sex_config_file, sex_param_file, sex_catalogfile, fitsfile, sex_logfile)
        print(sex_cmd)
        stdout_write("Running SExtractor to search for stars, be patient (logfile: %s) ..." % (sex_logfile))
        if ((not cmdline_arg_isset("-skip_sex")) or (not os.path.isfile(sex_catalogfile))):
            os.system(sex_cmd)
        else:
            stdout_write("Re-using previous source catalog\n")
        stdout_write(" done!\n")

        # Read the Sextractor output catalog
        sex_cat = numpy.loadtxt(sex_catalogfile)
        stdout_write("Reading %d stars from Sextractor catalog\n" % sex_cat.shape[0])

        # Eliminate all stars with flags
        sex_cat = sex_cat[sex_cat[:,10] == 0]
        stdout_write("%d stars left after eliminating flags\n" % sex_cat.shape[0])

        sex_dumpfile = fitsfile[:-5]+".sex.coord"
        sex_dump = open(sex_dumpfile, "w")
        for i in range(sex_cat.shape[0]):
            print(sex_cat[i, 6], sex_cat[i, 7], file=sex_dump)
        sex_dump.close()

        odi_ra = sex_cat[:,6]
        odi_dec = sex_cat[:,7]
        odi_mag = sex_cat[:,2]
        odi_ota = sex_cat[:,1]
    else:
        full_source_cat = None
        for ext in range(1, len(hdulist)):
            source_cat = podi_findstars.find_stars(hdulist[ext], binning=4, boxsize=24, dumpfile=None, verbose=False,
                                    detect_threshold=1.5, detect_minarea=6, roundness_limit=[-0.2,+0.2])
            source_cat[:,7] = ext
            if (full_source_cat is None):
                full_source_cat = source_cat
            else:
                full_source_cat = numpy.append(full_source_cat, source_cat, axis=0)
            #print source_cat.shape, full_source_cat.shape
            
        # Now all OTAs have been searched        

        odi_ra = full_source_cat[:,0]
        odi_dec = full_source_cat[:,1]
        odi_mag = -2.5 * numpy.log10(full_source_cat[:,6]) + 30
        odi_ota = full_source_cat[:,7]


    #
    # Now go through each of the extension
    # Improve: Change execution to parallel !!!
    #
    frame_shift = numpy.ones(shape=(len(hdulist)-1, 2)) * -999
    global_number_matches = None
    for ext in range(1, len(hdulist)):
        #
        # Select a sub-catalog with stars in this OTA
        #

        this_ota = odi_ota == ext
        ota_odi_ra = odi_ra[this_ota]
        ota_odi_dec = odi_dec[this_ota]
        ota_odi_mag = odi_mag[this_ota]
        
        
        # Select the N brightest stars to make alignment faster 
        # in the case of many stars, i.e. close to the milky way plane 
        # or in the case of resolved external galaxies
        if (ota_odi_ra.shape[0] > N_brightest_ota):
            ota_odi_ra, ota_odi_dec, ota_odi_mag = pick_brightest(ota_odi_ra, ota_odi_dec, ota_odi_mag, N_brightest_ota)
        elif (ota_odi_ra.shape[0] < 5):
            stdout_write("Insufficient number of ODI stars (%d) for alignment!\n" % (ota_odi_ra.shape[0]))
            continue
        
        #
        # Compute the median central position of this OTA
        #
        #center_ra, center_dec
        centerx, centery = hdulist[ext].header['NAXIS1']/2, hdulist[ext].header['NAXIS2']/2
        ota_center_ra  = (centerx-hdulist[ext].header['CRPIX1'])*hdulist[ext].header['CD1_1'] \
            + (centery-hdulist[ext].header['CRPIX2'])*hdulist[ext].header['CD1_2'] \
            + hdulist[ext].header['CRVAL1']
        ota_center_dec = (centerx-hdulist[ext].header['CRPIX1'])*hdulist[ext].header['CD2_1'] \
            + (centery-hdulist[ext].header['CRPIX2'])*hdulist[ext].header['CD2_2'] \
            + hdulist[ext].header['CRVAL2']
        print("center:", ra, dec)


        #
        # Now select stars that are within a given distance of this OTA
        # Assume a (half-)width of 6 arcmin, that's quite a bit larger than the 4armin of one OTA
        #
        width_ra = 0.1/math.cos(numpy.radians(dec))
        width_dec = 0.1
        print("RA-Range:",ra-width_ra,ra+width_ra)
        print("DEC-Range:",dec-width_dec,dec+width_dec)

        min_ra = ota_center_ra - width_ra
        max_ra = ota_center_ra + width_ra
        min_dec = ota_center_dec - width_dec
        max_dec = ota_center_dec + width_dec
        # Add here: Treatment in the case of ra <~ 0

        ref_this_ota = (ref_ra > min_ra) & (ref_ra < max_ra) & (ref_dec > min_dec) & (ref_dec < max_dec)
        ota_ref_ra = ref_ra[ref_this_ota]
        ota_ref_dec = ref_dec[ref_this_ota]
        ota_ref_mag = ref_mag[ref_this_ota]

        #
        # Make sure we have not too few and not too many reference stars to work with
        #
        if (ota_ref_ra.shape[0] > N_brightest_ref):
            ota_ref_ra, ota_ref_dec, ota_ref_mag = pick_brightest(ota_ref_ra, ota_ref_dec, ota_ref_mag, N_brightest_ref)
        if (ota_ref_ra.shape[0] < 5 or ota_odi_ra.shape[0] < 5):
            stdout_write("Insufficient number of stars (ODI: %d, REF: %d) for alignment!\n" % (ota_odi_ra.shape[0], ota_ref_ra.shape[0]))
            continue

        #
        # Now find the number of matches between the ODI and REF catalogs
        #
        dx, dy, n, number_matches = shift_align_wcs(ota_odi_ra, ota_odi_dec, ota_ref_ra, ota_ref_dec)

        frame_shift[ext-1,:] = [dx, dy]
        
        tmp = numpy.ones(shape=(number_matches.shape[0], number_matches.shape[1]+1))
        tmp[:,:-1] = number_matches
        tmp[:,-1] *= ota_id
        if (global_number_matches is None):
            global_number_matches = tmp
        else:
            global_number_matches = numpy.append(global_number_matches, tmp, axis=0)
                
        print("#######################")
        print("##")
        print("## OTA %d" % (ext))
        print("## dx,dy = %.5f, %.5f" % (dx, dy))
        print("## #matches = %d" % (n))
        print("##")
        print("#######################")

    #
    # Now that the individual OTAs have returned their shift positions, 
    # compute the mean shift for the full frame
    # Also do some 2-step iterative sigma clipping to get rid of outliers
    #
    #best_guess = get_overall_best_guess(frame_shift)
    best_guess_x, contrast, drxy = podi_fixwcs.optimize_shift(global_number_matches)
    best_guess = best_guess_x[0:2]

    #
    # Now that we have a pretty good guess on what the offsets are, 
    # extend the search on the full catalog and determine a better 
    # and more accurate solution!
    #


    #
    # In an ideal world, the remainder should yield the same results for all OTAs, 
    # but for now keep them separate and check if that's the case
    #

    #ref_x = usno[:,0]
    #ref_y = usno[:,1]
    #print ref_x.shape, ref_y.shape

    alignment_checkfile = cmdline_arg_set_or_default("-checkalign", None)
    #alignment_check = None
    print("ali-check=",alignment_checkfile)
    if (alignment_checkfile != None):
        alignment_check = open(alignment_checkfile, "w")
        print("Opened alignment check file")

    if (cmdline_arg_isset("-singleota_shift")):
        for ext in range(1, len(hdulist)):

            # Extract just the source catalog for this OTA
            ota_ra = odi_ra[odi_ota ==ext]
            ota_dec = odi_dec[odi_ota ==ext]

            # Compute the average shift
            median = refine_wcs_shift(ref_ra, ref_dec, ota_ra, ota_dec, best_guess, alignment_check)

            # Add the previous (best-guess) shift and the new refinement
            full_shift = best_guess + median

            # And apply this offset to the CRVAL headers
            hdulist[ext].header['CRVAL1'] += full_shift[0]
            hdulist[ext].header['CRVAL2'] += full_shift[1]

            # Add some HISTORY comments to keep track of what's happening
            hdulist[ext].header.add_history("GLOBAL SHIFT: dRA,dDEC=")
            hdulist[ext].header.add_history("dRA=%.6f dDEC=%.6f" % (best_guess[0], best_guess[1]))
            hdulist[ext].header.add_history("LOCAL SHIFT: dRA,dDEC=")
            hdulist[ext].header.add_history("dRA=%.6f dDEC=%.6f" % (median[0], median[1]))
    else:

        # Compute the refinement for the full source catalog across all OTAs
        median = refine_wcs_shift(ref_ra, ref_dec, odi_ra, odi_dec, best_guess, alignment_check)

        # Add the previous (best-guess) shift and the new refinement
        full_shift = best_guess + median

        # Apply refinement to all OTAs
        for ext in range(1, len(hdulist)):
            apply_wcs_shift(full_shift, hdulist[ext].header)
            
            
    # All work done, close the check-file
    if (alignment_checkfile != None):
        alignment_check.close()

    if (not cmdline_arg_isset('-computeonly')):
        stdout_write("Writing output file...\n")
        hdulist.writeto(output_filename, clobber=True)
    else:
        stdout_write("Skipping the output file, you requested -computeonly\n")



def apply_wcs_shift(shift, hdr, fixrot_trans=None, verbose=False):

    hdr['CRVAL1'] += shift[0]
    hdr['CRVAL2'] += shift[1]

    hdr['WCSOF_RA'] = (shift[0], "WCS Zeropoint offset RA")
    hdr['WCSOFDEC'] = (shift[1], "WCS Zeropoint offset DEC")

    if (fixrot_trans != None):

        angle = fixrot_trans[2]
        
        # Compute the new CRVAL1/2 values
        cr1 = hdr['CRVAL1']
        cr2 = hdr['CRVAL2']
        hdr['CRVAL1'] =  math.cos(angle)*cr1 + math.sin(angle)*cr2 + fixrot_trans[0]
        hdr['CRVAL2'] =  math.cos(angle)*cr2 - math.sin(angle)*cr1 + fixrot_trans[1] 
        if (verbose):
            print("CR1/2 before=",cr1, cr2)
            print("CR1/2 after =",hdr['CRVAL1'],  hdr['CRVAL2'])

        cd11 = hdr['CD1_1']
        cd12 = hdr['CD1_2']
        cd21 = hdr['CD2_1']
        cd22 = hdr['CD2_2']

        co = math.cos(angle)
        si = math.sin(angle)
        hdr['CD1_1'] = cd11*co - cd12*si
        hdr['CD1_2'] = cd11*si + cd12*co
        hdr['CD2_1'] = cd21*co - cd22*si
        hdr['CD2_2'] = cd21*si + cd22*co
        
        hdr['WCSMROT1'] = fixrot_trans[0]
        hdr['WCSMROT2'] = fixrot_trans[1]
        hdr['WCSMROT3'] = math.degrees(fixrot_trans[2])*60.

    # Make sure the RA range stays in valid ranges
    if (hdr['CRVAL1'] < 0): hdr['CRVAL1'] += 360 
    if (hdr['CRVAL1'] >360): hdr['CRVAL1'] -= 360 

    # Also add some history 
    #hdr.add_history("LOCAL SHIFT: dRA,dDEC=")
    #hdr.add_history("dRA=%.6f dDEC=%.6f" % (median[0], median[1]))

    #hdr.add_history("GLOBAL SHIFT: dRA,dDEC=")
    #hdr.add_history("dRA=%.6f dDEC=%.6f" % (best_guess[0], best_guess[1]))

    return




def apply_fit_params(wcsout, p):

    theta = numpy.radians(p[0])
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    cd11, cd12, cd21, cd22 = wcsout['CD1_1'], wcsout['CD1_2'], wcsout['CD2_1'], wcsout['CD2_2']
    wcsout['CD1_1'] = cd11 * cos_theta + cd12 * sin_theta
    wcsout['CD1_2'] = cd12 * cos_theta - cd11 * sin_theta
    wcsout['CD2_1'] = cd21 * cos_theta + cd22 * sin_theta
    wcsout['CD2_2'] = cd22 * cos_theta - cd21 * sin_theta

    wcsout['CRVAL1'] += p[1]
    wcsout['CRVAL2'] += p[2]

    return

def rotate_optimize(ota_list, extensions, ref_ra, ref_dec, odi_x, odi_y):

    wcs = [None] * len(ota_list)
    for ext in range(len(ota_list)):
        wcs[ext] = astWCS.WCS(ota_list[ext].header, mode="pyfits")

    
    def odi2wcs(odi_x, odi_y, extensions, wcs):
        ra, dec = numpy.zeros(odi_x.shape[0]), numpy.zeros(odi_y.shape[0])
        for src in range(odi_x.shape[0]):
            ext = int(extensions[src])
            ra[src], dec[src] = wcs[ext].pix2wcs(odi_x[src], odi_y[src])
        return ra, dec


    def d_radec(p, odi_x, odi_y, ref_ra, ref_dec, extensions, wcs):
        # Update the WCS by applying the shift and rotation
        wcs_thistime = [None] * len(wcs)
        for i in range(1, len(wcs)):
            if (wcs[i] != None):
                wcs_thistime[i] = wcs[i].copy()
                #apply_fit_params(wcs[i], p)
                #try:
                #    print i, wcs_thistime[i].header['CRPIX1']
                #except:
                #    print i, "no header"
                #continue
                apply_fit_params(wcs_thistime[i].header, p)
                #apply_fit_params(wcs_thistime[i].header, wcs[i].header, p)
                wcs_thistime[i].updateFromHeader()

        print(p)

        # Compute the updated RA/DEC from the pixel X/Ys
        ra, dec = odi2wcs(odi_x, odi_y, extensions, wcs_thistime)
        d_ra = ref_ra - ra
        d_dec = ref_dec - dec
        d_radec = d_ra**2 + d_dec**2

        # Mask out drastic outliers
        return numpy.sqrt(d_radec)


    df = open("dump", "w")
    ra, dec = odi2wcs(odi_x, odi_y, extensions, wcs)

    pinit = [0, 0, 0]
    out = scipy.optimize.leastsq(d_radec, pinit, 
                                 args=(odi_x, odi_y, ref_ra, ref_dec, extensions, wcs), 
                                 full_output=1)

    p_bestfit = out[0]
    
    #xxx = d_radec(pinit, odi_x, odi_y, ref_ra, ref_dec, extensions, wcs)
    #print xxx.shape, "\n",xxx

    for src in range(odi_x.shape[0]):
        print(ref_ra[src], ref_dec[src], odi_x[src], odi_y[src], ra[src], dec[src], file=df)

    df.close()


    # Compute Ra/Dec with original WCS
    ra_orig, dec_orig = odi2wcs(odi_x, odi_y, extensions, wcs)

    # Compute Ra/Dec with new optimized WCS
    wcs_thistime = [None] * len(wcs)
    for i in range(1, len(wcs)):
        wcs_thistime[i] = wcs[i].copy()
        apply_fit_params(wcs_thistime[i].header, p_bestfit)
        wcs_thistime[i].updateFromHeader()
    ra_new, dec_new = odi2wcs(odi_x, odi_y, extensions, wcs_thistime)
    
    old_new = numpy.ones(shape=(odi_x.shape[0],8))
    old_new[:,0] = ref_ra[:]
    old_new[:,1] = ref_dec[:]
    old_new[:,2] = ra_orig[:]
    old_new[:,3] = dec_orig[:]
    old_new[:,4] = ra_new[:]
    old_new[:,5] = dec_new[:]
    old_new[:,6] *= wcs[2].header['CRVAL1']
    old_new[:,7] *= wcs[2].header['CRVAL2']

    return p_bestfit, old_new

    

def optimize_shift(n_matches, verbose=False, declination=0, debuglogfile=None, sum_radius=5):

    # Sort out all matches that only have one match
    not_just_one = n_matches[:,0] > 1
    
    n_multiple = n_matches[not_just_one]

    # print "\n\n\nReducing: ", n_matches.shape, n_multiple.shape,"\n\n\n\n\n"

    # Create memory to hold the alignment data
    # First two columns are dx, dy, third will be #matches
    total_sum = numpy.zeros(shape=(n_multiple.shape[0],4))
    total_sum[:,0:2] = n_multiple[:,1:3]

    if (debuglogfile != None):
        debuglog = open(debuglogfile, "w")
        numpy.savetxt(debuglog, n_matches)
        print("\n\n\n\n\n", file=debuglog)

        numpy.savetxt(debuglog, n_multiple)
        print("\n\n\n\n\n", file=debuglog)

        numpy.savetxt(debuglog, total_sum)
        print("\n\n\n\n\n", file=debuglog)


    # Go through the list of possible shifts and sum up all shifts in the neighborhood
    d_neighbor = sum_radius / 3600. # i.e. 5 arcseconds
    d_ra_max = d_neighbor / math.cos(math.radians(declination))
    d_dec_max = d_neighbor

    if (verbose): print("\n\n\nXXX Search-radius for match-conglomeration (arcsec): ",\
            d_ra_max*3600, d_dec_max*3600,"\n\n\n")

    for i in range(total_sum.shape[0]):
        d_ra = numpy.fabs(total_sum[i,0] - total_sum[:,0])
        d_dec = numpy.fabs(total_sum[i,1] - total_sum[:,1])
        
        nearby = (d_ra < d_ra_max) & (d_dec < d_dec_max)

        #nearby = (total_sum[:,0] > total_sum[i,0]-d_ra) & (total_sum[:,0] < total_sum[i,0]+d_ra) \
        #    & (total_sum[:,1] > total_sum[i,1]-d_dec) & (total_sum[:,1] < total_sum[i,1]+d_dec)
        total_sum[i,2] = numpy.sum(n_multiple[:,0][nearby])
        total_sum[i,3] = numpy.sum(nearby)

    if (debuglogfile != None):
        numpy.savetxt(debuglog, total_sum)
        print("\n\n\n\n\n", file=debuglog)

    # Sort the values to make it easy to find the maximum
    #si = numpy.argsort(total_sum[:,2])
    best = numpy.argmax(total_sum[:,2])

    best_guess = total_sum[best,:]
    if (verbose): print("max positions",best_guess)

    max_matches = best_guess[2]

    # To find some level of uncertainty, determine the maximum shift positions with 
    # a match-count above 10% of the maximum 
    
    within_fmax = total_sum[:,2] > 0.5*max_matches #0.1*max_matches
    dxmin, dxmax = numpy.min(total_sum[:,0][within_fmax]), numpy.max(total_sum[:,0][within_fmax])
    dymin, dymax = numpy.min(total_sum[:,1][within_fmax]), numpy.max(total_sum[:,1][within_fmax])
    dr = numpy.sqrt( (dxmax-dxmin)**2 + (dymax-dymin)**2 )

    if (verbose): 
        stdout_write("best offset: %.3f deg    %.3f deg\n" % (best_guess[0], best_guess[1]))
        stdout_write("Range: %.3f...%.3f deg    %.3f...%.3f deg\n" % (dxmin, dxmax, dymin, dymax))
        stdout_write("uncertainty: %.3f'' %.3f'' (%.3f)\n" % ( 0.5*3600.*(dxmax-dxmin), 0.5*3600*(dymax-dymin), 0.5*3600*dr ) )

    # determine some contrast level: peak number of matches in this 
    # peak relative to the maximum number of matches in the next-highest peak

    outside_peak = (total_sum[:,0] < best_guess[0]-2*d_neighbor) | (total_sum[:,0] > best_guess[0]+2*d_neighbor) | \
        (total_sum[:,1] < best_guess[1]-2*d_neighbor) | (total_sum[:,1] > best_guess[1]+2*d_neighbor)
    secondary_solutions = total_sum[outside_peak]
    second_peak_idx = numpy.argmax(secondary_solutions[:,2])
    second_peak = secondary_solutions[second_peak_idx,:]

    contrast = 1.0 * float(best_guess[2]) / float(second_peak[2])

    if (debuglogfile is not None):
        numpy.savetxt(debuglog, secondary_solutions)
  
    if (verbose): 
        print("secondary peak:", second_peak)
        print("contrast:", contrast)


    if (debuglogfile is not None):
        debuglog.close()

    return best_guess, contrast, [dr, dxmin, dxmax, dymin, dymax]



if __name__ == "__main__":

    starfinder = cmdline_arg_set_or_default("-starfind", "findstars")
    refcatalog = cmdline_arg_set_or_default("-refcatalog", "ippref")
    
    if (cmdline_arg_isset("-multi")):
        for infile in get_clean_cmdline()[1:]:
            outfile = infile[0:-5]+".wcs.fits"
            fixwcs(infile, outfile, starfinder, refcatalog)
    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        fixwcs(fitsfile, output_filename, starfinder, refcatalog)
