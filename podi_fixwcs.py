#!/usr/bin/env python


import sys
import numpy
import os
from query_usno import query_usno
from podi_definitions import *
import pyfits
#import date
import datetime
#import pywcs  
#from astLib import astWCS
import pdb
import scipy
import scipy.stats


arcsec = 1./3600.
number_bright_stars = 100
max_offset = 0.1


N = 200
N_brightest_ota = 70
N_brightest_ref = 150

# Compute matches in smaller blocks to not blow up memory
# Improve: Change execution to parallel !!!
blocksize = 100


def shift_align_wcs(ota_x, ota_y, ref_x, ref_y):

    ota_count = ota_x.shape[0]
    ref_count = ref_x.shape[0]
    print ota_count, ref_count

    matching_radius = 3.6 * arcsec
    max_d2 = matching_radius * matching_radius

    print "Matching radius",matching_radius, max_d2
    
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

    for o in range(ota_count):
        for r in range(ref_count):
            
            # Take some random shift
            shift_dx = ref_x[r] - ota_x[o]
            shift_dy = ref_y[r] - ota_y[o]
            if (math.fabs(shift_dx) > max_offset or math.fabs(shift_dy) > max_offset):
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

            d2 = dx*dx + dy*dy
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
                    print "MATCH_%d"%(matched_pair),ota_x[idx_ota[match]], ota_y[idx_ota[match]], ref_x[idx_ref[match]], ref_y[idx_ref[match]], shift_dx, shift_dy

                matched_pair += 1
                pass

    #
    # Now that we have the counts for all pairs, find out some 
    # median background level of the number of (random) matches
    #
    good_datapoints = numpy.isfinite(match_results[:,:,0])

    selected_matches = match_results[good_datapoints]
    print "Selected matches=",selected_matches.shape

    if (selected_matches.shape[0] <= 1):
        # No matches found, return no shift and signal that the solution is invalid (Nmatch<0) 
        return 0, 0, -1
        
    median = numpy.median(selected_matches[:,0])
    std = numpy.std(selected_matches[:,0])
    sigma = scipy.stats.scoreatpercentile(selected_matches[:,0], 84) - median

    print "MATCH-RESULTS:"
    print "Random matches:",median,"+/-",std,"   (",sigma,")"
    most_matches_idx = numpy.argmax(selected_matches[:,0])
    print "max position", most_matches_idx
    most_matches = selected_matches[most_matches_idx,0]
    print "maximum matches:", most_matches

    #max_index = numpy.unravel_index(max_position, (match_results.shape[0:2]))
    #print max_index
    #print match_results[max_index[0],max_index[1],:]
    print selected_matches[most_matches_idx]

    #
    # Now to be sure, make sure the peak we found is significant
    #
    if (most_matches >= median + 3*sigma):
        print "Found seemingly significant match:"
    else:
        print "Found match, but not sure how reliable it is:"
    print selected_matches[most_matches_idx]
    print "Random matches: %d +/- %d (1sigma) (3-sigma: %d)" % (median, sigma, median+3*sigma)

    return selected_matches[most_matches_idx,1], selected_matches[most_matches_idx,2], selected_matches[most_matches_idx,0]


def refine_wcs_shift(ref_x, ref_y, ota_x, ota_y, best_guess, alignment_checkfile=None):

    matches = numpy.zeros(shape=(ota_x.shape[0],2))
    matched = numpy.zeros(shape=(ota_x.shape[0],4))
        #for i in range(ref_x.shape[0]):
        #    print "REF",ref_x[i], ref_y[i]
        #for i in range(ota_cat.shape[0]):
        #    print "SEX",ota_cat[i,6], ota_cat[i,7]

    start = 0
    while(start < ota_x.shape[0]):
        aligned_x = ota_x[start:start+blocksize] + best_guess[0]
        aligned_y = ota_y[start:start+blocksize] + best_guess[1]

        dx1 = aligned_x.reshape((aligned_x.shape[0], 1)).repeat(ref_x.shape[0], axis=1)
        dx2 = ref_x.reshape((1, ref_x.shape[0])).repeat(aligned_x.shape[0], axis=0)
        dx = dx1 - dx2

        dy1 = aligned_y.reshape((aligned_y.shape[0], 1)).repeat(ref_y.shape[0], axis=1)
        dy2 = ref_y.reshape((1, ref_y.shape[0])).repeat(aligned_y.shape[0], axis=0)
        dy = dy1 - dy2

        d2 = dx*dx + dy*dy

        closest = numpy.argmin(d2, axis=1)

        for i in range(closest.shape[0]):
            matches[start+i,0] = dx[i, closest[i]]
            matches[start+i,1] = dy[i, closest[i]]
            #print matches[start+i,:]

            matched[start+i,:] = [aligned_x[i], aligned_y[i], ref_x[closest[i]], ref_y[closest[i]]]

        start += blocksize

    # Now with all the matches for this OTA, left's determine the best solution
    # Again, use the 3-sigma clipping as above

    dmin, dmax = [-1, -1], [1,1]
    clipping = matches.copy()
    for rep in range(2):
        valid_range = (clipping[:,0] < dmax[0]) & (clipping[:,0] > dmin[0]) & (clipping[:,1] > dmin[1]) & (clipping[:,1] < dmax[1])
        print "valid_range.shape=",valid_range.shape
        print "count=",numpy.sum(valid_range)

        median = numpy.median(clipping[valid_range], axis=0)
        sigma_p = scipy.stats.scoreatpercentile(clipping[valid_range], 84)
        sigma_m = scipy.stats.scoreatpercentile(clipping[valid_range], 16)
        sigma = 0.5*(sigma_p - sigma_m)
        dmin = median - 3*sigma
        dmax = median + 3*sigma
        print "MEDIAN/SIGMA=",median, sigma, dmin, dmax

    # Now add all the offsets to the header

    if (alignment_checkfile != None):
        for i in range(matched.shape[0]):
            print >>alignment_check, \
                matched[i,0], matched[i,1], matched[i,2], matched[i,3],\
                ext,\
                matched[i,0]-median[0],matched[i,1]-median[1],matches[i,0], matches[i,1]
        print >>alignment_check, "\n\n\n\n\n"

    return median


def fixwcs(fitsfile, output_filename):

    tmp, dummy = os.path.split(sys.argv[0])
    dot_config_dir = tmp + "/.config/"
    print dot_config_dir
    
    hdulist = pyfits.open(fitsfile)
    ra, dec = hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']

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
    print "USNO dump-file:", usno_dumpfile
    usno_dump = open(usno_dumpfile, "w")
    for i in range(usno.shape[0]):
        print >>usno_dump, usno[i, 0], usno[i, 1]
    usno_dump.close()
    
    #
    # Run Sextractor on the input frame
    #
    sex_catalogfile = fitsfile[:-5]+".source.cat"
    sex_config_file = dot_config_dir+"fixwcs.conf"
    sex_param_file = dot_config_dir+"fixwcs.param"
    sex_logfile = fitsfile[:-5]+".sextractor.log"
    sex_cmd = "sex -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s >& %s" % (sex_config_file, sex_param_file, sex_catalogfile, fitsfile, sex_logfile)
    print sex_cmd
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
        print >>sex_dump, sex_cat[i, 6], sex_cat[i, 7]
    sex_dump.close()

    #
    # Now go through each of the extension
    # Improve: Change execution to parallel !!!
    #
    frame_shift = numpy.zeros(shape=(len(hdulist)-1, 2))
    for ext in range(1, len(hdulist)):
        #
        # Select a sub-catalog with stars in this OTA
        #
        ota_cat = sex_cat[sex_cat[:,1] == ext]
        print ota_cat.shape[0],"stars on OTA",ext

        # Select the N brightest stars to make alignment faster 
        # in the case of many stars, i.e. close to the milky way plane 
        # or in the case of resolved external galaxies
        use_N_brightest_ota = numpy.min([ota_cat.shape[0], N_brightest_ota])
        ota_magsort = numpy.argsort(ota_cat[:,2])
        ota_use = numpy.zeros(shape=(use_N_brightest_ota, ota_cat.shape[1]))
        for i in range(use_N_brightest_ota):
            ota_use[i,:] = ota_cat[ota_magsort[i], :]


        # Compute the median central position of this OTA
        #inherit_headers(hdulist[ext].header, hdulist[0].header)
        centerx, centery = hdulist[ext].header['NAXIS1']/2, hdulist[ext].header['NAXIS2']/2
        ra  = (centerx-hdulist[ext].header['CRPIX1'])*hdulist[ext].header['CD1_1'] + (centery-hdulist[ext].header['CRPIX2'])*hdulist[ext].header['CD1_2'] + hdulist[ext].header['CRVAL1']
        dec = (centerx-hdulist[ext].header['CRPIX1'])*hdulist[ext].header['CD2_1'] + (centery-hdulist[ext].header['CRPIX2'])*hdulist[ext].header['CD2_2'] + hdulist[ext].header['CRVAL2']
        print "center:", ra, dec
        #wcs = astWCS.WCS(hdulist[ext].header, mode="pyfits")
        #ra, dec = wcs.getCentreWCSCoords()
        #print "center of this ota:", ra,dec

        #
        # Now select stars that are within a given distance of this OTA
        # Assume a (half-)width of 6 arcmin, that's quite a bit larger than the 4armin of one OTA
        #
        width_ra = 0.1/math.cos(numpy.radians(dec))
        width_dec = 0.1
        #print width_ra, width_dec
        #print "RA=",usno[0:15,0]
        #print "DEC=",usno[0:15,1]
        print "RA-Range:",ra-width_ra,ra+width_ra
        print "DEC-Range:",dec-width_dec,dec+width_dec
        nearby_ra  = (usno[:,0] > ra-width_ra) & (usno[:,0] < ra+width_ra) 
        nearby_dec = (usno[:,1] > dec-width_dec) & (usno[:,1] < dec+width_dec)
        print numpy.sum(nearby_ra), numpy.sum(nearby_dec)
        nearby = nearby_ra & nearby_dec
        nearby_usno = usno[nearby]
        print "Full USNO:",usno.shape[0],", nearby:",nearby_usno.shape[0]

        if (nearby_usno.shape[0] <= 5):
            dx, dy, n = 0, 0, 0
        else:    
            #
            # In the catalog of nearby reference stars, select the N brightest stars
            #
            use_N_brightest_ref = numpy.min([nearby_usno.shape[0], N_brightest_ref])
            ref_magsort = numpy.argsort(nearby_usno[:, 4])
            ref_use = numpy.zeros(shape=(use_N_brightest_ref, nearby_usno.shape[1]))
            for i in range(use_N_brightest_ref):
                ref_use[i,:] = nearby_usno[ref_magsort[i], :]

            #for i in range(ota_cat.shape[0]):
            #    print "@@@@Sex",ota_cat[i,4], ota_cat[i,5],ota_cat[i,6], ota_cat[i,7]

            #for i in range(usno.shape[0]):
            #    print "@@@@Usno",usno[i,0], usno[i,1]

            #ota_x, ota_y = ota_cat[:,6], ota_cat[:,7]
            #ref_x, ref_y = nearby_usno[:,0], nearby_usno[:,1]

            ota_x, ota_y = ota_use[:,6], ota_use[:,7]
            ref_x, ref_y = ref_use[:,0], ref_use[:,1]

            dx, dy, n = shift_align_wcs(ota_x, ota_y, ref_x, ref_y)

        frame_shift[ext-1,:] = [dx, dy]

        print "#######################"
        print "##"
        print "## OTA %d" % (ext)
        print "## dx,dy = %.5f, %.5f" % (dx, dy)
        print "## #matches = %d" % (n)
        print "##"
        print "#######################"

    #
    # Now that the individual OTAs have returned their shift positions, 
    # compute the mean shift for the full frame
    #
    full_shift = numpy.median(frame_shift, axis=0)
    full_shift_std  = numpy.std(frame_shift, axis=0)
    print full_shift.shape

    print frame_shift
    print full_shift,full_shift_std

    #
    # Now that we have a pretty good guess on what the offsets are, 
    # extend the search on the full catalog and determine a better 
    # and more accurate solution!
    #


    #
    # Now do some 2-step iterative sigma clipping to get rid of outliers
    #
    dmin, dmax = [-1, -1], [1,1]
    clipping = frame_shift.copy()
    print full_shift,full_shift_std,dmin,dmax
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
        print median, sigma, dmin, dmax

    best_guess = median
    print "Best shift solution so far:",best_guess

    #
    # Now we have a pretty darn good estimate on the shift, and we also 
    # excluded outliers. Use this knowledge to derive an even better estimate 
    # by using the full input and reference catalogs
    #

    #
    # In an ideal world, the remainder should yield the same results for all OTAs, 
    # but for now keep them separate and check if that's the case
    #

    ref_x = usno[:,0]
    ref_y = usno[:,1]
    print ref_x.shape, ref_y.shape

    alignment_checkfile = cmdline_arg_set_or_default("-checkalign", None)
    alignment_check = None
    if (alignment_checkfile != None):
        alignment_check = open(alignment_checkfile, "w")

    if (cmdline_arg_isset("-singleota_shift")):
        for ext in range(1, len(hdulist)):

            # Extract just the source catalog for this OTA
            ota_cat = sex_cat[sex_cat[:,1]==ext]

            # Compute the average shift
            median = refine_wcs_shift(ref_x, ref_y, ota_cat[:,6], ota_cat[:,7], best_guess, alignment_check)

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
        median = refine_wcs_shift(ref_x, ref_y, sex_cat[:,6], sex_cat[:,7], best_guess, alignment_check)

        # Add the previous (best-guess) shift and the new refinement
        full_shift = best_guess + median

        # Apply refinement to all OTAs
        for ext in range(1, len(hdulist)):
            hdulist[ext].header['CRVAL1'] += full_shift[0]
            hdulist[ext].header['CRVAL2'] += full_shift[1]

            if (hdulist[ext].header['CRVAL1'] < 0): hdulist[ext].header['CRVAL1'] += 360 
            if (hdulist[ext].header['CRVAL1'] >360): hdulist[ext].header['CRVAL1'] -= 360 

            # Also add some history 
            hdulist[ext].header.add_history("LOCAL SHIFT: dRA,dDEC=")
            hdulist[ext].header.add_history("dRA=%.6f dDEC=%.6f" % (median[0], median[1]))

            hdulist[ext].header.add_history("GLOBAL SHIFT: dRA,dDEC=")
            hdulist[ext].header.add_history("dRA=%.6f dDEC=%.6f" % (best_guess[0], best_guess[1]))

    # All work done, close the check-file
    if (alignment_checkfile != None):
        alignment_check.close()

    if (not cmdline_arg_isset('-computeonly')):
        stdout_write("Writing output file...\n")
        hdulist.writeto(output_filename, clobber=True)
    else:
        stdout_write("Skipping the output file, you requested -computeonly\n")




if __name__ == "__main__":

    if (cmdline_arg_isset("-multi")):
        for infile in get_clean_cmdline()[1:]:
            outfile = infile[0:-5]+".wcs.fits"
            fixwcs(infile, outfile)
    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        fixwcs(fitsfile, output_filename)
