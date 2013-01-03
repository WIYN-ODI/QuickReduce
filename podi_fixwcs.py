#!/usr/bin/env python


import sys
import numpy
import os
from query_usno import query_usno
from podi_definitions import *
import pyfits
#import date
import datetime
import pywcs  
from astLib import astWCS
import pdb
import scipy
import scipy.stats


arcsec = 1./3600.

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
            if (math.fabs(shift_dx) > 0.05 or math.fabs(shift_dy) > 0.05):
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
    #print "best position", match_results[argmax

            #sys.exit(0)

    #for o in range(ota_count):
    #    for r in range(ref_count):
    #        print "RESULT",o,r,number_of_matches[o,r]

            
    


if __name__ == "__main__":

    fitsfile = sys.argv[1]
    hdulist = pyfits.open(fitsfile)
    ra, dec = hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']

    output_filename = sys.argv[2]

    # Obtain catalog listing from USNO-B1
    #usno = query_usno(10.74225833333333, 41.37409777777778, 45, 1000, "usno.cat")
    usno = query_usno(ra, dec, 45, 10000, "usno.cat", download=False)
    #usno = numpy.loadtxt("test.debug_ref")
    stdout_write("Read %s stars from USNO-B1\n" % (usno.shape[0]))
    

    # Run Sextractor on the datafile
    sex_catalogfile = "test.cat"
    sex_config_file = "sex.conf"
    sex_cmd = "sex -c %s %s" % (sex_config_file, fitsfile)
    print sex_cmd
    #os.system(sex_cmd)

    # Read the Sextractor output catalog
    #sex_catalogfile = "test.debug_ota"
    sex_cat = numpy.loadtxt(sex_catalogfile)
    stdout_write("Reading %d stars from Sextractor catalog\n" % sex_cat.shape[0])

    # Eliminate all stars with flags
    sex_cat = sex_cat[sex_cat[:,10] == 0]
    stdout_write("%d stars left after eliminating flags\n" % sex_cat.shape[0])

    #print "Read",sex_cat.shape[0],"detections from Sextractor\n"


    # Now go through each of the extensions
    for ext in range(1, len(hdulist)):
    #for ext in range(7, 8): 
        ota_cat = sex_cat[sex_cat[:,1] == ext]
        
        # Compute the median central position of this OTA
        #inherit_headers(hdulist[ext].header, hdulist[0].header)
        wcs = astWCS.WCS(hdulist[ext].header, mode="pyfits")
        ra, dec = wcs.getCentreWCSCoords()
        print ra,dec

        # Assume a (half-)width of 6 arcmin, that's quite a bit larger than the 4armin of one OTA
        width_ra = 0.1/math.cos(numpy.radians(dec))
        width_dec = 0.1
        print width_ra, width_dec

        print "RA=",usno[0:15,0]
        print "DEC=",usno[0:15,1]

        print "RA-Range:",ra-width_ra,ra+width_ra
        print "DEC-Range:",dec-width_dec,dec+width_dec
        nearby_ra  = (usno[:,0] > ra-width_ra) & (usno[:,0] < ra+width_ra) 
        nearby_dec = (usno[:,1] > dec-width_dec) & (usno[:,1] < dec+width_dec)
        
        print numpy.sum(nearby_ra), numpy.sum(nearby_dec)

        nearby = nearby_ra & nearby_dec
        nearby_usno = usno[nearby]
        print "Full USNO:",usno.shape[0],", nearby:",nearby_usno.shape[0]

        print ota_cat.shape[0],"stars on OTA",ext

        #for i in range(ota_cat.shape[0]):
        #    print "@@@@Sex",ota_cat[i,4], ota_cat[i,5],ota_cat[i,6], ota_cat[i,7]

        #for i in range(usno.shape[0]):
        #    print "@@@@Usno",usno[i,0], usno[i,1]

        ota_x, ota_y = ota_cat[:,6], ota_cat[:,7]
        ref_x, ref_y = nearby_usno[:,0], nearby_usno[:,1]

        dx, dy, n = shift_align_wcs(ota_x, ota_y, ref_x, ref_y)
        print "#######################"
        print "##"
        print "## OTA %d" % (ext)
        print "## dx,dy = %.5f, %.5f" % (dx, dy)
        print "## #matches = %d" % (n)
        print "##"
        print "#######################"

        hdulist[ext].header['CRVAL1'] += dx
        hdulist[ext].header['CRVAL2'] += dy

    hdulist.writeto(output_filename, clobber=True)
