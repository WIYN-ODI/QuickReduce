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

import podi_matchcatalogs

arcsec = 1./3600.
number_bright_stars = 100
max_offset = 0.1


N = 200
N_brightest_ota = 70
N_brightest_ref = 150

# Compute matches in smaller blocks to not blow up memory
# Improve: Change execution to parallel !!!
blocksize = 100


def load_catalog_from_stripe82cat(ra, dec, calib_directory, sdss_filter):

    # Load the standard star catalog
    sdss_cat_filename = "%s/stripe82__%s.cat.npy" % (calib_directory, sdss_filter)
    stdout_write("Loading standard stars (%s)..." % (sdss_cat_filename))
    sdss_cat = numpy.load(sdss_cat_filename)
    stdout_write(" %d stars loaded!\n" % (sdss_cat.shape[0]))

    # Now take take of pointings close to RA of 0 hours
    if (ra < 5):
        stdout_write("Close to RA=0h, wrapping coordinates ...\n")
        sdss_cat[:,0][sdss_cat[:,0] > 180] -= 360.

    d_dec = 0.6
    d_ra  = d_dec / math.cos(math.radians(dec))
    
    # Break down the full catalog to stars nearby
    nearby = (sdss_cat[:,0] > ra-d_ra) & (sdss_cat[:,0] < ra+d_ra) & (sdss_cat[:,1] > dec-d_dec) & (sdss_cat[:,1] < dec+d_dec) 
    std_stars = sdss_cat[nearby]
    stdout_write("Found %d nearby (RA=%.1f, DEC=%.1f, +/- %.1fdeg) standard stars!\n" % (std_stars.shape[0], ra, dec, d_dec))

    """
    Output format is:
    Ra, Dec, mag_median, mag_mean, mag_std, mag_var
    """
    
    return std_stars


def load_catalog_from_sdss(ra, dec, sdss_filter, verbose=False):

    #import sqlcl

    #ra = 0
    
    if (numpy.array(ra).ndim > 0):
        min_ra = ra[0]
        max_ra = ra[1]
    else:
        min_ra = ra - 0.6/math.cos(math.radians(dec))
        max_ra = ra + 0.6/math.cos(math.radians(dec))
        
    if (min_ra < 0):
        ra_query = "ra > %(min_ra)f or ra < %(max_ra)f" % {"min_ra": min_ra+360, "max_ra": max_ra,} 
    else:
        ra_query = "ra BETWEEN %(min_ra)f and %(max_ra)f" % {"min_ra": min_ra, "max_ra": max_ra,} 
        
    if (numpy.array(dec).ndim > 0):
        min_dec = dec[0]
        max_dec = dec[1]
    else:
        min_dec = dec - 0.6
        max_dec = dec + 0.6

    #
    # This query is taken from the SDSS website and selects stars with clean photometry
    # --> http://skyserver.sdss3.org/dr8/en/help/docs/realquery.asp#cleanStars
    #
    sql_query = """\
SELECT ra,dec, u, err_u, g, err_g, r, err_r, i, err_i, z, err_z
FROM Star 
WHERE 
%(ra_query)s AND dec BETWEEN %(min_dec)f and %(max_dec)f
AND ((flags_r & 0x10000000) != 0)
-- detected in BINNED1
AND ((flags_r & 0x8100000c00a4) = 0)
-- not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
-- SATURATED, or BAD_COUNTS_ERROR
AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2))
-- not DEBLEND_NOPEAK or small PSF error
-- (substitute psfmagerr in other band as appropriate)
AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)
-- not INTERP_CENTER or not COSMIC_RAY
""" % {"filter": sdss_filter,
       "min_ra": min_ra, "max_ra": max_ra,
       "min_dec": min_dec, "max_dec": max_dec,
       "ra_query": ra_query,
       }

    # print sql_query
    
    # Taken from Tomas Budavari's sqlcl script
    # see http://skyserver.sdss3.org/dr8/en/help/download/sqlcl/default.asp 
    import urllib
    # Filter out comments starting with "--"
    fsql = ""
    for line in sql_query.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep;
    params = urllib.urlencode({'cmd': fsql, 'format': 'csv'})
    url = 'http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp'
    sdss = urllib.urlopen(url+'?%s' % params)
    # Budavari end


    answer = []
    for line in sdss:
        answer.append(line)
        if (((len(answer)-1)%10) == 0):
            stdout_write("\rFound %d stars so far ..." % (len(answer)-1))
    #answer = sdss.readlines()
    if (answer[0].strip() == "No objects have been found"):
        return numpy.zeros(shape=(0,0))
    else:
        stdout_write("\rFound a total of %d stars in SDSS catalog\n" % (len(answer)-1))
        
    if (verbose):
        print "Returned from SDSS:"
        print "####################################################"
        print ''.join(answer)
        print "####################################################"

    # If we are here, then the query returned at least some results.
    # Dump the first line just repeating what the output was
    del answer[0]

    
    print "Found %d results" % (len(answer))
    results = numpy.zeros(shape=(len(answer),12))
    # Results are comma-separated, so split them up and save as numpy array
    for i in range(len(answer)):
        items = answer[i].split(",")
        for col in range(len(items)):
            results[i, col] = float(items[col])
        #ra, dec = float(items[0]), float(items[1])
        #mag, mag_err =  float(items[2]), float(items[3])
        #results[i, :] = [ra, dec, mag, mag, mag_err, mag_err]
    
    return results
    
    
def photcalib_old(fitsfile, output_filename, calib_directory, overwrite_cat=None):

    tmp, dummy = os.path.split(sys.argv[0])
    dot_config_dir = tmp + "/.config/"
    print dot_config_dir
    
    hdulist = pyfits.open(fitsfile)
    ra, dec = hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']

    #
    # Run Sextractor on the input frame
    #
    sex_catalogfile = fitsfile[:-5]+".photcalib.cat"
    sex_config_file = dot_config_dir+"fixwcs.conf"
    sex_param_file = dot_config_dir+"fixwcs.param"
    sex_logfile = fitsfile[:-5]+".sextractor.log"
    sex_cmd = "sex -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s >& %s" % (sex_config_file, sex_param_file, sex_catalogfile, fitsfile, sex_logfile)
    print sex_cmd
    stdout_write("Running SExtractor to search for stars, be patient (logfile: %s) ..." % (sex_logfile))
    if ((not cmdline_arg_isset("-skip_sex")) or (not os.path.isfile(sex_catalogfile)) or cmdline_arg_isset("-forcesex")):
        os.system(sex_cmd)
    else:
        stdout_write("Re-using previous source catalog\n")
    stdout_write(" done!\n")

    stdout_write("\nPreparing work ...\n\n")
    
    # Read the Sextractor output catalog
    sex_cat = numpy.loadtxt(sex_catalogfile)
    stdout_write("Reading %d stars from Sextractor catalog\n" % sex_cat.shape[0])

    # Eliminate all stars with flags
    sex_cat = sex_cat[sex_cat[:,10] == 0]
    stdout_write("%d stars left after eliminating flags\n" % sex_cat.shape[0])


    # Figure out which SDSS to use for calibration
    filter = hdulist[0].header['FILTER']
    sdss_filter = sdss_equivalents[filter]
    
    std_stars = numpy.array([])
    if (overwrite_cat != None):
        # Read this catalog instead of any of the default ones:
        std_stars = numpy.loadtxt(overwrite_cat, delimiter=";")
    else:
        if (calib_directory != None):
            std_stars = load_catalog_from_stripe82cat(ra, dec, calib_directory, sdss_filter)
            print std_stars.shape
        if (std_stars.shape[0] <= 0):
            if (calib_directory != None):
                stdout_write("Couldn't find any stars in the Stripe82 Standard star catalog :(\n")

            stdout_write("Trying to get one directly from SDSS, please wait!\n\n")
            std_stars = load_catalog_from_sdss(ra, dec, sdss_filter)

            if (std_stars.shape[0] <= 0):
                stdout_write("No stars not found - looks like this region isn't covered by SDSS - sorry!\n\n")
                sys.exit(0)
                return
        
    dump = open("dump.cat", "w")
    for i in range(std_stars.shape[0]):
        print >>dump, std_stars[i,0], std_stars[i,1]
    print >>dump, "\n\n\n\n\n\n\n\n"
    for i in range(sex_cat.shape[0]):
        print >>dump, sex_cat[i,6], sex_cat[i,7]
    
    dump.close()

    #
    # Now go through each of the extension
    # Improve: Change execution to parallel !!!
    #
    stdout_write("\nStarting work, results in %s ...\n\n" % output_filename)
    results = open(output_filename, "w")

    # Now go through the reference catalog and search for matches.
    # To speed things up and avoid looking for matches between distant stars, break down the area into
    # a number of smaller blocks (scanning the full area in 5x5 blocks)
    ref_ra_min, ref_ra_max = numpy.min(std_stars[:,0]), numpy.max(std_stars[:,0])
    ref_dec_min, ref_dec_max = numpy.min(std_stars[:,1]), numpy.max(std_stars[:,1])

    ra_ranges = numpy.arange(ref_ra_min, ref_ra_max, 5)
    dec_ranges = numpy.arange(ref_dec_min, ref_dec_max, 5)

    dummy, ra_ranges = numpy.histogram([], bins=5, range=(ref_ra_min, ref_ra_max))
    dummy, dec_ranges = numpy.histogram([], bins=5, range=(ref_dec_min, ref_dec_max))

    # Reorganize the SExtractor oiutput file to have the RA/DEC values be in the first two columns
    # New order is then Ra, Dec, mag, mag_err, x, y, ota
    sex_reorg = numpy.zeros(shape=(sex_cat.shape[0], 7))
    sex_reorg[:,0:2] = sex_cat[:,6:8]
    sex_reorg[:,2:6] = sex_cat[:,2:6]
    sex_reorg[:,6] = sex_cat[:,1]

    # Select one of the ranges and hand the parameters off to the matching routine
    matched_cat = podi_matchcatalogs.match_catalogs(std_stars, sex_reorg)
            
    if (matched_cat != None):
        numpy.savetxt(results, matched_cat, delimiter=" ")
                
    results.close()






def photcalib(source_cat, output_filename, filtername, diagplots=True, calib_directory=None, overwrite_cat=None):

    # Figure out which SDSS to use for calibration
    sdss_filter = sdss_equivalents[filtername]
    if (sdss_filter == None):
        # This filter is not covered by SDSS, can't perform photometric calibration
        return None

    pc = sdss_photometric_column[sdss_filter]


    # Eliminate all stars with flags
    flags = source_cat[:,7]
    no_flags_set = (flags == 0)

    source_cat = source_cat[no_flags_set]

    ra_min = numpy.min(source_cat[:,0])
    ra_max = numpy.max(source_cat[:,0])
                       
    dec_min = numpy.min(source_cat[:,1])
    dec_max = numpy.max(source_cat[:,1])
                       
    stdout_write("\nPreparing work ...\n\n")
    

    std_stars = numpy.array([])
    if (overwrite_cat != None):
        # Read this catalog instead of any of the default ones:
        if (os.path.isfile(overwrite_cat)):
            print "reading stdstars from file"
            std_stars = numpy.loadtxt(overwrite_cat) #, delimiter=";")

    if (std_stars.shape[0] <= 0):
        if (calib_directory != None):
            std_stars = load_catalog_from_stripe82cat(ra, dec, calib_directory, sdss_filter)
            print std_stars.shape
        if (std_stars.shape[0] <= 0):
            if (calib_directory != None):
                stdout_write("Couldn't find any stars in the Stripe82 Standard star catalog :(\n")

            stdout_write("Trying to get one directly from SDSS, please wait!\n\n")
            #std_stars = load_catalog_from_sdss(ra, dec, sdss_filter)
            std_stars = load_catalog_from_sdss([ra_min, ra_max], [dec_min, dec_max], sdss_filter)

            if (std_stars.shape[0] <= 0):
                stdout_write("No stars not found - looks like this region isn't covered by SDSS - sorry!\n\n")
                sys.exit(0)
                return
        
            numpy.savetxt("stdstars", std_stars)

#    dump = open("dump.cat", "w")
#    for i in range(std_stars.shape[0]):
#        print >>dump, std_stars[i,0], std_stars[i,1]
#    print >>dump, "\n\n\n\n\n\n\n\n"
#    for i in range(sex_cat.shape[0]):
#        print >>dump, sex_cat[i,6], sex_cat[i,7]
#    
#    dump.close()

    #
    # Now go through each of the extension
    # Improve: Change execution to parallel !!!
    #
    stdout_write("\nStarting work, results in %s ...\n\n" % output_filename)
    results = open(output_filename, "w")

    odi_sdss_matched = podi_matchcatalogs.match_catalogs(source_cat, std_stars)

    if (odi_sdss_matched != None):
        numpy.savetxt(results, odi_sdss_matched, delimiter=" ")

    # Stars without match in SDSS have RA=-9999, let's sort them out
    found_sdss_match = odi_sdss_matched[:,2] >= 0
    
    odi_sdss_matched = odi_sdss_matched[found_sdss_match]

    odi_ra, odi_dec = odi_sdss_matched[:,0], odi_sdss_matched[:,1]
    sdss_ra, sdss_dec = odi_sdss_matched[:,2], odi_sdss_matched[:,3]

    # Use photometry for the 3'' aperture
    odi_mag, odi_magerr = odi_sdss_matched[:,17], odi_sdss_matched[:,25]
    sdss_mag = odi_sdss_matched[:,(source_cat.shape[1]+pc)]
    sdss_magerr = odi_sdss_matched[:,(source_cat.shape[1]+pc+1)]

    # Determine the zero point
    zp = sdss_mag - odi_mag
    zperr = numpy.hypot(sdss_magerr, odi_magerr)
    import podi_collectcells
    zp_clipped = podi_collectcells.three_sigma_clip(zp)
    zp_median = numpy.median(zp_clipped)
    zp_std = numpy.std(zp_clipped)
    print "zeropoint (clipped)",zp_median," +/-", zp_std

    zp_upper1sigma = scipy.stats.scoreatpercentile(zp_clipped, 84)
    zp_lower1sigma = scipy.stats.scoreatpercentile(zp_clipped, 16)
    print zp_lower1sigma, zp_upper1sigma, 0.5*(zp_upper1sigma-zp_lower1sigma)

    zp_median_ = numpy.median(zp)
    zp_std_ = numpy.std(zp)
    print "zeropoint (un-clipped)",zp_median_," +/-", zp_std_

    # Make plots
    if (diagplots):
        import podi_diagnosticplots

        zp_calib_plot = output_filename[:-5]+".photZP.png"
        podi_diagnosticplots.photocalib_zeropoint(odi_mag, odi_magerr, sdss_mag, sdss_magerr,
                                                  zp_calib_plot,
                                                  zp_median, zp_std,
                                                  "r", "odi_r",
                                                  title="this is a test"
                                                  )

                
    results.close()

    return zp_median, zp_std, odi_sdss_matched



    
if __name__ == "__main__":

    
    if (cmdline_arg_isset("-multi")):
        calibdir = get_cmdline_arg("-calib")
        for infile in get_clean_cmdline()[1:]:
            output_filename = infile[0:-5]+".photcalib.dat"
            if (os.path.isfile(output_filename)):
                continue
            photcalib(infile, output_filename, calibdir)

    elif (cmdline_arg_isset("-new")):
        catalogfile = get_clean_cmdline()[1]
        source_cat = numpy.loadtxt(catalogfile)

        filtername = "odi_r"

        output_filename = get_clean_cmdline()[2]
        photcalib(source_cat, output_filename, filtername, diagplots=True, calib_directory=None, overwrite_cat="stdstars")
    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        calibdir = get_clean_cmdline()[3]

        photcalib(fitsfile, output_filename, calibdir)
