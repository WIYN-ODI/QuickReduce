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
    
    min_ra = ra - 0.6/math.cos(math.radians(dec))
    max_ra = ra + 0.6/math.cos(math.radians(dec))
    if (min_ra < 0):
        ra_query = "ra > %(min_ra)f or ra < %(max_ra)f" % {"min_ra": min_ra+360, "max_ra": max_ra,} 
    else:
        ra_query = "ra BETWEEN %(min_ra)f and %(max_ra)f" % {"min_ra": min_ra, "max_ra": max_ra,} 
        
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
    
    
def photcalib(fitsfile, output_filename, calib_directory, overwrite_cat=None):

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

    for cur_ra in range(ra_ranges.shape[0]-1):
        for cur_dec in range(dec_ranges.shape[0]-1):

            # Select one of the ranges and hand the parameters off to the matching routine
            matched_cat = match_catalogs(std_stars, sex_reorg, ra_ranges[cur_ra:cur_ra+2], dec_ranges[cur_dec:cur_dec+2])
            
            if (matched_cat != None):
                numpy.savetxt(results, matched_cat, delimiter=" ")
                
    results.close()


def match_catalogs(ref_full, odi_full, ra_ranges, dec_ranges, matching_radius=2):


    # First of all extract only the reference stars in the selected box
    ref_select = (ref_full[:,0] >  ra_ranges[0]) & (ref_full[:,0] <  ra_ranges[1]) \
               & (ref_full[:,1] > dec_ranges[0]) & (ref_full[:,1] < dec_ranges[1])

    ref = ref_full[ref_select]

    if (ref.shape[0] <= 0):
        return None

    # Now do the same, just with ranges extended by the matching radius / cos(dec)
    max_cos_dec = numpy.max(numpy.fabs(ref[:,1]))
    d_ra  = matching_radius * arcsec * math.cos(math.radians(max_cos_dec))
    d_dec = matching_radius * arcsec
    
    odi_select = (odi_full[:,0] >  ra_ranges[0]-d_ra ) & (odi_full[:,0] <  ra_ranges[1]+d_ra ) \
               & (odi_full[:,1] > dec_ranges[0]-d_dec) & (odi_full[:,1] < dec_ranges[1]+d_dec)
    odi = odi_full[odi_select]

    print "Searching for reference stars: RA: %3.1f - %3.1f,  DEC: %+4.1f - %+4.1f" % (ra_ranges[0]-d_ra, ra_ranges[1]+d_ra,dec_ranges[0]-d_dec,dec_ranges[1]+d_dec),
    # Create the output array, and set all ODI position to the "not found" standard value
    output_array = numpy.zeros(shape=(ref.shape[0], ref.shape[1]+odi.shape[1]))
    output_array[:, 2:4] = -99999.9

    # Copy all values from the reference catalog
    odi_data_start = 4 + ref.shape[1] - 2 
    output_array[:, 0:2] = ref[:, 0:2]
    output_array[:, 4:odi_data_start] = ref[:, 2:]

    matches_found = 0
    if (odi.shape[0] > 0):
        #
        # Now for each star, find the closest match in the reference catalog
        #
        matching_radius = 2 * arcsec
        matching_radius_squared = matching_radius * matching_radius
        for star in range(ref.shape[0]):

            d_ra  = (ref[star,0] - odi[:,0]) * math.cos(math.radians(ref[star,1]))
            d_dec = (ref[star,1] - odi[:,1])
            sq_d = d_ra * d_ra + d_dec * d_dec

            si = numpy.argsort(sq_d)
            if (sq_d[si[0]] < matching_radius_squared):
                # The closest match is a valid one, with distance < matching_radius
                # Copy the ODI coordinates and catalog data into the output format
                output_array[star, 2:4] = odi[si[0], 0:2]
                output_array[star, odi_data_start:] = odi[si[0], 2:]
                matches_found += 1

    print "  -->    %5d matches (%5d vs %5d)" % (matches_found, ref.shape[0], odi.shape[0])
    
    return output_array

    
if __name__ == "__main__":

    
    if (cmdline_arg_isset("-multi")):
        calibdir = get_cmdline_arg("-calib")
        for infile in get_clean_cmdline()[1:]:
            output_filename = infile[0:-5]+".photcalib.dat"
            if (os.path.isfile(output_filename)):
                continue
            photcalib(infile, output_filename, calibdir)
    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        calibdir = get_clean_cmdline()[3]

        photcalib(fitsfile, output_filename, calibdir)
