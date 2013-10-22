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
from astLib import astWCS
import pdb
import scipy
import scipy.stats

import podi_matchcatalogs
import podi_sitesetup

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


def load_catalog_from_sdss(ra, dec, sdss_filter, verbose=False, return_query=False, max_catsize=-1):

    #import sqlcl

    #ra = 0
    #print "# Loading catalog from SDSS online..."
    
    if (numpy.array(ra).ndim > 0):
        min_ra = ra[0]
        max_ra = ra[1]
    else:
        min_ra = ra - 0.6/math.cos(math.radians(dec))
        max_ra = ra + 0.6/math.cos(math.radians(dec))
        
    if (min_ra < 0):
        ra_query = "( ra > %(min_ra)f or ra < %(max_ra)f )" % {"min_ra": min_ra+360, "max_ra": max_ra,} 
    elif (max_ra > 360):
        ra_query = "( ra > %(min_ra)f or ra < %(max_ra)f )" % {"min_ra": min_ra, "max_ra": max_ra-360.,} 
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

    if (verbose): print sql_query
    
    stdout_write("Downloading catalog from SDSS ...")

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
        if (max_catsize > 0 and len(answer) >= max_catsize):
            break
        answer.append(line)
        if (((len(answer)-1)%10) == 0):
            if (verbose): stdout_write("\rFound %d stars so far ..." % (len(answer)-1))
    #answer = sdss.readlines()
    if (answer[0].strip() == "No objects have been found"):
        stdout_write(" nothing found\n")
        if (return_query):
            return numpy.zeros(shape=(0,12)), fsql #sql_query
        return numpy.zeros(shape=(0,12))

    stdout_write(" found %d stars!\n" % (len(answer)-1))
        
    if (verbose):
        print "Returned from SDSS:"
        print "####################################################"
        print ''.join(answer)
        print "####################################################"

    # If we are here, then the query returned at least some results.
    # Dump the first line just repeating what the output was
    del answer[0]

    
    if (verbose): print "Found %d results" % (len(answer))
    results = numpy.zeros(shape=(len(answer),12))
    # Results are comma-separated, so split them up and save as numpy array
    for i in range(len(answer)):
        items = answer[i].split(",")
        for col in range(len(items)):
            results[i, col] = float(items[col])
        #ra, dec = float(items[0]), float(items[1])
        #mag, mag_err =  float(items[2]), float(items[3])
        #results[i, :] = [ra, dec, mag, mag, mag_err, mag_err]
    
    if (return_query):
        return results, fsql #sql_query
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



def load_sdss_catalog_from_fits(sdss_ref_dir, ra_range, dec_range, verbose=True):

    #print "# Loading catalog from SDSS offline fits catalog..."

    # Load the SkyTable so we know in what files to look for the catalog"
    skytable_filename = "%s/SkyTable.fits" % (sdss_ref_dir)
    skytable_hdu = pyfits.open(skytable_filename)

    skytable = skytable_hdu['SKY_REGION'].data
    
    # Select entries that match our list
    if (verbose): print "# Searching for stars in sky with ra=%s, dec=%s" % (ra_range, dec_range)

    min_dec = dec_range[0]
    max_dec = dec_range[1]
    min_ra = ra_range[0]
    max_ra = ra_range[1]

    if (verbose): print min_ra, max_ra, min_dec, max_dec

    if (max_ra > 360.):
        # This wraps around the high end, shift all ra values by -180
        # Now all search RAs are ok and around the 180, next also move the catalog values
        selected = skytable['R_MIN'] < 180
        skytable['R_MAX'][selected] += 360
        skytable['R_MIN'][selected] += 360
    if (min_ra < 0):
        # Wrap around at the low end
        selected = skytable['R_MAX'] > 180
        skytable['R_MAX'][selected] -= 360
        skytable['R_MIN'][selected] -= 360

    if (verbose): print "# Search radius: RA=%.1f ... %.1f   DEC=%.1f ... %.1f" % (min_ra, max_ra, min_dec, max_dec)
    
    needed_catalogs = (skytable['R_MAX'] > min_ra)  & (skytable['R_MIN'] < max_ra) & \
                      (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)

    if (verbose): print skytable[needed_catalogs]
    
    files_to_read = skytable['NAME'][needed_catalogs]

    # Now we are with the skytable catalog, so close it
    skytable_hdu.close()
    del skytable

    #
    # Load all frames, one by one, and select all stars in the valid range.
    # Then add them to the catalog with RAs and DECs
    #
    catalog_columns = ['RA', 'DEC',
                       'MAG_U', 'MAGERR_U',
                       'MAG_G', 'MAGERR_G',
                       'MAG_R', 'MAGERR_R',
                       'MAG_I', 'MAGERR_I',
                       'MAG_Z', 'MAGERR_Z',
                       ]

    full_catalog = numpy.zeros(shape=(0,len(catalog_columns)))
    for catalogname in files_to_read:

        if (verbose): print "Reading 2mass catalog"
        catalogfile = "%s/%s.fits" % (sdss_ref_dir, catalogname)
        hdu_cat = pyfits.open(catalogfile)
        if (hdu_cat[1].header['NAXIS2'] <= 0):
            hdu_cat.close()
            continue

        # Read the RA and DEC values
        cat_ra  = hdu_cat[1].data.field('RA')
        cat_dec = hdu_cat[1].data.field('DEC')

        # To select the right region, shift a temporary catalog
        cat_ra_shifted = cat_ra
        if (max_ra > 360.):
            cat_ra_shifted[cat_ra < 180] += 360
        elif (min_ra < 0):
            cat_ra_shifted[cat_ra > 180] -= 360

        in_search_range = (cat_ra_shifted > min_ra) & (cat_ra_shifted < max_ra ) & (cat_dec > min_dec) & (cat_dec < max_dec)

        array_to_add = numpy.zeros(shape=(numpy.sum(in_search_range),len(catalog_columns)))
        if (verbose): print catalogfile, numpy.sum(in_search_range)
        for col in range(len(catalog_columns)):
            array_to_add[:,col] = hdu_cat[1].data.field(catalog_columns[col])[in_search_range]

        full_catalog = numpy.append(full_catalog, array_to_add, axis=0)
        
    if (verbose): print "# Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read))
    return full_catalog




def query_sdss_catalog(ra_range, dec_range, sdss_filter, verbose=False):

    #
    # Get a photometric reference catalog one way or another
    # 
    std_stars = None

    if (podi_sitesetup.sdss_ref_type == "stripe82"):
        std_stars = numpy.zeros(shape=(0,0)) #load_catalog_from_stripe82cat(ra, dec, calib_directory, sdss_filter)
        print std_stars.shape

    elif (podi_sitesetup.sdss_ref_type == 'local'):
        
        std_stars = load_sdss_catalog_from_fits(podi_sitesetup.sdss_ref_dir, ra_range, dec_range,
                                                verbose=verbose)
        
    elif (podi_sitesetup.sdss_ref_type == 'web'):
        #stdout_write("Trying to get one directly from SDSS, please wait!\n\n")
        #std_stars = load_catalog_from_sdss(ra, dec, sdss_filter)
        std_stars = load_catalog_from_sdss(ra_range, dec_range, sdss_filter,
                                           verbose=verbose)

    #print "returning",std_stars.shape[0],"stars..."
    return std_stars





def photcalib(source_cat, output_filename, filtername, exptime=1, 
              diagplots=True, calib_directory=None, overwrite_cat=None,
              plottitle=None, otalist=None,
              options=None):

    error_return_value = (99.9, 99.9, None, 99.9)

    # Figure out which SDSS to use for calibration
    sdss_filter = sdss_equivalents[filtername]
    if (sdss_filter == None):
        # This filter is not covered by SDSS, can't perform photometric calibration
        return error_return_value

    pc = sdss_photometric_column[sdss_filter]


    # Eliminate all stars with flags
    flags = source_cat[:,7]
    no_flags_set = (flags == 0)

    source_cat = source_cat[no_flags_set]

    ra_min = numpy.min(source_cat[:,0])
    ra_max = numpy.max(source_cat[:,0])
    # Make sure we deal with RAs around 0 properly
    if (math.fabs(ra_max - ra_min) > 180):
        ra_min, ra_max = ra_max, ra_min+360.

    dec_min = numpy.min(source_cat[:,1])
    dec_max = numpy.max(source_cat[:,1])
                 
      
    stdout_write("\nPreparing work ...\n\n")
    
    ra_range  = [ra_min, ra_max]
    dec_range = [dec_min, dec_max]

    std_stars = query_sdss_catalog(ra_range, dec_range, sdss_filter)

    if (std_stars.shape[0] <= 0):
        stdout_write("No stars not found - looks like this region isn't covered by SDSS - sorry!\n\n")
        return error_return_value


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

    zp_correction_exptime = -2.5 * math.log10(exptime)
    odi_mag -= zp_correction_exptime

    # Determine the zero point
    zp = sdss_mag - odi_mag 
    zperr = numpy.hypot(sdss_magerr, odi_magerr)
    #import podi_collectcells
    if (zp.shape[0] <= 0):
        zp_median = 99
        zp_std = 1
        zp_exptime = 99
        diagplots = False
    else:
        if (zp.shape[0] <=5):
            zp_clipped = zp
        else:
            zp_clipped = three_sigma_clip(zp)

        zp_upper1sigma = scipy.stats.scoreatpercentile(zp_clipped, 84)
        zp_lower1sigma = scipy.stats.scoreatpercentile(zp_clipped, 16)
        print zp_lower1sigma, zp_upper1sigma, 0.5*(zp_upper1sigma-zp_lower1sigma)

        zp_median = numpy.median(zp_clipped)
        zp_std = numpy.std(zp_clipped)
        print "zeropoint (clipped)",zp_median," +/-", zp_std

        zp_median_ = numpy.median(zp)
        zp_std_ = numpy.std(zp)
        zp_exptime = zp_median + zp_correction_exptime

    print "zeropoint (un-clipped)",zp_median_," +/-", zp_std_

    # Make plots
    if (diagplots):
        import podi_diagnosticplots

        # zp_calib_plot = output_filename[:-5]+".photZP"
        zp_calib_plot = create_qa_filename(output_filename, "photZP", options)
        podi_diagnosticplots.photocalib_zeropoint(odi_mag, odi_magerr, sdss_mag, sdss_magerr,
                                                  zp_calib_plot,
                                                  zp_median, zp_std,
                                                  sdss_filter, filtername, #"r", "odi_r",
                                                  title=plottitle,
                                                  options=options,
                                                  also_plot_singleOTAs=options['otalevelplots'])

        print odi_sdss_matched[0,:]

        ota = odi_sdss_matched[:,10]
        ra = odi_sdss_matched[:,0]
        dec = odi_sdss_matched[:,1]
        # plotfilename = output_filename[:-5]+".photZP_map"
        plotfilename = create_qa_filename(output_filename, "photZP_map", options)

        ota_outlines = None
        if (otalist != None):
            ota_outlines = derive_ota_outlines(otalist)

        podi_diagnosticplots.photocalib_zeropoint_map(odi_mag, sdss_mag, ota, ra, dec,
                                                      output_filename=plotfilename,
                                                      sdss_filtername=sdss_filter, odi_filtername=filtername,
                                                      title=plottitle,
                                                      ota_outlines=ota_outlines,
                                                      options=options,
                                                      also_plot_singleOTAs=options['otalevelplots'])

    results.close()

    return zp_median, zp_std, odi_sdss_matched, zp_exptime



    
if __name__ == "__main__":

    if (cmdline_arg_isset("-multi")):
        calibdir = get_cmdline_arg("-calib")
        for infile in get_clean_cmdline()[1:]:
            output_filename = infile[0:-5]+".photcalib.dat"
            if (os.path.isfile(output_filename)):
                continue
            photcalib(infile, output_filename, calibdir)

    elif (cmdline_arg_isset("-new")):
        fitsfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(fitsfile)
        catalogfile = fitsfile+".src.cat"
        source_cat = numpy.loadtxt(catalogfile)
        
        filtername = hdulist[0].header['FILTER']
        print catalogfile,"-->",source_cat.shape

        output_filename = get_clean_cmdline()[2]

        from podi_collectcells import read_options_from_commandline
        options = read_options_from_commandline(None)

        photcalib(source_cat, output_filename, filtername, exptime=100,
                  diagplots=True, calib_directory=None, overwrite_cat="stdstars",
                  otalist=hdulist,
                  options=options)

    elif (cmdline_arg_isset("-querysdss")):
        ra_min = float(get_clean_cmdline()[1])
        ra_max = float(get_clean_cmdline()[2])
        dec_min = float(get_clean_cmdline()[3])
        dec_max = float(get_clean_cmdline()[4])

        ra_range = [ra_min, ra_max]
        dec_range = [dec_min, dec_max]
        catalog = query_sdss_catalog(ra_range, dec_range, "r", verbose=False)

        if (not catalog == None):
            #print catalog.shape
            pass
        else:
            print "there was a problem..."

        if (cmdline_arg_isset("-print")):
            numpy.savetxt(sys.stdout, catalog)
        else:
            print "Found",catalog.shape[0],"results"
        #print catalog

    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        calibdir = get_clean_cmdline()[3]

        photcalib(fitsfile, output_filename, calibdir)
