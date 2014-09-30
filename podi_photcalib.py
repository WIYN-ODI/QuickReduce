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
This module contains all functionality for performing a photometric calibration
for a given frame. It does so by matching the ODI source catalog to a catalog
of stars from the SDSS. The comparison between instrumental and reference 
magnitudes then yields the calibration zeropoint.

If requested, podi_photcalib also creates diagnostic plots that allow to easily
judge the quality of the obtained calibration.

The module also contains all required functionality to access local SDSS star
catalogs or, if no local copies are available, query the SDSS online.

Command line options
--------------------

* **-multi**

  Deprecated for now, do not use

* **-new**

* **-querysdss**

  Test-routine that queries the SDSS and prints the returned catalog

  Run as:

      ``podi_photcalib -querysdss ra_min ra_max dec_min dec_max (-print)``

  If the -print option is given as well print the entire catalog to stdout,
  otherwise only print the number of returned stars.


* **-swarp**

  Compute the photometric zeropoint for a frame created with swarp.

  Run as:

      ''podi_photcalib.py -swarp swarped_frame.fits swarp_input_frame.fits``

  swarp_input_frame here is one of the frames that went into the swarped frame.
  This is necessary to obtain the FILTER and EXPTIME keywords that are required
  for the photometric calibration, but are typically not contained in the 
  swarped frame any longer.

* **-aucap**

  Compute the photometric zeropoint for a frame created with the Automatic
  Calibration Pipeline (AuCaP).

  Run as:
 
    ``podi_photcalib.py -aucap file1.fits file2.fits file3.fits ...``

  In this mode, podi_photcalib can run on any number of input frames in sequence.

Note for the -aucap and -swarp modi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To save time during programming, two additional flags are supported in these 
modi. By default, if the source catalog for the specified input file already 
exists, it is simply re-used and not created. To force re-running SExtractor,
add the -resex flag. Additionally, in the ``-aucap`` mode, if the source catalog
exists, podi_photcalib also does not re-create the diagnostic plots. To change
this behavior and re-create the plots, add the ``-replot`` flag.

Methods
-------

"""


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

import podi_matchcatalogs
import podi_sitesetup as sitesetup
import podi_search_ipprefcat

import podi_logging
import logging
import time
import subprocess

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
    """
    Load the source catalog from the custom pODI Stripe82 catalog.

    Parameters
    ----------
    ra, dec -- float

        Center position of the search area. The search radius is set to 0.6 deg.

    calib_directory - string

        directory that contains the Stripe82 catalog

    sdss_filter - string

        Name of the filter for which the magnitudes should be returned.


    Returns
    -------
    A catalog of sources within the search area, stored in an array. Columns
    are as follows. Ra, Dec, mag_median, mag_mean, mag_std, mag_var

    """

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


    
    
def photcalib_old(fitsfile, output_filename, calib_directory, overwrite_cat=None):
    """
    Deprecated, do not use.
    """

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
    """
    Retrieve the SDSS catalog for the specified area, using the local FITS-based
    catalog of the SDSS.

    Parameters
    ----------

    sdss_ref_dir : string

        Directory containing the local copy of the SDSS star catalog. This 
        directory has to contain the catalog index in a file named 
        "SkyTable.fits"

    ra_range : float[2] (e.g. [165.5, 165.9])

        Minimum and maximum RA range, in degrees

    dec_range : float[2]

        as ra_range, just for declination

    verbose : Bool

        add some debugging output useful for error tracking


    Returns
    -------
   
    source-catalog

        contains the following columns: ra, dec, mag_u, magerr_u, g, r, i, z

    """
    #print "# Loading catalog from SDSS offline fits catalog..."
    logger = logging.getLogger("ReadSDSSCatalogFromFits")

    # Load the SkyTable so we know in what files to look for the catalog"
    skytable_filename = "%s/SkyTable.fits" % (sdss_ref_dir)
    skytable_hdu = pyfits.open(skytable_filename)

    skytable = skytable_hdu['SKY_REGION'].data
    
    # Select entries that match our list
    
    logger.debug("# Searching for stars in sky with ra=%s, dec=%s" % (ra_range, dec_range))
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
    logger.debug(files_to_read)

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
    logger.debug("# Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read)))

    return full_catalog




def query_sdss_catalog(ra_range, dec_range, sdss_filter, verbose=False):
    """
    High-level interface to the SDSS catalog search. Depending on the settings
    in podi_sitesetup (see sdss_ref_type) this forwards the query to search 
    either the local Stripe82 catalog, the online SDSS catalog or the local FITS
    catalog.

    Parameters
    ----------

    sdss_ref_dir : string

        Directory containing the local copy of the SDSS star catalog. This 
        directory has to contain the catalog index in a file named 
        "SkyTable.fits"

    ra_range : float[2] (e.g. [165.5, 165.9])

        Minimum and maximum RA range, in degrees

    dec_range : float[2]

        as ra_range, just for declination

    verbose : Bool

        add some debugging output useful for error tracking


    Returns
    -------
   
    source-catalog

    """

    #
    # Get a photometric reference catalog one way or another
    # 
    std_stars = None

    if (sitesetup.sdss_ref_type == "stripe82"):
        std_stars = numpy.zeros(shape=(0,0)) #load_catalog_from_stripe82cat(ra, dec, calib_directory, sdss_filter)
        print std_stars.shape

    elif (sitesetup.sdss_ref_type == 'local'):
        
        std_stars = load_sdss_catalog_from_fits(sitesetup.sdss_ref_dir, ra_range, dec_range,
                                                verbose=verbose)
        
    elif (sitesetup.sdss_ref_type == 'web'):
        #stdout_write("Trying to get one directly from SDSS, please wait!\n\n")
        #std_stars = load_catalog_from_sdss(ra, dec, sdss_filter)
        std_stars = load_catalog_from_sdss(ra_range, dec_range, sdss_filter,
                                           verbose=verbose)

    #print "returning",std_stars.shape[0],"stars..."
    return std_stars





def photcalib(source_cat, output_filename, filtername, exptime=1, 
              diagplots=True, calib_directory=None, overwrite_cat=None,
              plottitle=None, otalist=None,
              options=None,
              verbose=False,
              eliminate_flags=True,
              matching_radius=3,
              detailed_return={}):

    """
    Perform the photometric calibration, create the diagnostic plots and return
    the results.

    Parameters
    ----------
    source_cat : ndarray

        Input source array containing the catalog output from the SExtractor 
        run.

    output_filename : string

        Filename of the file containing the resulting image file. In the case of
        collectcells, this is the output filename, hence the name. However, if 
        photcalib is run on existing images, this is the image INPUT filename.

    filtername : string

        Name of the ODI filter. Based on this filter name photcalib will decide
        the corresponding SDSS filter.

    exptime : float

        Exposure time of the frame, required to normalize the photometric zero-
        point.

    diagplots : Bool

        If True, create the diagnostic plots for the photometric calibration.

    calib_directory : string

    overwrite_cat : Bool

    plottitle : string

        Title of the diagnostic plot, e.g. containing the filename, or some
        information about the Object, etc.

    otalist : HDUList

        If given, photcalib uses the WCS information from the individual
        extensions to add OTA outlines to the global, focal-plane wide diagnostic
        plot.

    options : dictionary

        General options, containing information about some command line flags.

    eliminate_flags : Bool

        If set to True, only sources with no SourceExtractor flags (that might
        indicate problems with the photometry) are used for the calibration. If 
        set, fewer stars are used but these stars are of the best possible 
        quality.

    matching_radius : float

        Matching radius to be used when matching the ODI source catalog to the
        SDSS reference catalog.


    Returns
    -------

    ZP-median : float

        median zeropoint across the entire focal plane

    ZP-Std : float

        standard deviation of the photometric zeropoint

    ODI-SDSS-matched : ndarray

        matched catalog of ODI sources and SDSS reference stars.

    ZP-normalized : float

        photometric zeropoint, corrected and normalized to an exposure time
        of 1 second.

    """

    logger = logging.getLogger("PhotCalib")

    error_return_value = (99.9, 99.9, None, 99.9)
    detailed_return['median'] = -99.
    detailed_return['std'] = -99.
    detailed_return['zp_exptime'] = -99.
    detailed_return['stderrofmean'] = -99.
    detailed_return['zp_upper1sigma'] = -99.
    detailed_return['zp_lower1sigma'] = -99.
    detailed_return['n_clipped'] = -1
    detailed_return['n_raw'] = -1
    detailed_return['colorterm'] = None
    detailed_return['colorcorrection'] = None
    detailed_return['radialZPfit'] = None
    detailed_return['radialZPfit_error'] = None
    detailed_return['zp_restricted'] = None
    detailed_return['zp_magnitude_slope'] = None
    detailed_return['aperture_mode'] = None
    detailed_return['aperture_size'] = -99.
    detailed_return['aperture_mag'] = None
    detailed_return['aperture_magerr'] = None
    detailed_return['catalog'] = "none"
    detailed_return['reference_filter'] = "none"

    # Eliminate all stars with flags
    if (eliminate_flags):
        flags = source_cat[:,SXcolumn['flags']]
        no_flags_set = (flags == 0)
        source_cat = source_cat[no_flags_set]

    # numpy.savetxt("photcat", source_cat)

    ra_min = numpy.min(source_cat[:,SXcolumn['ra']])
    ra_max = numpy.max(source_cat[:,SXcolumn['ra']])
    logger.debug("ra-range: %f --- %f" % (ra_min, ra_max))
    # Make sure we deal with RAs around 0 properly
    if (math.fabs(ra_max - ra_min) > 180):
        ra_fixed = source_cat[:,SXcolumn['ra']].copy()
        ra_fixed[ra_fixed > 180] -= 360
        ra_min = numpy.min(ra_fixed)
        ra_max = numpy.max(ra_fixed)
        logger.debug("ra-range (fixed): %f --- %f" % (ra_min, ra_max))

    dec_min = numpy.min(source_cat[:,SXcolumn['dec']])
    dec_max = numpy.max(source_cat[:,SXcolumn['dec']])
                 
    logger.debug("Starting photometric calibration!")
    #stdout_write("\nPreparing work ...\n\n")
    
    ra_range  = [ra_min, ra_max]
    dec_range = [dec_min, dec_max]

    sdss_filter = sdss_equivalents[filtername]
    logger.debug("Translating filter: %s --> %s" % (filtername, sdss_filter))
    if (sdss_filter == None):
        # This filter is not covered by SDSS, can't perform photometric calibration
        return error_return_value

    std_stars = None
    for catname in sitesetup.photcalib_order:

        if (catname == "sdss"):
            # Figure out which SDSS to use for calibration
            logger.info("Querying the SDSS catalog")
            _std_stars = podi_search_ipprefcat.get_reference_catalog(ra_range, dec_range, 
                                                                    radius=None,
                                                                    basedir=sitesetup.sdss_ref_dir,
                                                                    cattype="sdss")
            if (not _std_stars == None and _std_stars.shape[0] > 0):
                detailed_return['catalog'] = "SDSS"
                # Found some SDSS stars
                std_stars = _std_stars
                logger.debug("Found %d stars in the SDSS" % (std_stars.shape[0]))
                break
            else:
                logger.debug("No stars not found - looks like this region isn't covered by SDSS - sorry!")

        elif (catname == "ucac4"):

            detailed_return['catalog'] = "UCAC4"

            # Try to get some data from UCAC4 instead
            if (not sitesetup.ucac4_ref_dir == None and 
                not sitesetup.ucac4_ref_dir == "none" and
                os.path.isdir(sitesetup.ucac4_ref_dir)):
                logger.debug("Running UCAC4 query (%s)" % (sitesetup.ucac4_ref_dir))
                _std_stars = podi_search_ipprefcat.get_reference_catalog(
                    ra=ra_range, 
                    dec=dec_range,
                    radius=None,
                    basedir=sitesetup.ucac4_ref_dir,
                    cattype='ucac4',
                    verbose=False
                )


                # Sort out all stars with invalid magnitudes
                valid = (_std_stars[:,2] < 20) & (numpy.fabs(_std_stars[:,3]) <= 0.9)
                logger.debug("Number of UCAC4 stars: %s" % (str(_std_stars.shape)))
                # numpy.savetxt("ucac.dump", _std_stars)

                std_stars = _std_stars[valid]
                # numpy.savetxt("ucac.dump2", std_stars)
                logger.debug("Found a valid UCAC4 source catalog (%s)" % (str(std_stars.shape)))
                break
            else:
                logger.warning("Problem with querying the UCAC4 catalog")

        elif (catname == "ippref"):
            logger.debug("Running IPPRef query (%s)" % (sitesetup.ippref_ref_dir))
            _std_stars = podi_search_ipprefcat.get_reference_catalog(ra_range, dec_range, 
                                                                     radius=None,
                                                                     basedir=sitesetup.ippref_ref_dir,
                                                                     cattype="IPPRef")
            if (not _std_stars == None and _std_stars.shape[0] > 0):
                detailed_return['catalog'] = "IPPRef"
                std_stars = _std_stars
                logger.debug("Found %d reference sources in IPPRef catalog" % (std_stars.shape[0]))
                break
            else:
                logger.warning("Problem with querying the IPPRef catalog")
                
        else:
            logger.error("Illegal catalog name in the photcalib order: %s" % (catname))
            return error_return_value


    #
    # Now go through each of the extension
    # Improve: Change execution to parallel !!!
    #
    logger.debug("Starting work, results in %s ..." % output_filename)
    # results = open(output_filename+".photcal", "w")

    odi_sdss_matched = podi_matchcatalogs.match_catalogs(source_cat, std_stars, matching_radius=matching_radius)
    # numpy.savetxt("odisource.dump", source_cat)
    # numpy.savetxt("matched.dump", odi_sdss_matched)

    # Stars without match in SDSS have RA=-9999, let's sort them out
    found_sdss_match = odi_sdss_matched[:,2] >= 0
    odi_sdss_matched = odi_sdss_matched[found_sdss_match]
    n_matches = numpy.sum(found_sdss_match)

    # numpy.savetxt("matched.dump2", odi_sdss_matched)
    logger.debug("Found %d matches between ODI source catalog and photometric reference catalog" % (
            n_matches))

    if (n_matches <= 0):
        logger.warning("No overlap between source and photometric reference catalog found!")
        return error_return_value

    odi_ra, odi_dec = odi_sdss_matched[:,0], odi_sdss_matched[:,1]
    sdss_ra, sdss_dec = odi_sdss_matched[:,2], odi_sdss_matched[:,3]

    #
    # Use photometry for the 3'' aperture
    #
    # The +2 in the index is because match_catalog adds the reference 
    # coordinates as columns 3 & 4
    
    # NEW: use a user-definable aperture
    photcalib_odi_aperture = "fwhm_x_3" #"auto"

    try:
        if (photcalib_odi_aperture == "fwhm_x_3"):
            logger.debug("Using seeing-based photcalib aperture")
            # If set to auto, base the magnitude on the seeing to be > 2.5 * FWHM
            seeing = numpy.median(odi_sdss_matched[:, SXcolumn['fwhm_world']+2]) * 3600.
            logger.debug("Found seeing = %.2f" % (seeing))
            detailed_return['aperture_mode'] = "fwhm_x_3"
            for i in range(len(SXapertures)):
                if (SXapertures[i] > seeing * 3 or i == (len(SXapertures)-1)):
                    col_mag = 'mag_aper_%0.1f' % (SXapertures[i])
                    col_magerr = 'mag_err_%0.1f' % (SXapertures[i])
                    detailed_return['aperture_size'] = SXapertures[i]
                    break
        elif (photcalib_odi_aperture == "auto"):
            logger.debug("Choosing AUTO_MAG for photometric calibration")
            col_mag = 'mag_auto'
            col_magerr = 'mag_err_auto'
            detailed_return['aperture_mode'] = "mag_auto"
            detailed_return['aperture_size'] = 0.
        elif (photcalib_odi_aperture in SXapertures):
            logger.debug("Choosing a numeric aperture size: %.1f" % (photcalib_odi_aperture))
            col_mag = 'mag_aper_%0.1f' % (photcalib_odi_aperture)
            col_magerr = 'mag_err_%0.1f' % (photcalib_odi_aperture)
            detailed_return['aperture_mode'] = "manual"
            detailed_return['aperture_size'] = photcalib_odi_aperture
        elif (type(photcalib_odi_aperture) == tuple and len(photcalib_odi_aperture) == 2):
            logger.debug("Using user-defined aperture/error for phot.calib: %s" % (str(photcalib_odi_aperture)))
            col_mag, col_magerr = photcalib_odi_aperture
            detailed_return['aperture_mode'] = "manual"
            detailed_return['aperture_size'] = 0
        else:
            logger.debug("Using default aperture of 6.0 arcsec")
            col_mag, col_magerr = 'mag_aper_6.0', 'mag_err_6.0'
            detailed_return['aperture_mode'] = "default"
            detailed_return['aperture_size'] = 6.0
    except:
        podi_logging.log_exception()
        logger.error("Problem with configuring the phot.calib aperture, using 6.0'' diameter instead")
        detailed_return['aperture_mode'] = "default-error"
        detailed_return['aperture_size'] = 6.0
        col_mag, col_magerr = 'mag_aper_6.0', 'mag_err_6.0'

    logger.debug("Aperture for photcalib: %s %s" % (col_mag, col_magerr))
    detailed_return['aperture_mag'] = col_mag
    detailed_return['aperture_magerr'] = col_magerr

    # only select stars with properly detected magnitudes in the chosen aperture
    valid_odi_mag = (odi_sdss_matched[:,SXcolumn[col_mag]+2] < 75)
    odi_sdss_matched = odi_sdss_matched[valid_odi_mag]

    odi_mag = odi_sdss_matched[:,SXcolumn[col_mag]+2]
    odi_magerr = odi_sdss_matched[:,SXcolumn[col_magerr]+2]

    detailed_return['odi_sdss_matched_raw'] = odi_sdss_matched.copy()

    photref_col_mag = photref_col_err = -1

    if (detailed_return['catalog'] == "SDSS"):
        # Compute the calibration magnitude from SDSS, 
        # accounting for color-terms if needed
        pc = sdss_photometric_column[sdss_filter]
        detailed_return['reference_filter'] = sdss_filter

        photref_col_mag = source_cat.shape[1]+pc
        photref_col_err = source_cat.shape[1]+pc+1

        sdss_mag = odi_sdss_matched[:,(source_cat.shape[1]+pc)]
        sdss_magerr = odi_sdss_matched[:,(source_cat.shape[1]+pc+1)]

        if (filtername in photzp_colorterms):
            logger.debug("Found color-term definition for this filter (%s)" % (filtername))
            colorterm, filter1, filter2 = photzp_colorterms[filtername]

            col1 = sdss_photometric_column[filter1]
            col2 = sdss_photometric_column[filter2]
            sdss_color = odi_sdss_matched[:,(source_cat.shape[1]+col1)] - odi_sdss_matched[:,(source_cat.shape[1]+col2)]
            color_correction = colorterm * sdss_color

            odi_sdss_matched[:,(source_cat.shape[1]+pc)] -= color_correction
            
            detailed_return['colorterm'] = colorterm
            detailed_return['colorcorrection'] = "sdss_%s - sdss_%s" % (filter1, filter2)
        else:    
            detailed_return['colorterm'] = None
            logger.debug("No color-term definition for this filter (%s)" % (filtername))
    elif (detailed_return['catalog'] == "UCAC4"):
        
        logger.info("Using UCAC for photometric calibration.")
        # For test-purposes, always assume the g-band filter
        pc = 2

        photref_col_mag = source_cat.shape[1]+pc
        photref_col_err = source_cat.shape[1]+pc+1

        detailed_return['reference_filter'] = "UCAC-Red"
        sdss_mag = odi_sdss_matched[:,(source_cat.shape[1]+pc)]
        sdss_magerr = odi_sdss_matched[:,(source_cat.shape[1]+pc+1)]

    elif (detailed_return['catalog'] == "IPPRef"):
        
        logger.info("Using IPPRef for photometric calibration.")
        ipp_ref_columns = {"g": 2, "r": 4, "i": 6, "z": 8}
        if (sdss_filter in ipp_ref_columns):
            pc = ipp_ref_columns[sdss_filter]
            sdss_mag = odi_sdss_matched[:,(source_cat.shape[1]+pc)]
            sdss_magerr = odi_sdss_matched[:,(source_cat.shape[1]+pc+1)]

            photref_col_mag = source_cat.shape[1]+pc
            photref_col_err = source_cat.shape[1]+pc+1
        else:
            logger.debug("No reference photometry for this filter (%s -> %s) found in IPPRef" % (
                filtername, sdss_filter))
    else:
        return error_return_value

    detailed_return['odi_sdss_matched'] = odi_sdss_matched.copy()

    if (photref_col_mag < 0):
        logger.debug("No Valid reference photometry found!")
        return error_return_value
        
    logger.debug("ODI/SDSS: %s %s" % (str(odi_mag.shape), str(sdss_mag.shape)))

    zp_correction_exptime = -2.5 * math.log10(exptime)
    odi_mag -= zp_correction_exptime

    # Determine the zero point
    zp = (sdss_mag - odi_mag)
    zperr = numpy.hypot(sdss_magerr, odi_magerr)

    #
    # If requested, remove radial ZP trend
    # 
    if (not otalist == None
        and options['fitradialZP']):

        logger.debug("Starting to fit the radial ZP gradient")

        # get Ra/Dec center position
        center_dec = otalist[1].header['CRVAL2'] - (240./3600)
        cos_declination = math.cos(math.radians(center_dec))
        center_ra  = otalist[1].header['CRVAL1'] + ((240./3600) / cos_declination)
        logger.debug("Using center position %f %f" % (center_ra, center_dec))

        # Convert all positions into radial coordinates
        r = numpy.sqrt(
            ((odi_sdss_matched[:,0] - center_ra) * cos_declination)**2 + 
            (odi_sdss_matched[:,1] - center_dec)**2 )

        # Now we have a and we have ZP
        median_ZP = numpy.median(zp)
        if (numpy.sum(r > 0.3) > 15):
            median_ZP = numpy.median(zp[r > 0.3])

        # select all valid ZP points within 0.3 degrees of center
        valid_points = (r < 0.3) & (numpy.fabs(zp-median_ZP) < 0.3)

        #
        # Now we have radius r and we have ZP
        #

        def zptrend(p, r):
            return p[0] + p[1]*r

        def zptrend_err(p, r, zp, zp_err):
            correction = zptrend(p, r)
            diff = zp - correction
            return diff/zp_err

        # Perform a linear least-square fit to the data
        pinit = [0.,0.]
        args = (r[valid_points], (zp-median_ZP)[valid_points], zperr[valid_points])
        fit = scipy.optimize.leastsq(zptrend_err, pinit, args=args, full_output=1)
        pfit = fit[0]
        uncert = numpy.sqrt(numpy.diag(fit[1]))

        logger.debug("leastsq results: ZP = %.4f + %.4f * r" % (pfit[0], pfit[1]))
        logger.debug("lastesq uncerts: %.4f / %.4f" % (uncert[0], uncert[1]))

        # Create the best-fit polynomial
        # if (numpy.sum(valid_points) > 3):
        #     p4 = numpy.polyfit (r[valid_points], (zp-median_ZP)[valid_points], 1)
        # else:
        #     p4 = (0, 0)
        # poly = numpy.poly1d(p4)
        # logger.debug("Found offset = %.3f mag, slope = %.4f mag/deg" % (p4[0], p4[1])) 
        # print p4

        # And subtract the fit from the data
        # zp[valid_points] -= poly(r[valid_points])
        
        detailed_return['radialZPfit'] = pfit
        detailed_return['radialZPfit_error'] = uncert

        # Add some plotting here if needed
        import podi_plotting, matplotlib, matplotlib.pyplot
        fig = matplotlib.pyplot.figure()
        ax = [fig.add_subplot(211), fig.add_subplot(212)]
        for i in range(2):
            ax[i].set_ylim((median_ZP-0.3, median_ZP+0.3))
            ax[i].set_xlim((-0.02, 0.6))
            ax[i].set_ylabel("photometric ZP")
        ax[1].set_xlabel("distance from center [deg]")
        ax[0].set_title("Photometric zeropoint gradient -- before & after")
        # plot original data as red + signs
        ax[0].scatter(r, zp, c='red', marker='+')
        # plot the corrected data as blue circles
        # zp_fixed = zp[valid_points]-poly(r[valid_points])

        r_fixed = r[valid_points]
        zp_fixed = zp[valid_points]
        zp_fixed -= zptrend(pfit, r_fixed)
        ax[1].scatter(r_fixed,  zp_fixed, c='blue', marker='o', s=10, linewidths=0)
        ax[1].scatter(r[r >= 0.3],  zp[r >= 0.3], c='blue', marker='o', s=10, linewidths=0)

        # Add fit as black, solid line
        fit_line = numpy.linspace(0,0.3,200)
        # ax.plot(fit_line, poly(fit_line)+median_ZP, 'k-')
        ax[0].plot(fit_line, zptrend(pfit, fit_line)+median_ZP, 'k-', linewidth=2)

        median_line_x = numpy.linspace(-0.1, 1., 200)
        median_line_y = numpy.zeros_like(median_line_x) + median_ZP
        ax[0].plot(median_line_x, median_line_y, 'g-')
        ax[1].plot(median_line_x, median_line_y, 'g-')
        zp_radial_plot = create_qa_filename(output_filename, "photZP_radial", options)+".png"
        logger.debug("Saving plot to %s" % (zp_radial_plot))
        fig.savefig(zp_radial_plot)


    # Set some default values for most return results
    zp_median = 99
    zp_std = 1
    zp_exptime = 99
    zp_stderrofmean = 0.0
    zp_upper1sigma = 99.
    zp_lower1sigma = 99.
    n_clipped = 0

    zp = (sdss_mag - odi_mag)
    zperr = numpy.hypot(sdss_magerr, odi_magerr)

    detailed_return['odi_col_mag'] = SXcolumn[col_mag]+2
    detailed_return['odi_col_err'] = SXcolumn[col_magerr]+2
    detailed_return['photref_col_mag'] = photref_col_mag
    detailed_return['photref_col_err'] = photref_col_err

    detailed_return['odi_sdss_matched'] = odi_sdss_matched
    detailed_return['odi_sdss_matched_smallerrors'] = odi_sdss_matched
    detailed_return['odi_sdss_matched_clipped'] = odi_sdss_matched

    # print "\n"*10,detailed_return,"\n"*5


    if (odi_sdss_matched.shape[0] > 0):

        #
        # Now select only sources with reference uncertainties <~ 0.02 mag
        # The detailed cutoff can be configured via the sitesetup.py configuration
        #
        small_reference_errors = odi_sdss_matched[:, photref_col_err] \
                            < sitesetup.photcalib_error_cutoff[detailed_return['catalog']]
        if (numpy.sum(small_reference_errors) > 5):
            small_errors = odi_sdss_matched[small_reference_errors].copy()
            large_errors = odi_sdss_matched[~small_reference_errors].copy()
        else:
            small_errors = odi_sdss_matched.copy()
            large_errors = numpy.zeros((0,odi_sdss_matched.shape[1]))

        detailed_return['odi_sdss_matched_smallerrors'] = small_errors
        detailed_return['odi_sdss_matched_largeerrors'] = large_errors

        #
        # From the optimized catalog, re-extract magnitudes and erors for ODI 
        # and the photometric reference catalog
        #
        sdss_mag = small_errors[:, photref_col_mag]
        sdss_magerr = small_errors[:, photref_col_err]
        odi_mag = small_errors[:, SXcolumn[col_mag]+2]
        odi_magerr = small_errors[:, SXcolumn[col_magerr]+2]
        zp = small_errors[:, photref_col_mag] - small_errors[:, SXcolumn[col_mag]+2]
        zp_err = numpy.hypot(small_errors[:, photref_col_err], small_errors[:, SXcolumn[col_magerr]+2])

        #
        # Reject outliers by iterative sigma-clipping
        #
        detailed_return['odi_sdss_matched_clipped'] = small_errors
        if (zp.shape[0] <=5):
            clipping_mask = numpy.ones((zp.shape), dtype=numpy.bool)
        else:
            _dum, clipping_mask = three_sigma_clip(zp, return_mask=True)

        zp_clipped = zp[clipping_mask]
        zperr_clipped = zperr[clipping_mask]
        sdss_mag_clipped = sdss_mag[clipping_mask]

        detailed_return['odi_sdss_matched_clipped'] = small_errors[clipping_mask].copy()
        detailed_return['odi_sdss_matched_ref'] = small_errors[clipping_mask].copy()
        detailed_return['odi_sdss_matched_outlier'] = small_errors[~clipping_mask].copy()

        zp_upper1sigma = scipy.stats.scoreatpercentile(zp_clipped, 84)
        zp_lower1sigma = scipy.stats.scoreatpercentile(zp_clipped, 16)
        # print zp_lower1sigma, zp_upper1sigma, 0.5*(zp_upper1sigma-zp_lower1sigma)

        zp_median = numpy.median(zp_clipped)
        zp_std = numpy.std(zp_clipped)
        # print "zeropoint (clipped)",zp_median," +/-", zp_std

        zp_median_ = numpy.median(zp)
        zp_std_ = numpy.std(zp)
        zp_exptime = zp_median - zp_correction_exptime

        # compute the r.m.s. value of the distribution as well
        final_cat = small_errors[clipping_mask]
        zp_uncert = numpy.hypot(zp_err, 0.001)
        rms2 = numpy.sum((zp - zp_median)**2/zp_uncert) / numpy.sum(1./zp_uncert)
        detailed_return['rms'] = numpy.sqrt(rms2)
        detailed_return['sem'] = scipy.stats.sem(zp)

        # print "zeropoint (un-clipped)",zp_median_," +/-", zp_std_

        #
        # Now try to only use the top 100 brightest stars with errors < 0.1 mag
        #
        good_photometry = odi_magerr < 0.1
        if (numpy.sum(good_photometry) > 3):
            brightest_star_mag = numpy.min(odi_mag[odi_magerr < 0.1])
            si = numpy.argsort(odi_mag)
            number_within_3mags = numpy.sum( ((odi_mag < brightest_star_mag+3) & (odi_magerr < 0.1)) )
            sel_median, sel_std, sel_psigma, sel_msigma, sel_n = -99., -99., -99., -99., -1
            if (number_within_3mags > 3):
                if (number_within_3mags > 100): number_within_3mags = 100
                zp_sel = numpy.zeros( (number_within_3mags, 4) )
                for i in range(number_within_3mags):
                    zp_sel[i,:] = [odi_mag[si[i]], sdss_mag[si[i]], zp[si[i]], zperr[si[i]]]
                # Now we have only a subset of all points
                sel_median = numpy.median(zp_sel[:,2])
                sel_std = numpy.std(zp_sel[:,2])
                sel_psigma = scipy.stats.scoreatpercentile(zp_sel[:,2], 84)
                sel_msigma = scipy.stats.scoreatpercentile(zp_sel[:,2], 16)
                sel_medodimag = numpy.median(zp_sel[:,0])
                sel_maxodimag = numpy.max(zp_sel[:,0])
                sel_minodimag = numpy.min(zp_sel[:,0])
                sel_n = number_within_3mags
                detailed_return['zp_restricted'] = (sel_median, sel_std, sel_psigma, sel_msigma, sel_n, sel_medodimag, sel_maxodimag, sel_minodimag)
        else:
            detailed_return['zp_restricted'] = None

        #
        # Also fit a slope to the full data set. This way, if the slope is
        # significantly larger than zero this means trouble
        #
        def linear_fit(p, odi_mag):
            return p[0] + p[1] * odi_mag
        def linear_fit_err(p, odi_mag, zp, zp_err):
            linear = linear_fit(p, odi_mag)
            return (linear - zp) / zp_err

        if (zp_clipped.shape[0] > 10):
            p_init = [zp_median, 0]

            # enlarge the error bars a bit to prevent single points with 
            # unrealistically small uncertainties to dominate teh fit
            zperr_larger = numpy.hypot(zperr_clipped, 0.05)

            args = (sdss_mag_clipped, zp_clipped, zperr_larger)
            fit = scipy.optimize.leastsq(linear_fit_err, p_init, args=args, full_output=1)
            pfit = fit[0]
            uncert = numpy.sqrt(numpy.diag(fit[1]))
            detailed_return['zp_magnitude_slope'] = (pfit, uncert)
            logger.debug("Clipped zp-slope: %.3f (+/- %.3f) + (%.5f +/- %.5f) * sdss_mag" % (pfit[0], uncert[0], pfit[1], uncert[1]))

            zp_stderrofmean = scipy.stats.sem(zp_clipped)
            n_clipped = zp_clipped.shape[0]
        else:
            detailed_return['zp_magnitude_slope'] = None

    detailed_return['median'] = zp_median
    detailed_return['std'] = zp_std
    detailed_return['zp_exptime'] = zp_exptime
    detailed_return['stderrofmean'] = zp_stderrofmean
    detailed_return['zp_upper1sigma'] = zp_upper1sigma
    detailed_return['zp_lower1sigma'] = zp_lower1sigma
    detailed_return['n_clipped'] = n_clipped
    detailed_return['n_raw'] = zp.shape[0]

    # Make plots
    logger.debug("Creating phot-ZP diagnostic plots? %s" % (str(diagplots)))
    if (diagplots):
        import podi_diagnosticplots

        # zp_calib_plot = output_filename[:-5]+".photZP"
        logger.debug("Preparing the photZP plots")
        zp_calib_plot = create_qa_filename(output_filename, "photZP", options)
        try:
            podi_diagnosticplots.photocalib_zeropoint(output_filename=zp_calib_plot,
                                                      sdss_filtername=sdss_filter, 
                                                      odi_filtername=filtername,
                                                      title=plottitle,
                                                      options=options,
                                                      also_plot_singleOTAs=options['otalevelplots'],
                                                      details=detailed_return)
        except:
            podi_logging.log_exception()
            logger.error("Problem with creating the photometric calibration diagnostic plot")

        # print odi_sdss_matched[0,:]

        ota = odi_sdss_matched[:,10]
        ra = odi_sdss_matched[:,0]
        dec = odi_sdss_matched[:,1]
        # plotfilename = output_filename[:-5]+".photZP_map"
        plotfilename = create_qa_filename(output_filename, "photZP_map", options)

        ota_outlines = None
        if (otalist != None):
            ota_outlines = derive_ota_outlines(otalist)

        logger.debug("Preparing the photZP_map plots")
        try:
            podi_diagnosticplots.photocalib_zeropoint_map(
                details=detailed_return,
                output_filename=plotfilename,
                sdss_filtername=sdss_filter, odi_filtername=filtername,
                title=plottitle,
                ota_outlines=ota_outlines,
                options=options,
                also_plot_singleOTAs=options['otalevelplots'])

        except:
            podi_logging.log_exception()
            logger.error("Problem with creating the photometric zeropoint map")

    # results.close()

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

        sys.exit(0)

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

    elif (cmdline_arg_isset("-swarp")):

        print """\

   ************************************************
   * photcalib for - Swarp'ed QR data             *
   *               - single 
   * part of the QuickReduce pipeline package     *
   * (c) 2014, Ralf Kotulla and WIYN Observatory  *
   ************************************************
"""
        import podi_collectcells
        options = podi_collectcells.read_options_from_commandline()
        # Setup everything we need for logging
        podi_logging.setup_logging(options)

        logger = logging.getLogger("SwarpPhotCalib")
        recreate_catalogs = cmdline_arg_isset("-resex")

        # Disable OTA-level plots
        options['otalevelplots'] = False

        for inputfile in get_clean_cmdline()[1:]:
            try:
                logger.info("Starting phot-calib (%s)..." % (inputfile))

                # Run SourceExtractor
                logger.info("Starting source-extractor")
                sex_config_file = "%s/.config/wcsfix.sex" % (options['exec_dir'])
                parameters_file = "%s/.config/wcsfix.sexparam" % (options['exec_dir'])
                catfile = "%s.cat" % (inputfile[:-5])
                sexcmd = "%(sex)s -c %(config)s -PARAMETERS_NAME %(params)s -CATALOG_NAME %(cat)s %(fits)s" % {
                    'sex': sitesetup.sextractor, 
                    'config': sex_config_file,
                    'params': parameters_file,
                    'cat': catfile, 
                    'fits': inputfile,
                }

                weight_file = inputfile[:-5]+".weight.fits"
                if (os.path.isfile(weight_file)):
                    logger.info("Using weight file for source extraction (%s)" % (weight_file))
                    sexcmd += " -WEIGHT_TYPE MAP_WEIGHT -RESCALE_WEIGHTS Y -WEIGHT_IMAGE %s" % (weight_file)

                logger.debug(sexcmd)
                if (os.path.isfile(catfile) and not recreate_catalogs):
                    logger.info("catalog exists, re-using it")
                else:
                    # logger.info("\n\n\n%s\n\n\n" % (sexcmd))
                    start_time = time.time()
                    try: 
                        ret = subprocess.Popen(sexcmd.split(), 
                                               stdout=subprocess.PIPE, 
                                               stderr=subprocess.PIPE)
                        (sex_stdout, sex_stderr) = ret.communicate()
                        #os.system(sexcmd)
                        if (ret.returncode != 0):
                            logger.warning("Sextractor might have a problem, check the log")
                            logger.debug("Stdout=\n"+sex_stdout)
                            logger.debug("Stderr=\n"+sex_stderr)
                    except OSError as e:
                        podi_logging.log_exception()
                        print >>sys.stderr, "Execution failed:", e
                    end_time = time.time()
                    logger.debug("SourceExtractor returned after %.3f seconds" % (end_time - start_time))
                    # os.system(sexcmd)
                    logger.debug("Sextractor complete")

                # Load the catalog
                src_catalog = numpy.loadtxt(catfile)
                logger.info("Found %d stars in frame" % (src_catalog.shape[0]))

                # Now eliminate all frames with magnitude 99
                good_photometry = src_catalog[:,SXcolumn['mag_aper_3.0']] < 0
                src_catalog = src_catalog[good_photometry]
                logger.info("%s stars have good photometry (no flags)" % (src_catalog.shape[0]))

                hdulist = pyfits.open(inputfile) #raw_frame)

                try:
                    filter_name = hdulist[0].header['FILTER']
                except:
                    filter_name = cmdline_arg_set_or_default('-filter', 'odi_r')

                exptime = 1.0

                try:
                    plottitle = "%(file)s\n%(object)s -- %(exptime).1f sec -- %(filter)s" % {
                        'file': inputfile,
                        'object': hdulist[0].header['OBJECT'],
                        'exptime': hdulist[0].header['EXPTIME'],
                        'filter': hdulist[0].header['FILTER'],
                    }
                except:
                    plottitle = "raw file"

                hdulist.close()

                pc =  photcalib(src_catalog, 
                                inputfile, 
                                filtername=filter_name, 
                                exptime=exptime,
                                diagplots=True,
                                plottitle=plottitle,
                                otalist=None,
                                options=options,
                                verbose=False,
                                eliminate_flags=True)
                zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = pc
                logger.info("Zeropoint: %.4f +/- %f (#stars: %d)\n" % (
                    zeropoint_median, zeropoint_std, odi_sdss_matched.shape[0])
                )

                numpy.savetxt(inputfile[:-5]+".matchedcat", odi_sdss_matched)

            except:
                podi_logging.log_exception()
                pass

        # Shutdown logging to shutdown cleanly
        podi_logging.shutdown_logging(options)

    elif (cmdline_arg_isset("-aucap")):

        import podi_collectcells
        options = podi_collectcells.read_options_from_commandline()

        print "Starting phot-calib..."
        for inputfile in get_clean_cmdline()[1:]:

            print "Working on",inputfile

            if (inputfile[-3:] == ".fz"):
                if (os.path.isfile(inputfile[-3:])):
                    print "given fz-compressed file, but uncompressed file also exists..."
                    inputfile = inputfile[:-3]
                else:
                    print "found fz-compressed file, unpacking..."
                    os.system("funpack -v "+inputfile)
                    inputfile = inputfile[:-3]

                print "continuing work on",inputfile
                if (not os.path.isfile(inputfile)):
                    print "file not found, something must have gone wrong with funpack"
                    continue
                          
            # Run SourceExtractor
            print "Running source-extractor"
            sex_config_file = "%s/.config/wcsfix.sex" % (options['exec_dir'])
            parameters_file = "%s/.config/wcsfix.sexparam" % (options['exec_dir'])
            catfile = "%s.cat" % (inputfile[:-5])
            sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s %s" % (
                        sitesetup.sextractor, sex_config_file, parameters_file, catfile, 
                        inputfile, sitesetup.sex_redirect)
            print sexcmd
            if (os.path.isfile(catfile) and not cmdline_arg_isset("-resex")):
                print "catalog exists, re-using it"
                if (not cmdline_arg_isset("-replot")):
                    continue

            else:
                if (options['verbose']): print sexcmd
                os.system(sexcmd)

            hdulist = pyfits.open(inputfile)
            filter_name = hdulist[0].header['FILTER']
            exptime = hdulist[0].header['EXPTIME']

            plottitle = "%s\nMAGZERO=%.3f/%.3f MAGZSIG=%.3f MAGZERR=%.3f" % (inputfile,
                                                                        hdulist[0].header['MAGZERO'],
                                                                        hdulist[0].header['MAGZERO']-2.5*math.log10(exptime),
                                                                        hdulist[0].header['MAGZSIG'],
                                                                        hdulist[0].header['MAGZERR'])
            hdulist.close()

            # Load the catalog
            src_catalog = numpy.loadtxt(catfile)
            print "Found",src_catalog.shape[0],"stars in frame"

            # Now eliminate all frames with magnitude 99
            good_photometry = src_catalog[:,SXcolumn['mag_aper_3.0']] < 0
            src_catalog = src_catalog[good_photometry]
            print src_catalog.shape[0],"stars with good photometry"


            outputfile = inputfile #[:-5]+".photcalib"
            pc =  photcalib(src_catalog, 
                            output_filename=outputfile, 
                            filtername=filter_name, 
                            exptime=exptime,
                            diagplots=True,
                            plottitle=plottitle,
                            otalist=None,
                            options=options,
                            verbose=False,
                            eliminate_flags=True)
            zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = pc
            print zeropoint_median
            print zeropoint_std

    elif (cmdline_arg_isset("-zpmapcombine")):

        from podi_commandline import *
        import podi_logging

        options = read_options_from_commandline(None)
        podi_logging.setup_logging(options)

        plotfile = get_clean_cmdline()[1]
        logger = logging.getLogger("ZPMap")

        import matplotlib, matplotlib.pyplot
        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)

        ota_pixelsize = 4100.

        #colors = np.r_[np.linspace(0.1, 1, 5), np.linspace(0.1, 1, 5)] 
        mymap = matplotlib.pyplot.get_cmap("spectral")

        # get the colors from the color map
        #my_colors = mymap(colors)
        # here you give floats as color to scatter and a color map
        #         # scatter "translates" this
        # axes[0].scatter(data[:, 0], data[:, 1], s=40,
        #                 c=colors, edgecolors='None',
        #                 cmap=mymap)
        # for n in range(5):
        #     # here you give a color to scatter
        #     axes[1].scatter(data[n, 0], data[n, 1], s=40,
        #                     color=my_colors[n], edgecolors='None',
        #                     label="point %d" %(n))
        # # by default legend would show multiple scatterpoints (as you would normally
        # # plot multiple points with scatter)
        # # I reduce the number to one here
        # plt.legend(scatterpoints=1)
        # plt.tight_layout()
        # plt.show()

        zp_range = [-0.2, 0.2]

        ota_xmin, ota_xmax, ota_ymin, ota_ymax = 2,4,2,4

        all_x = numpy.array([])
        all_y = numpy.array([])
        all_zp = numpy.array([])
        all_err = numpy.array([])

        for fitsfile in get_clean_cmdline()[2:]:
            logger.info("Reading %s ..." % (fitsfile))
            try:
                hdulist = pyfits.open(fitsfile)
                tbhdu = hdulist['CAT.PHOTCALIB']
            except:
                continue

            filtername = hdulist[0].header['FILTER']
            sdss_eq = sdss_equivalents[filtername]
            # print sdss_eq

            sdss_mag_col = 'SDSS_MAG_%s' % (sdss_eq)
            sdss_mag = tbhdu.data.field(sdss_mag_col)
            
            odi_mag_col = 'ODI_MAG_D%02d' % (10*hdulist[0].header['MAG0SIZE'])
            odi_mag = tbhdu.data.field(odi_mag_col)

            odi_err_col = 'ODI_ERR_D%02d' % (10*hdulist[0].header['MAG0SIZE'])
            odi_err = tbhdu.data.field(odi_err_col)

            ota = tbhdu.data.field('ODI_OTA')
            ota_x = numpy.floor(ota/10.)
            ota_y = numpy.fmod(ota, 10.)
            x = tbhdu.data.field('ODI_X') + ota_x * ota_pixelsize
            y = tbhdu.data.field('ODI_Y') + ota_y * ota_pixelsize

            zp = sdss_mag - odi_mag# - hdulist[0].header['PHOTZP']
            zp_median = numpy.median(zp)
            zp -= zp_median

            # print zp[:10]

            good_points = (odi_err < 0.05) & \
                          (ota_x >= ota_xmin) & (ota_x <= ota_xmax) & \
                          (ota_y >= ota_ymin) & (ota_y <= ota_ymax)

            all_x = numpy.append(all_x, x[good_points])
            all_y = numpy.append(all_y, y[good_points])
            all_zp = numpy.append(all_zp, zp[good_points])
            all_err = numpy.append(all_err, odi_err[good_points])

            # matplotlib.pyplot.scatter(x,y, c=zp, 
            #                           cmap=mymap,
            #                           vmin=zp_range[0], vmax=zp_range[1], 
            #                           edgecolor='none',
            #                           s=9,
            #                       )

            # ax.scatter(x,y, c=zp, 
            #                           cmap=mymap,
            #                           vmin=zp_range[0], vmax=zp_range[1], 
            #                           edgecolor='none',
            #                           s=9,
            #                       )

        import scipy.interpolate 

        logger.info("Plotting...")

        minx, maxx = ota_xmin*ota_pixelsize, (ota_xmax+1)*ota_pixelsize
        miny, maxy = ota_ymin*ota_pixelsize, (ota_ymax+1)*ota_pixelsize

        if (cmdline_arg_isset('-interpolate')):

            method = get_cmdline_arg_set_or_default("-interpolate", "linear")

            grid_x, grid_y = numpy.mgrid[minx:maxx:100j, miny:maxy:100j]

            xy = numpy.zeros(shape=(all_x.shape[0], 2))
            xy[:,0] = all_x
            xy[:,1] = all_y

            grid_z0 = scipy.interpolate.griddata(xy, all_zp, (grid_x, grid_y), method=method)

            matplotlib.pyplot.imshow(grid_z0.T, extent=(minx,maxx,miny,maxy), 
                                     origin='lower',
                                     cmap=mymap,
                                     vmin=zp_range[0], vmax=zp_range[1],)


        if (cmdline_arg_isset("-spline")):

            err = numpy.hypot(all_err, 0.01)
            splinemode = cmdline_arg_set_or_default("-spline", 'nw,5')
            items = splinemode.split(',')

            weight = numpy.ones(err.shape)
            weight_txt = ''
            if items[0] == 'w':
                weight = 0.01/err
                weight_txt = "weighted "

            logger.info("Fitting %sspline surface to residuals" % (weight_txt))


            if (len(items) <= 2):
                nodes = int(items[1]) if len(items) == 2 else 5
                xcoord = numpy.linspace(minx, maxx, nodes)
                ycoord = numpy.linspace(miny, maxy, nodes)
            elif (len(items) >= 3):
                xcoord = numpy.linspace(minx, maxx, int(items[1]))
                ycoord = numpy.linspace(miny, maxy, int(items[2]))
                
            # print all_err[:10], err[:10]
            grid_x, grid_y = numpy.mgrid[minx:maxx:100j, miny:maxy:100j]
            #grid_x, grid_y = numpy.mgrid[minx:maxx:24j, miny:maxy:24j]
            # bispline = scipy.interpolate.bisplrep(x=all_x, y=all_y, z=all_zp,
            #                                       w=weight,
            #                                       # xb=minx, xe=maxx, yb=miny, ye=maxy,
            #                                       kx=3, ky=3,
            #                                       s=0.2, #all_x.shape[0],
            #                                   )

            # znew = scipy.interpolate.bisplev(grid_x[:,0],grid_y[0,:],bispline)

            # Create the knots (10 knots in each direction, making 100 total
            #xcoord = numpy.linspace(numpy.min(all_x), numpy.max(all_x), 5)
            #ycoord = numpy.linspace(numpy.min(all_y), numpy.max(all_y), 5)
            spl = scipy.interpolate.LSQBivariateSpline(x=all_x, y=all_y, z=all_zp, 
                                                       tx=xcoord, ty=ycoord,
                                                       w=weight, 
                                                       kx=3, ky=3,
                                                       #bbox=[minx,maxx,miny,maxy]
            )
            znew = spl(grid_x[:,0],grid_y[0,:])

            matplotlib.pyplot.imshow(znew.T, extent=(minx,maxx,miny,maxy), 
                                     origin='lower',
                                     cmap=mymap,
                                     vmin=zp_range[0], vmax=zp_range[1],
                                     interpolation='none')


            # labels_x, labels_y = numpy.mgrid[minx:maxx:10j, miny:maxy:10j]
            # zlabels = scipy.interpolate.bisplev(labels_x[:,0],labels_y[0,:],bispline)
            # print zlabels.shape
            # for x in range(zlabels.shape[0]):
            #     for y in range(zlabels.shape[1]):
            #         print x,y,labels_x[x,0], labels_y[y,0],"%.3f" % zlabels[x,y]
            #         ax.text(labels_x[x], labels_y[y], 
            #                                "%.3f" % zlabels[x,y], 
            #                                ha='center', va='center', color='black')

            #grid_x, grid_y = numpy.mgrid[minx:maxx:24j, miny:maxy:24j]
 
        matplotlib.pyplot.scatter(all_x,all_y, c=all_zp, 
                                  cmap=mymap,
                                  vmin=zp_range[0], vmax=zp_range[1], 
                                  edgecolor='none',
                                  s=9,
                                  )

        ax.scatter(all_x,all_y, c='black', marker='.', 
                                  edgecolor='none', s=3,
                                  )


        #
        # Draw the limits of the OTAs
        #
        for i in range(ota_ymin, ota_ymax+1):
            ax.axhline(y=i*ota_pixelsize, linewidth=1, color='white', ls='-')
            #ax.axhline(y=i*ota_pixelsize, linewidth=1, color='black', ls=':')
        for i in range(ota_xmin, ota_xmax+1):
            ax.axvline(x=i*ota_pixelsize, linewidth=1, color='white', ls='-')
            # ax.axvline(x=i*ota_pixelsize, linewidth=1, color='black', ls=':')

        # print ax.get_xticklabels()
        # for x in range(len(ax.get_xticklabels())): 
        #     print ax.get_xticklabels()[x]
        #     del ax.get_xticklabels()[x]
        # ax.set_xticklabels([])
        # ax.set_xticks([])

        matplotlib.pyplot.xlim((minx, maxx))
        matplotlib.pyplot.ylim((miny, maxy))
#        matplotlib.pyplot.ylim((0, 7*ota_pixelsize))

        matplotlib.pyplot.xlabel("OTA")
        matplotlib.pyplot.ylabel("OTA")
        #label_text = ["%d" % ota for ota in range(minx,7)]
        #label_pos = [(x+0.5)*ota_pixelsize for x in range(0,7)]
        matplotlib.pyplot.xticks([(ota+0.5)*ota_pixelsize for ota in range(ota_xmin,ota_xmax+1)],
                                 ["%d" % ota for ota in range(ota_xmin, ota_xmax+1)])
        matplotlib.pyplot.yticks([(ota+0.5)*ota_pixelsize for ota in range(ota_ymin, ota_ymax+1)],
                                 ["%d" % ota for ota in range(ota_ymin, ota_ymax+1)])
#        matplotlib.pyplot.yticks(label_pos, label_text)

        #ax.set_xticks(label_pos, label_text)
        #ax.set_yticks(label_pos, label_text)
        #ax.set_xlim((0, 7*ota_pixelsize))
        #ax.set_ylim((0, 7*ota_pixelsize))

        colorbar = matplotlib.pyplot.colorbar(cmap=mymap)
        colorbar.set_label("phot. zeropoint deviation from median")

        # fig.colorbar(ax)
        #matplotlib.pyplot.show()
        #fig.show()
        # fig.savefig(plotfile)
        matplotlib.pyplot.savefig(plotfile)
        logger.info("Results written to %s!\n" % (plotfile))

        podi_logging.shutdown_logging(options)

    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        calibdir = get_clean_cmdline()[3]

        photcalib(fitsfile, output_filename, calibdir)
