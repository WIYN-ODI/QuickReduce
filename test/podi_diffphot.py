#!/usr/bin/env python

import os, sys
d,_=os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append("%s/../"%d)

# import podi_asteroids
import podi_swarpstack
from podi_commandline import *
from podi_definitions import *
import podi_logging
import logging
import astropy.io.fits
import astropy.io.votable
import astropy.time
import math
import numpy
import ephem
import time
import scipy.stats
import scipy.spatial
import bottleneck
import matplotlib.pyplot

import podi_sitesetup as sitesetup


def make_catalogs(inputlist):


    for fitsfile in inputlist:
        logger.info("Creating source catalog for %s" % (fitsfile))
        catfile = "%s.cat" % (fitsfile[:-5])

        if (os.path.isfile(catfile)):
            # Don't do anything if the catalog already exists
            continue

        sex_config_file = "%s/.config/diffphot.sex" % (sitesetup.exec_dir)
        parameters_file = "%s/.config/diffphot.sexparam" % (sitesetup.exec_dir)
        sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s" % (
            sitesetup.sextractor, sex_config_file, parameters_file, catfile, 
            fitsfile)
        
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

    


#
#
# To-do list
#
# - make sure every star is only part of a single tracklet
#
#




def differential_photometry(inputlist, source_coords):

    logger = logging.getLogger("DiffPhot")

    # interpret source_coords
    sc_items = source_coords.split(",")
    src_ra = float(sc_items[0])
    src_dec = float(sc_items[1])

    src_dra, src_ddec, src_dmjd = 0., 0., 0.

    if (len(sc_items) == 5):
        src_dra = float(sc_items[2])
        src_ddec = float(sc_items[3])
        if (os.path.isfile(sc_items[4])):
            hdu = astropy.io.fits.open(sc_items[4])
            src_mjd = hdu[0].header['MJD-OBS']
        else:
            src_mjd = float(sc_items[4])

    print src_ra, src_dec, src_dra, src_ddec, src_mjd

    # Internally, all coordinates are corrected for cos(dec)
    cos_declination = math.cos(math.radians(src_dec))
    src_ra *= cos_declination

    make_catalogs(inputlist)
    
    catalog_filelist = []
    catalog = []
    mjd_hours = []
    source_id = []
    filters = []

    phot_ZP = 25.

    # set apertures to use
    aperture_size = '4.0'
    aperture_column = SXcolumn["mag_aper_%s" % aperture_size]
    aperture_error_column = SXcolumn["mag_err_%s" % aperture_size]

    # Now open all files, get time-stamps etc
    for fitsfile in inputlist:

        catalog_filename = fitsfile[:-5]+".cat"
        dumpfile = fitsfile[:-5]+".catdump"
        if (os.path.isfile(dumpfile)):
            cat_data = numpy.load(dumpfile)
        else:
            cat_data = numpy.loadtxt(catalog_filename)
            cat_data.dump(dumpfile)

        # Correct for the absolute zeropoint
        magzero = 25.
        hdu = astropy.io.fits.open(fitsfile)
        magzero = hdu[0].header['MAGZERO']
        magzero = hdu[0].header['PHOTZP_X']

        cat_data[:, SXcolumn['mag_aper_2.0']:SXcolumn['mag_aper_12.0']+1] += (magzero - 25.0)

        # Account for cos(dec) distortion
        cat_data[:, SXcolumn['ra']] *= cos_declination

        hdulist = astropy.io.fits.open(fitsfile)
        obs_mjd = hdulist[0].header['MJD-OBS']
        exptime = hdulist[0].header['EXPTIME']
        mjd_middle = (obs_mjd * 24.) + (exptime / 3600.)
        filter_name = hdulist[0].header['FILTER']

        # Do some preparing of the catalog
        good_photometry = cat_data[:, aperture_column] < 0
        cat_data = cat_data[good_photometry]

        
        # Apply photometric zeropoint to all magnitudes
        for key in SXcolumn_names:
            if (key.startswith("mag_aper")):
                cat_data[:, SXcolumn[key]] += phot_ZP

        catalog_filelist.append(catalog_filename)
        catalog.append(cat_data)
        mjd_hours.append(mjd_middle)
        source_id.append(numpy.arange(cat_data.shape[0]))
        filters.append(filter_name)

        print catalog_filename, cat_data.shape


    # Now find sources that have good photometry in all frames that 
    # we can use as reference stars
    ref_stars = catalog[0]
    small_errors = (ref_stars[:, aperture_error_column] < 0.05) & \
                   (ref_stars[:, aperture_column] > 16)
    no_flags = ref_stars[:, SXcolumn['flags']] == 0
    ref_stars = ref_stars[(small_errors & no_flags)]

    # Sort by magnitude and pick the 10 brightest stars
    si = numpy.argsort(ref_stars[:, aperture_column])
    ref_sorted = ref_stars[si]
    print ref_sorted[:, aperture_column]

    n_ref = 25
    phot_refs = ref_sorted[0:n_ref,0:2]
    print phot_refs[:,0:2]

    # Now go through all catalogs and find the source and the reference 
    # magnitudes

    output_file = open('output.test', 'w')
    match_radius = 2./3600. # 2 arcsec

    all_cats = []
    target_source = []
    target_mjd = []

    #
    # Create source catalogs for the source and n reference sources from each catalog
    #
    
    for i_cat in range(len(catalog)):
        print "\n"*5

        src_cat = catalog[i_cat]

        # turn the catalog coordinates into a KDTree
        cat_tree = scipy.spatial.cKDTree(catalog[i_cat][:,0:2])
        # print catalog[i_cat][:,0:2]

        # Search for the actual source
        # add here: fudge with coordinates to account for motion
        src_radec = numpy.array([[src_ra, src_dec]]) + (mjd_hours[i_cat] - src_mjd*24.) * numpy.array([src_dra, src_ddec]) / 3600.
        # print src_radec * numpy.array([1./cos_declination, 1.])
        
        # print src_radec.shape, phot_refs.shape

        #########################################################################
        #
        # Determine the magnitude of the source
        #
        #########################################################################

        # Search for the source we are interested in
        src_counterpart = cat_tree.query_ball_point(src_radec, r=match_radius, p=2)
        src_data = None
        if (len(src_counterpart) <= 0):
            # We couldn't find any source counterpart
            continue
        elif (len(src_counterpart) > 1):
            # Found multiple counterparts
            # add some brain to pick the most likely one
            continue
        else:
            # 
            try:
                src_data = (catalog[i_cat][src_counterpart[0]])[0]
            except:
                podi_logging.log_exception()
                continue
            print 'src-data=',src_data.shape

        #########################################################################
        #
        # Determine the aperture correction from the source magnitude 
        # to the 12'' aperture we consider for "total intensity"
        #
        #########################################################################

        # Find the aperture magnitude that is valid
        numpy.savetxt(sys.stdout, [src_data[SXcolumn['mag_aper_2.0']:SXcolumn['mag_err_12.0']+1]], " %8.4f")
        good_source_aperture = '4.0'
        src_aper_mag = "mag_aper_%s" % (good_source_aperture)
        src_aper_err = "mag_err_%s" % (good_source_aperture)

        # By default, use the 4.0'' aperture, but aperture-correct it to 12.0 arcsec
        # Use the weighted difference between the 4 and 12'' magnitude 
        valid_photometry = (catalog[i_cat][:,SXcolumn[src_aper_err]] < 0.1) & \
                           (catalog[i_cat][:,SXcolumn['mag_err_12.0']] < 0.1) & \
                           (catalog[i_cat][:, SXcolumn['flags']] == 0)
        aperture_diffs = (catalog[i_cat][:,SXcolumn['mag_aper_12.0']] 
                          - catalog[i_cat][:, SXcolumn[src_aper_mag]])[valid_photometry]
        aperture_diff_filtered = three_sigma_clip(aperture_diffs)
        aperture_correction = numpy.median(aperture_diff_filtered)
        aperture_correction_std = numpy.std(aperture_diff_filtered)

        numpy.savetxt("apcorr%d.raw" % (i_cat), aperture_diffs)
        numpy.savetxt("apcorr%d.filter" % (i_cat), aperture_diff_filtered)
        print aperture_correction, aperture_correction_std


        #########################################################################
        #
        # Now determine all magnitudes for all the reference sources
        #
        #########################################################################

        # all_sources = numpy.append(src_radec, phot_refs, axis=0)
        # print all_sources.shape, match_radius
        src_tree = scipy.spatial.cKDTree(phot_refs) 
        counterparts = src_tree.query_ball_tree(cat_tree, r=match_radius, p=2)

        phot_data = numpy.ones((phot_refs.shape[0], catalog[i_cat].shape[1])) 
        print phot_data.shape
        phot_data[:] = numpy.NaN

        print >>output_file, mjd_hours[i_cat], mjd_hours[i_cat]/24., \
            src_radec[0,0]/cos_declination, src_radec[0,1], \
            inputlist[i_cat],
        for i_src in range(len(counterparts)):
            src = counterparts[i_src]
            if (len(src) == 1):
                # Found a unique counterpart
                # remember the magnitude and error
                phot_data[i_src,:] = src_cat[src[0], :]
            if (len(src) == 1):
                print >>output_file, \
                    src_cat[src[0], aperture_column],\
                    src_cat[src[0], aperture_error_column],\
                    src_cat[src[0], SXcolumn['ra']]/cos_declination, \
                    src_cat[src[0], SXcolumn['dec']],
            else:
                print >>output_file, -99, -99, "x", "x",
            print >>output_file, len(src),

        print >>output_file

        #
        # Use the mag_auto data column to hold the aperture-corrected magnitude
        #
        phot_data[:, SXcolumn['mag_auto']] = phot_data[:, SXcolumn[src_aper_mag]] + aperture_correction
        phot_data[:, SXcolumn['mag_err_auto']] = phot_data[:, SXcolumn[src_aper_err]]
        numpy.savetxt("photcat%d" % (i_cat), phot_data)
        #
        # also apply aperture correction to target source
        #
        src_data[SXcolumn['mag_auto']] = src_data[SXcolumn[src_aper_mag]] + aperture_correction
        src_data[SXcolumn['mag_err_auto']] = src_data[SXcolumn[src_aper_err]]


        #
        # Combine all catalogs for later
        #
        all_cats.append(phot_data)
        target_source.append(src_data)
        print "Added target source", len(target_source)
        target_mjd.append(mjd_hours[i_cat])


    print "\n"*10

    # print target_source
    # print "shapes:",[i.shape for i in target_source]

    #
    # Convert catalogs into proper numpy catalog
    #
    all_cats = numpy.array(all_cats)
    target_source = numpy.array(target_source)
    target_mjd = numpy.array(target_mjd)

    print target_source.shape
    print all_cats.shape

    #############################################################################
    #
    # Now compute the median and filtered magnitude of each reference star
    #
    #############################################################################


    # Iteratively eliminate outliers
    for i in range(3):
        median_mags = bottleneck.nanmedian(all_cats[:,:, SXcolumn['mag_auto']], axis=0)
        std_mags = bottleneck.nanstd(all_cats[:,:, SXcolumn['mag_auto']], axis=0)
        numpy.savetxt("photcat_median_%d" % i, median_mags)
        numpy.savetxt("photcat_std_%d" % i, std_mags)
        outlier = (all_cats[:,:, SXcolumn['mag_auto']] > median_mags + 0.25) | \
                  (all_cats[:,:, SXcolumn['mag_auto']] < median_mags - 0.25)
        #print outlier

        if (numpy.sum(outlier) <= 0):
            break

        all_cats[outlier, :] = numpy.NaN

    
    #############################################################################
    #
    # Compute the differential photometry correction
    # 
    #############################################################################

    # Now we have the median magnitudes of each reference star.
    # Additionally, all outliers have already been flagged as NaNs
    print median_mags

    # Compute the deviation of each source from its median magnitude
    mag_deviation = all_cats[:,:, SXcolumn['mag_auto']] - median_mags
    print mag_deviation

    # Compute the median deviation across all stars for each of the source catalogs
    diffphot_correction = numpy.median(mag_deviation, axis=1)
    print diffphot_correction
    print diffphot_correction.shape

    numpy.save("correction", diffphot_correction)
    numpy.save("all_cats", all_cats)

    # As verification, apply the correction to each of the reference star catalogs
    print all_cats.shape
    all_cats[:,:, SXcolumn['mag_auto']] -= diffphot_correction.reshape((diffphot_correction.shape[0],1))
    for i in range(all_cats.shape[0]):
        numpy.savetxt("photcat%d.corr2" % (i), all_cats[i,:,SXcolumn['mag_auto']])

    #############################################################################
    #
    # Now apply the differential photometry correction to the actual source photometry
    #
    #############################################################################

    numpy.savetxt("mjd", target_mjd)
    numpy.savetxt("target_source.raw", target_source)
    print target_source.shape
    print target_source
    target_source[:, SXcolumn['mag_auto']] -= diffphot_correction
    numpy.savetxt("target_source.fixed", target_source)

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    time = target_mjd  /24.
    ax.errorbar(time, target_source[:, SXcolumn['mag_auto']], 
                yerr=target_source[:, SXcolumn['mag_err_auto']], 
                fmt='o', c='blue')
    ax.scatter(time, target_source[:, SXcolumn['mag_auto']],
               c='blue')
    ax.set_xlabel("MJD")
    ax.set_ylabel("apparent magnitude [12'' aperture']")
    # ax.set_ylim((-0.5, 0.5))
    matplotlib.pyplot.show()
    fig.show()
    fig.savefig('test.png')



    return



if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)
    logger = logging.getLogger("Asteroidmetry: Main")

    source_data = get_clean_cmdline()[1]
    catalog_dir = None #get_clean_cmdline()[2]
    inputlist = get_clean_cmdline()[2:]

    print "source data", source_data
    print "filelist",inputlist

    data = differential_photometry(inputlist, source_coords=source_data)

    logger.info("Done, shutting down")
    podi_logging.shutdown_logging(options)


    sys.exit(0)
