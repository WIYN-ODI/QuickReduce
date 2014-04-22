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

        cat_data[:, aperture_column] += (magzero - 25.0)

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

    #
    # Create source catalogs for the source and n reference sources from each catalog
    #

    for i_cat in range(len(catalog)):
        
        src_cat = catalog[i_cat]

        # turn the catalog coordinates into a KDTree
        cat_tree = scipy.spatial.cKDTree(catalog[i_cat][:,0:2])
        # print catalog[i_cat][:,0:2]

        # Search for the actual source
        # add here: fudge with coordinates to account for motion
        src_radec = numpy.array([[src_ra, src_dec]]) + (mjd_hours[i_cat] - src_mjd*24.) * numpy.array([src_dra, src_ddec]) / 3600.
        print src_radec * numpy.array([1./cos_declination, 1.])
        
        # print src_radec.shape, phot_refs.shape

        all_sources = numpy.append(src_radec, phot_refs, axis=0)
        # print all_sources.shape, match_radius
        src_tree = scipy.spatial.cKDTree(all_sources) #catalog[i_cat][:,0:2])

        # Search for the source we are interested in

        counterparts = src_tree.query_ball_tree(cat_tree, r=match_radius, p=2)
        # print counterparts

        phot_data = numpy.ones((all_sources.shape[0], 2)) 
        phot_data[:] = numpy.NaN

        print >>output_file, mjd_hours[i_cat], mjd_hours[i_cat]/24., \
            src_radec[0,0]/cos_declination, src_radec[0,1], \
            inputlist[i_cat],
        for i_src in range(len(counterparts)):
            src = counterparts[i_src]
            if (len(src) == 1):
                # Found a unique counterpart
                # remember the magnitude and error
                phot_data[i_src,0] = src_cat[src[0], aperture_column]
                phot_data[i_src,1] = src_cat[src[0], aperture_error_column]
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

        all_cats.append(phot_data)

    all_cats = numpy.array(all_cats)
    #print all_cats
    #print all_cats.shape
    # numpy.savetxt('all_cats.dat', all_cats, "%6.3f")


    #
    # Compute the average magnitude of each star to isolate just the variance
    #
    weights = 1./all_cats[:,:,1]
    weighted = all_cats[:,:,0] / all_cats[:,:,1]
    avg_mags = bottleneck.nansum(weighted, axis=0) / bottleneck.nansum(weights, axis=0)
    # avg_mags = numpy.average(all_cats[:,:,0], weights=1./all_cats[:,:,1], axis=1)
#    print avg_mags
#    print all_cats.shape

    # 
    # Add some filtering to eliminate outliers in the reference star magnitudes
    # These can be due to edge-effects, cosmics, or ...
    #
    median = numpy.zeros(all_cats.shape[1])
    std = numpy.ones(all_cats.shape[1]) * 50

    # Extract only the reference stars
    ref_cat = all_cats[:,1:,:]
    print all_cats.shape, ref_cat.shape

    for loop in range(3):
        median_mags = numpy.median(ref_cat[:,:,0], axis=0)
        print median_mags
        sigmas = numpy.array(scipy.stats.scoreatpercentile(ref_cat[:,:,0], [16,50,84], axis=0))

        upper_limit = sigmas[1,:] + 3*(sigmas[2,:]-sigmas[1,:])
        lower_limit = sigmas[1,:] - 3*(sigmas[2,:]-sigmas[1,:])

        outlier = (ref_cat[:,:,0] > upper_limit) | (ref_cat[:,:,0] < lower_limit)
        ref_cat[outlier] = numpy.NaN

    # Re-insert the reference stars in the full catalog
    all_cats[:,1:,:] = ref_cat

    # subtract the average magnitude of each star to only derive the 
    # deviation from the average
    # don't subtract anything from source 0 as this is what we are interested in
    avg_mags[0] = 0.

    corr_cats = all_cats.copy()
    corr_cats[:,:,0] -= avg_mags

    print >>output_file, "\n"*10
    numpy.savetxt(output_file, all_cats[:,:,0])
    print >>output_file, "\n"*10

    #
    # compute the average deviation of all reference stars
    #
    weights = 1/corr_cats[:,1:,1]
    weighted = corr_cats[:,1:, 0] * weights
#    print weighted
    correction = bottleneck.nansum(weighted, axis=1) / bottleneck.nansum(weights, axis=1)
    d_corr = bottleneck.nanstd(corr_cats[:,1:, 0], axis=1)
    print correction.shape

    raw_mags = corr_cats[:,0,:]
    corr_mags = raw_mags[:,0] - correction

#    print raw_mags
#    print corr_mags

    print >>output_file, "\n"*10
    numpy.savetxt(output_file, avg_mags)

    print >>output_file, "\n"*10
    for i in range(len(mjd_hours)):
        print >>output_file, i, raw_mags[i,0], raw_mags[i,1], corr_mags[i], correction[i], d_corr[i]

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    print corr_cats.shape
    for i in range(corr_cats.shape[1]):
        ax.plot((numpy.array(mjd_hours)-mjd_hours[0])*60, all_cats[:,i,0])
    # ax.set_ylim((-0.5, 0.5))
    matplotlib.pyplot.show()
    fig.show()
    fig.savefig('test.png')

    for i in range(len(mjd_hours)):
        print >>output_file, "\n"*10
        numpy.savetxt(output_file, all_cats[i])
        
#     print all_cats


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
