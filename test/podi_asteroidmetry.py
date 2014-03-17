#!/usr/bin/env python

import os, sys
d,_=os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append("%s/../"%d)

import podi_asteroids
import podi_swarpstack
from podi_commandline import *
from podi_definitions import *
import podi_logging
import logging
import astropy.io.votable
import math
import numpy
import ephem
import time
import scipy.stats
import scipy.spatial

import podi_sitesetup as sitesetup


def make_catalogs(inputlist):


    for fitsfile in inputlist:
        logger.info("Creating source catalog for %s" % (fitsfile))
        catfile = "%s.cat" % (fitsfile[:-5])

        if (os.path.isfile(catfile)):
            # Don't do anything if the catalog already exists
            continue

        sex_config_file = "%s/.config/moving.sex" % (sitesetup.exec_dir)
        parameters_file = "%s/.config/moving.sexparam" % (sitesetup.exec_dir)
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

    

if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    logger = logging.getLogger("Astro(id)metry")

    inputlist = get_clean_cmdline()[1:]

    make_catalogs(inputlist)
    
    catalog_filelist = []
    catalog = []
    mjd_hours = []
    source_id = []

    # Now open all files, get time-stams etc
    for fitsfile in inputlist:

        catalog_filename = fitsfile[:-5]+".cat"
        cat_data = numpy.loadtxt(catalog_filename)
        
        hdulist = astropy.io.fits.open(fitsfile)
        obs_mjd = hdulist[0].header['MJD-OBS']
        exptime = hdulist[0].header['EXPTIME']
        mjd_middle = (obs_mjd * 24.) + (exptime / 3600.)

        # Do some preparing of the catalog
        good_photometry = cat_data[:, SXcolumn['mag_aper_2.0']] < 0

#        good_photometry = (cat_data[:, SXcolumn['x']] > 6100) & (cat_data[:, SXcolumn['x']] < 6300)
#                          (cat_data[:, SXcolumn['y']] > 5100) & (cat_data[:, SXcolumn['y']] < 5350)
        print numpy.sum(good_photometry)
        # if (numpy.sum(good_photometry) <= 0):
        cat_data = cat_data[good_photometry]

        catalog_filelist.append(catalog_filename)
        catalog.append(cat_data)
        mjd_hours.append(mjd_middle)
        source_id.append(numpy.arange(cat_data.shape[0]))

        print catalog_filename, cat_data.shape


    max_rate = 100 # arcsec / hour
    min_distance = 2. / 3600. # arcsec
    motion_rates = numpy.zeros((0,6))

    for ref_cat in range(len(catalog)-1):
        #motion_rates = numpy.zeros((0,6))
        for second_cat in range(ref_cat+1, len(catalog)):
            
            _motion_rates = numpy.zeros((0,6))
            print catalog_filelist[ref_cat], catalog_filelist[second_cat]

            d_hours = mjd_hours[second_cat] - mjd_hours[ref_cat]
            diff_to_rate = 3600 / d_hours

            for obj1 in range(catalog[ref_cat].shape[0]):
                # Take one source in first catalog; 
                # compute difference and rate from each source in the second 
                # catalog to the source in the first catalog.

                cos_declination = math.cos(math.radians(catalog[ref_cat][obj1, SXcolumn['dec']]))
                
                d_ra = catalog[second_cat][:,SXcolumn['ra']] - catalog[ref_cat][obj1,SXcolumn['ra']]
                d_ra *= cos_declination

                d_dec = catalog[second_cat][:,SXcolumn['dec']] - catalog[ref_cat][obj1,SXcolumn['dec']]

                d_total = numpy.hypot(d_ra, d_dec) 
                
                possible_match = (d_total < math.fabs(max_rate/diff_to_rate)) & (d_total > min_distance)
                # require a minimum distance of 2 arcsec
                
                matches = numpy.ones(shape=(numpy.sum(possible_match),6))
                matches[:,0] = ref_cat
                matches[:,1] = second_cat
                matches[:,2] = obj1
                matches[:,3] = source_id[second_cat][possible_match]
                matches[:,4] = d_ra[possible_match] * diff_to_rate
                matches[:,5] = d_dec[possible_match] * diff_to_rate

                _motion_rates = numpy.append(_motion_rates, matches, axis=0)

            motion_rates = numpy.append(motion_rates, _motion_rates, axis=0)

        #numpy.savetxt("motion_rates_%d" % (ref_cat), motion_rates)

    numpy.savetxt("motion_rates", motion_rates)

    
    logger.debug("Finding density")
    # Now that we have all potential rates, find density enhancements
    src_tree = scipy.spatial.cKDTree(motion_rates[:,4:6])
    match_tree = scipy.spatial.cKDTree(motion_rates[:,4:6])

    # Now select all neighbors within the search radius
    radius = 1.
    matches = match_tree.query_ball_tree(src_tree, radius, p=2)

    counts = numpy.zeros(shape=(motion_rates.shape[0], motion_rates.shape[1]+1))
    counts[:, :motion_rates.shape[1]] = motion_rates

    for i in range(len(matches)):
        counts[i,6] = len(matches[i])

    numpy.savetxt("motion_rates2", counts)

    median_count = numpy.median(counts[:,6])
    std_count = numpy.std(counts[:,6])
    print median_count, std_count

    logger.info("Filtering only significant matches")
    filtered, mask = three_sigma_clip(counts[:,6], return_mask=True)

    significant = counts[counts[:,6] > median_count+std_count]
    numpy.savetxt("motion_rates3", significant)

    #
    # Now we have some motion groups, look fro groups of at least 3 stars
    #
    candidate_tree = scipy.spatial.cKDTree(significant[:,4:6])
    full_tree = scipy.spatial.cKDTree(motion_rates[:,4:6])

    results = open("results", "w")
    chains = candidate_tree.query_ball_tree(full_tree, radius, p=2)
    for i in range(significant.shape[0]):
        for j in range(len(chains[i])):
            idx = chains[i][j]
            c1 = int(motion_rates[idx, 0])
            c2 = int(motion_rates[idx, 1])
            o1 = int(motion_rates[idx, 2])
            o2 = int(motion_rates[idx, 3])
            
            print i,j,idx, "   ", c1, c2, o1, o2

            ra1 = catalog[c1][o1, SXcolumn['ra']]
            dec1 = catalog[c1][o1, SXcolumn['dec']]
            ra2 = catalog[c2][o2, SXcolumn['ra']]
            dec2 = catalog[c2][o2, SXcolumn['dec']]
           
            mjd1 = mjd_hours[c1]
            mjd2 = mjd_hours[c2]
            print >>results, ra1, dec1, ra2, dec2, mjd1, mjd2
        print >>results, "\n\n\n\n\n"

    results.close()

            
    podi_logging.shutdown_logging(options)


