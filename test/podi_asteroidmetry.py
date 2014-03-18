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
#        print numpy.sum(good_photometry)
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

    # Now select all neighbors within the search radius
    radius = 1.

    results = open("results", "w")
    moving_object_id = 0

    for ref_cat in range(len(catalog)-1):
        logger.debug("Checking out catalog %d" % (ref_cat))

        for obj1 in range(catalog[ref_cat].shape[0]):
            # Take one source in first catalog; 
            # compute difference and rate from each source in the second 
            # catalog to the source in the first catalog.

            motion_rates = numpy.zeros((0,7))
            cos_declination = math.cos(math.radians(catalog[ref_cat][obj1, SXcolumn['dec']]))
            for second_cat in range(ref_cat+1, len(catalog)):
            
                # print catalog_filelist[ref_cat], catalog_filelist[second_cat]

                d_hours = mjd_hours[second_cat] - mjd_hours[ref_cat]
                diff_to_rate = 3600 / d_hours

                d_ra = catalog[second_cat][:,SXcolumn['ra']] - catalog[ref_cat][obj1,SXcolumn['ra']]
                d_ra *= cos_declination
                d_dec = catalog[second_cat][:,SXcolumn['dec']] - catalog[ref_cat][obj1,SXcolumn['dec']]
                d_total = numpy.hypot(d_ra, d_dec) 
                
                possible_match = (d_total < math.fabs(max_rate/diff_to_rate)) & (d_total > min_distance)
                # require a minimum distance of 2 arcsec
                
                matches = numpy.ones(shape=(numpy.sum(possible_match),7))
                matches[:,0] = ref_cat
                matches[:,1] = second_cat
                matches[:,2] = obj1
                matches[:,3] = source_id[second_cat][possible_match]
                matches[:,4] = d_ra[possible_match] * diff_to_rate
                matches[:,5] = d_dec[possible_match] * diff_to_rate

                motion_rates = numpy.append(motion_rates, matches, axis=0)


            # Now we have a whole set of potential counterparts for one star in the source catalog
            # Try to find significant matches
            if (motion_rates.shape[0] < 5):
                # not enough nearby stars found
                continue
            
            src_tree = scipy.spatial.cKDTree(motion_rates[:,4:6])
            match_tree = scipy.spatial.cKDTree(motion_rates[:,4:6])
            matches = match_tree.query_ball_tree(src_tree, radius, p=2)

            for i in range(len(matches)):
                motion_rates[i,6] = len(matches[i])

            filename = "motion_rates_%02d_%04d.dump" % (ref_cat, obj1)
            numpy.savetxt(filename, motion_rates)

            # numpy.savetxt(sys.stdout, motion_rates)

            #
            # Now we know how often each star appears
            #

            # check often the most frequent rate appears
            ratecount_max = numpy.max(motion_rates[:,6])

            if (ratecount_max > 4):
                # We have 4 consistent data points, this might be something

                # Find the rate that has the most matches
                frequent = (motion_rates[:,6] == ratecount_max)

                frequent_rates = motion_rates[:, 4:6][frequent]
                numpy.savetxt(sys.stdout, frequent_rates)

                # compute average rate as the average of all frequent rates
                avg_rate = numpy.median(frequent_rates, axis=0)
                print ratecount_max, avg_rate

                # Now identify all points with rates close to the average rate
                counterparts = src_tree.query_ball_point(avg_rate, r=1, p=2)


                position_in_sequence = 1
                moving_object_id += 1
                print >>results, "\n\n\n\n"
                print >>results, "%.7f %.7f %3d %3d %3d %4d " % (
                    catalog[ref_cat][obj1, SXcolumn['ra']],
                    catalog[ref_cat][obj1, SXcolumn['dec']],
                    moving_object_id,
                    position_in_sequence,
                    ref_cat,
                    obj1,
                )

                for p in range(len(counterparts)):
                    p_idx = counterparts[p]
                        
                    c2 = int(motion_rates[p_idx, 1])
                    o2 = int(motion_rates[p_idx, 3])
                    position_in_sequence += 1
                    print >>results, "%.7f %.7f %3d %3d %9.4f %9.4f %8.2f %8.2f %9.4f %9.4f" % (
                        catalog[c2][o2, SXcolumn['ra']],
                        catalog[c2][o2, SXcolumn['dec']],
                        moving_object_id,
                        position_in_sequence,
                        catalog[c2][o2, SXcolumn['ra']]-catalog[ref_cat][obj1, SXcolumn['ra']],
                        catalog[c2][o2, SXcolumn['dec']]-catalog[ref_cat][obj1, SXcolumn['dec']],
                        (catalog[c2][o2, SXcolumn['ra']]-catalog[ref_cat][obj1, SXcolumn['ra']])*3600./(mjd_hours[c2]-mjd_hours[ref_cat]),
                        (catalog[c2][o2, SXcolumn['dec']]-catalog[ref_cat][obj1, SXcolumn['dec']])*3600./(mjd_hours[c2]-mjd_hours[ref_cat]),
                        motion_rates[p_idx, 4],
                        motion_rates[p_idx, 5]
                    )

        #numpy.savetxt("motion_rates_%d" % (ref_cat), motion_rates)

    numpy.savetxt("motion_rates", motion_rates)

    logger.info("all work done!")

    podi_logging.shutdown_logging(options)


    sys.exit(0)
