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
import astropy.time
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

    


#
#
# To-do list
#
# - make sure every star is only part of a single tracklet
#
#


def format_source(src_entry, mjd, magname, filter_name='r', discovery=False):

    mag_name = "mag_aper_%s" % magname
    err_name = "mag_err_%s" % magname
    
    merr = src_entry[SXcolumn[err_name]]
    mag_digits = 2 if (merr < 0.02) else 1

    format_string = "%(id)5s%(designation)7s%(discovery)1s%(note1)1s%(note2)1s%(date)17s%(ra)12s%(dec)12s%(empty)9s%(mag).1f %(filter)1s%(empty)5s%(obscode)03d"

    j2k = ephem.Equatorial(math.radians(src_entry[SXcolumn['ra']]), 
                           math.radians(src_entry[SXcolumn['dec']]), 
                           epoch=ephem.J2000)

    t = astropy.time.Time([mjd], format='mjd', scale='utc')
    dt = t.datetime[0]
    formatted_date = "%04d %02d %09.6f" % (
        dt.year, dt.month, 
        dt.day + dt.hour/24. + dt.minute/1440. + (dt.second + dt.microsecond*1e-6)/86400.
    )
    
    str_ra = str(j2k.ra).replace(":", " ")[:12]
    str_dec = str(j2k.dec).replace(":", " ")[:12]
    ret = format_string % {
        'empty': "",
        'id': "",
        'designation': '',
        'discovery': '*' if discovery else '',
        'note1': '?',
        'note2': "C",
        'date': formatted_date,
        'ra': str_ra,
        'dec': str_dec,
        'mag': src_entry[SXcolumn[mag_name]],
        'filter': filter_name,
        'obscode': 695,
    }
    return ret

if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)

    logger = logging.getLogger("Astro(id)metry")

    sidereal_reference = get_clean_cmdline()[1]
    inputlist = get_clean_cmdline()[2:]

    make_catalogs(inputlist)
    
    catalog_filelist = []
    catalog = []
    mjd_hours = []
    source_id = []

    phot_ZP = 25.

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
        cat_data = cat_data[good_photometry]

        # Apply photometric zeropoint to all magnitudes
        for key in SXcolumn_names:
            if (key.startswith("mag_aper")):
                cat_data[:, SXcolumn[key]] += phot_ZP

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

    candidates = []

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
                
                # Only use sources with positive flags - we'll mask stars that 
                # are already matched to another tracklet by setting flags to -1 
                possible_match = (d_total < math.fabs(max_rate/diff_to_rate)) & \
                                 (d_total > min_distance) & \
                                 (catalog[ref_cat][obj1,SXcolumn['flags']] >= 0)
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

            source_data = []
            mjd_data = []
            formatted_data = []

            if (ratecount_max > 4):
                # We have 4 consistent data points, this might be something

                # Find the rate that has the most matches
                frequent = (motion_rates[:,6] == ratecount_max)

                frequent_rates = motion_rates[:, 4:6][frequent]
                # numpy.savetxt(sys.stdout, frequent_rates)

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
                formatted_data.append(
                    format_source(catalog[ref_cat][obj1], mjd_hours[ref_cat]/24., '3.0')
                )
                source_data.append(catalog[ref_cat][obj1])
                mjd_data.append(mjd_hours[ref_cat])

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
                    formatted_data.append(
                        format_source(catalog[c2][o2], mjd_hours[c2]/24., '3.0')
                    )
                    # Mark this star as used so we don;t try to re-use it for another tracklet
                    # This hopefully gets rid of duplicates.
                    catalog[c2][o2,SXcolumn['flags']] = -1.
                    source_data.append(catalog[c2][o2])
                    mjd_data.append(mjd_hours[c2])

                source_data = numpy.array(source_data)
                mjd_data = numpy.array(mjd_data) / 24.

                #numpy.savetxt(sys.stdout, source_data)
                #numpy.savetxt(sys.stdout, mjd_data)
                #print "\n\n\n\n"

                source_entry = {
                    "sex": source_data,
                    "mjd": mjd_data,
                    'rate': avg_rate,
                    'positions': source_data.shape[0],
                    'formatted': formatted_data,
                }
                candidates.append(source_entry)

        #numpy.savetxt("motion_rates_%d" % (ref_cat), motion_rates)
    logger.info("all search-related work done, associating with known objects !")

    # print candidates

    # Query MPC for all objects nearby
    logger.info("Getting list of known objects from MPC")
    import podi_mpchecker
    mpc_data = podi_mpchecker.get_mpc_catalog(sidereal_reference)
    print mpc_data['Name']

    reference_hdu = astropy.io.fits.open(sidereal_reference)
    mjd_reference = reference_hdu[0].header['MJD-OBS']
    print mjd_reference

    #print mpc_data['RA'].shape
    #print mpc_data['RA']

    mpc_ra = [None]*len(mpc_data['Name'])
    mpc_dec = [None]*len(mpc_data['Name'])
    # Convert the sexagesimal data from MPC into degrees
    ds9_known = open("asteroidmetry_known.reg", "w")
    print >>ds9_known, """\
# Region file format: DS9 version 4.1
global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5\
"""
    for i in range(len(mpc_data['Name'])):
        j2000 = ephem.Equatorial(mpc_data['RA'][i], mpc_data['DEC'][i], epoch=ephem.J2000)
        mpc_ra[i] = math.degrees(j2000.ra)
        mpc_dec[i] = math.degrees(j2000.dec)
        print >>ds9_known, 'circle(%f,%f,%f")' % (mpc_ra[i], mpc_dec[i], 3.) 
    ds9_known.close()

    #print mpc_ra
    #print mpc_dec

    #
    # Now go though each found asteroid, and see if it's known already
    #
    ds9_reg = open("asteroidmetry.reg", "w")
    print >>ds9_reg, """\
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5\
"""
    for cand in candidates:
        # Correct all coordinates to match the MJD of the reference frame

        # Get time difference in hours
        delta_mjd_hours = (cand['mjd'] - mjd_reference) * 24.
        # total motion
        delta_radec = numpy.array(cand['rate']).reshape((1,2)).repeat(delta_mjd_hours.shape[0], axis=0) * delta_mjd_hours.reshape((delta_mjd_hours.shape[0],1)).repeat(2, axis=1)
        # print delta_radec

        # Account for the cos(dec) factor
        cos_dec = numpy.cos(numpy.radians(cand['sex'][:,SXcolumn['dec']]))
        delta_radec[:,0] /= cos_dec

        ra_dec_corrected = cand['sex'][:,0:2] - delta_radec/3600.

        radec = numpy.median(ra_dec_corrected, axis=0)
        radec_std = numpy.std(ra_dec_corrected, axis=0)*3600.

        # print ra_dec_corrected
        print >>ds9_reg, 'circle(%f,%f,%f")' % (radec[0], radec[1], 10.) #13:22:17.482,-15:33:16.91,11.9449")
        # With this information, search for asteroids that are both 
        # nearby in Ra/Dec and in proper motion
        #
        # search radii: ra/dec: 5 arcsec, +/- 1''/hr

        d_ra_cat = (mpc_ra - radec[0]) * math.cos(math.radians(radec[1]))
        d_dec_cat = mpc_dec - radec[1]
        d_total = numpy.hypot(d_ra_cat, d_dec_cat) * 3600. # convert to arcsec
        # print d_total

        d_dra = mpc_data['dracosdec'] - cand['rate'][0]
        d_ddec = mpc_data['ddec'] - cand['rate'][1]
        d_dtotal = numpy.hypot(d_dra, d_ddec) # this already is in arcsec/hr
        # print d_dtotal

        likely_id = (d_total < 5.) & (numpy.fabs(d_dra) < 1.5) & (numpy.fabs(d_ddec) < 1.5)

        # numpy.savetxt(sys.stdout, likely_id)

        n_likelies = numpy.sum(likely_id)
        if (n_likelies == 0):
            logger.info("This is a new object")
            print >>ds9_reg, '# text(%f,%f) text={NEW: %.1f, %.1f}' % (radec[0], radec[1], cand['rate'][0], cand['rate'][1])

            # If this object is not known, draw markers at all detected positions
            for i in range(cand['sex'].shape[0]):
                print >>ds9_reg, 'point(%f,%f) # point=circle' % (cand['sex'][i,0], cand['sex'][i,1])

            print d_total
            print d_dtotal

        elif (n_likelies == 1):
            counterpart_name = (numpy.array(mpc_data['Name'])[likely_id])[0]
            logger.info("Found unique counterpart: %s" % (counterpart_name))
            print >>ds9_reg, '# text(%f,%f) text={%s (%.1f %.1f)}' % (radec[0], radec[1], counterpart_name, cand['rate'][0], cand['rate'][1])
        else:
            logger.info("Found more than one likely counterpart")
            print >>ds9_reg, '# text(%f,%f) text={%s}' % (radec[0], radec[1], "-?-?-")

        time.sleep(0.01)
        print radec
        print radec_std
        print cand['rate']
        print "\n".join(cand['formatted'])
        print "\n\n"


    ds9_reg.close()

    logger.info("Done, shutting down")
    podi_logging.shutdown_logging(options)


    sys.exit(0)
