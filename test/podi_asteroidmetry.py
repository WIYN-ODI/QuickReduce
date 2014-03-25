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


def format_source(src_entry, mjd, magname, filter_name='V', discovery=False, designation='ODI1'):

    mag_name = "mag_aper_%s" % magname
    err_name = "mag_err_%s" % magname
    
    merr = src_entry[SXcolumn[err_name]]
    mag_digits = 2 if (merr < 0.02) else 1

    format_string = "%(id)5s%(designation)7s%(discovery)1s%(note1)1s%(note2)1s%(date)16s%(ra)12s%(dec)12s%(empty)10s%(mag).1f %(filter)1s%(empty)6s%(obscode)03d"

    j2k = ephem.Equatorial(math.radians(src_entry[SXcolumn['ra']]), 
                           math.radians(src_entry[SXcolumn['dec']]), 
                           epoch=ephem.J2000)

    t = astropy.time.Time([mjd], format='mjd', scale='utc')
    dt = t.datetime[0]
    formatted_date = "%04d %02d %08.5f" % (
        dt.year, dt.month, 
        dt.day + dt.hour/24. + dt.minute/1440. + (dt.second + dt.microsecond*1e-6)/86400.
    )
    
    str_ra = str(j2k.ra).replace(":", " ")[:12]
    str_dec = str(j2k.dec).replace(":", " ")[:12]
    ret = format_string % {
        'empty': "",
        'id': "",
        'designation': designation,
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


def write_mpc_header(mpc):

    print >>mpc, "COD 695" # Kitt Peak
    print >>mpc, "CON R S. McMillan, Univ. of Arizona, 1629 E. Univ. Blvd, Tucson AZ 85721"
    print >>mpc, "CON [bob@lpl.arizona.edu]"
    print >>mpc, "OBS --- add info here ---"
    print >>mpc, "MEA --- add names here ---"
    print >>mpc, "TEL WIYN 3.5-m f/6.3 reflector + CCDs"
    print >>mpc, "NET 2MASS"
    print >>mpc, "COM Astrometric errors are generally 0.6-0.8 arcsec including WCS fitting"
    print >>mpc, "COM and centroiding."
    print >>mpc, "BND r"
    print >>mpc, "ACK WIYN 3.5-m followup"
    print >>mpc, "AC2 bob@lpl.arizona.edu"



def find_moving_objects(sidereal_reference, inputlist, min_count, min_rate, mpcfile=None):

    logger = logging.getLogger("Astro(id)metry")

    make_catalogs(inputlist)
    
    catalog_filelist = []
    catalog = []
    mjd_hours = []
    source_id = []
    filters = []

    phot_ZP = 25.

    # Now open all files, get time-stams etc
    for fitsfile in inputlist:

        catalog_filename = fitsfile[:-5]+".cat"
        cat_data = numpy.loadtxt(catalog_filename)
        
        hdulist = astropy.io.fits.open(fitsfile)
        obs_mjd = hdulist[0].header['MJD-OBS']
        exptime = hdulist[0].header['EXPTIME']
        mjd_middle = (obs_mjd * 24.) + (exptime / 3600.)
        filter_name = hdulist[0].header['FILTER']

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
        filters.append(filter_name)

        print catalog_filename, cat_data.shape


    max_rate = 100 # arcsec / hour
    min_distance = 2. / 3600. # arcsec
    motion_rates = numpy.zeros((0,6))

    # Now select all neighbors within the search radius
    radius = 1.

    results = open("results", "w")
    moving_object_id = 0

    candidates = []

    for ref_cat in range(len(catalog)-1-min_count):
        logger.debug("Checking out catalog %d" % (ref_cat))

        for obj1 in range(catalog[ref_cat].shape[0]):
            # Take one source in first catalog; 
            # compute difference and rate from each source in the second 
            # catalog to the source in the first catalog.

#            if (len(candidates) > 75):
#                break

            if (catalog[ref_cat][obj1,SXcolumn['flags']] < 0):
                continue

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
                                 (catalog[second_cat][:,SXcolumn['flags']] >= 0)
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

            # filename = "motion_rates_%02d_%04d.dump" % (ref_cat, obj1)
            # numpy.savetxt(filename, motion_rates)

            # numpy.savetxt(sys.stdout, motion_rates)

            #
            # Now we know how often each star appears
            #

            # check often the most frequent rate appears
            ratecount_max = numpy.max(motion_rates[:,6])

            source_data = []
            mjd_data = []
            formatted_data = []

            if (ratecount_max > min_count):
                # We have 4 consistent data points, this might be something

                # Save the motion rates and counts for later
                fn = "motion_rates_cand%d" % (len(candidates)+1)
                numpy.savetxt(fn, motion_rates)

                # Find the rate that has the most matches
                frequent = (motion_rates[:,6] == ratecount_max)

                frequent_rates = motion_rates[:, 4:6][frequent]
                # numpy.savetxt(sys.stdout, frequent_rates)

                # compute average rate as the average of all frequent rates
                avg_rate = numpy.median(frequent_rates, axis=0)
                logger.info("Found new candidate % 3d (#src=%d [% 2d/% 2d]): %.3f %.3f" % (
                    len(candidates)+1, ref_cat+1, len(catalog), ratecount_max+1, avg_rate[0], avg_rate[1])
                )

                # Now identify all points with rates close to the average rate
                counterparts = src_tree.query_ball_point(avg_rate, r=1, p=2)

                # Determine all source catalogs that contribute to this tracklets
                tracklet = motion_rates[counterparts]
                source_catid = tracklet[:, 1]
                source_cat_set = set(source_catid)
                tracklet_valid = source_catid > 0
                # Now make sure to exclude all sources in the tracklet that appear multiple times
                for srccat in source_cat_set:
                    from_this_cat = (source_catid == srccat)
                    if (numpy.sum(from_this_cat) > 1):
                        # multiple sources are a no-no
                        # logger.warning("Found multiple sources (%d) from the same catalog" % (numpy.sum(from_this_cat)))
                        tracklet_valid = tracklet_valid & (from_this_cat == False)

                # Now check again if we still have more sources than we require
                if (numpy.sum(tracklet_valid) < min_count):
                    logger.debug("Excluding tracklet with too many same-catalog sources (%d --> %d)" % (tracklet.shape[0], numpy.sum(tracklet_valid)))
                    continue

                tracklet = tracklet[tracklet_valid]

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

                # for p in range(len(counterparts)):
                #     p_idx = counterparts[p]
                #     c2 = int(motion_rates[p_idx, 1])
                #     o2 = int(motion_rates[p_idx, 3])

                for p in range(tracklet.shape[0]): 
                    c2 = int(tracklet[p, 1])
                    o2 = int(tracklet[p, 3])

                    position_in_sequence += 1
                    print >>results, "%.7f %.7f %3d %3d %9.4f %9.4f %8.2f %8.2f %9.4f %9.4f %.5f %.5f" % (
                        catalog[c2][o2, SXcolumn['ra']],
                        catalog[c2][o2, SXcolumn['dec']],
                        moving_object_id,
                        position_in_sequence,
                        catalog[c2][o2, SXcolumn['ra']]-catalog[ref_cat][obj1, SXcolumn['ra']],
                        catalog[c2][o2, SXcolumn['dec']]-catalog[ref_cat][obj1, SXcolumn['dec']],
                        (catalog[c2][o2, SXcolumn['ra']]-catalog[ref_cat][obj1, SXcolumn['ra']])*3600./(mjd_hours[c2]-mjd_hours[ref_cat]),
                        (catalog[c2][o2, SXcolumn['dec']]-catalog[ref_cat][obj1, SXcolumn['dec']])*3600./(mjd_hours[c2]-mjd_hours[ref_cat]),
                        tracklet[p, 4],
                        tracklet[p, 5],
                        mjd_hours[c2], mjd_hours[ref_cat]
                    )
                    formatted_data.append(
                        format_source(catalog[c2][o2], mjd_hours[c2]/24., '3.0')
                    )
                    # Mark this star as used so we don;t try to re-use it for another tracklet
                    # This hopefully gets rid of duplicates.
                    catalog[c2][o2,SXcolumn['flags']] = -1.
                    source_data.append(catalog[c2][o2])
                    mjd_data.append(mjd_hours[c2])

                # Collect additional information for the stars in the most frequent pattern
                dum = numpy.zeros(shape=(tracklet.shape[0], tracklet.shape[1]+6))
                dum[:,:tracklet.shape[1]] = tracklet
                for i in range(dum.shape[0]):
                    c1, o1 = int(dum[i, 0]), int(dum[i,2])
                    c2, o2 = int(dum[i, 1]), int(dum[i,3])
                    dum[i,7:9] = catalog[c1][o1,0:2]
                    dum[i,9:11] = catalog[c2][o2,0:2]
                    dum[i,11] = mjd_hours[c1]
                    dum[i,12] = mjd_hours[c2]
                fn += ".frequent"
                numpy.savetxt(fn, dum)

                # dump all source information for all sources in this tracklet
                fn = "motion_rates_cand%d.source" % (len(candidates)+1)
                dfn = open(fn, "w")
                print >>dfn, mjd_hours[ref_cat]," ",
                numpy.savetxt(dfn, [catalog[ref_cat][obj1]])
                cp = motion_rates[counterparts]
                for i in range(tracklet.shape[0]): #len(counterparts)):
                    c2, o2 = int(tracklet[i, 1]), int(tracklet[i,3])
                    print >>dfn, mjd_hours[c2]," ",
                    numpy.savetxt(dfn, [catalog[c2][o2]])
                dfn.close()

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
                    'tracklet': tracklet,
                }
                candidates.append(source_entry)

        #numpy.savetxt("motion_rates_%d" % (ref_cat), motion_rates)
    logger.info("all search-related work done, associating with known objects !")

    # print candidates

    # Query MPC for all objects nearby
    logger.info("Getting list of known objects from MPC")
    import podi_mpchecker
    mpc_data = podi_mpchecker.get_mpc_catalog(sidereal_reference)
    # print mpc_data['Name']

    reference_hdu = astropy.io.fits.open(sidereal_reference)
    mjd_reference = reference_hdu[0].header['MJD-OBS']
    logger.debug("MJD of reference frame: %f" % (mjd_reference))

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

    # If requested, also write the information in the MPC-style into a separate MPC file
    mpc = None
    if (not mpcfile == None):
        mpc = open(mpcfile, "w")
        
        # Write some basic info into MPC file
        write_mpc_header(mpc)

    new_discovery_number = 0
    for cand_id in range(len(candidates)):

        cand = candidates[cand_id]
        logger = logging.getLogger("Candidate(%d)" % (cand_id+1))

        # Check if the average rate is larger than the minimum rate
        combined_rate = numpy.hypot(cand['rate'][0], cand['rate'][1])
        if (combined_rate < min_rate):
            logger.info("Ignoring candidate that does not meet the minimum rate requirement")
            continue

        # Sort the source-data by time (MJD)
        si = numpy.argsort(cand['mjd'])
        mjd_sorted = cand['mjd'][si]
        sex_sorted = cand['sex'][si]
        #print cand['mjd']
        #print mjd_sorted
        #print sex_sorted

        radec_start = sex_sorted[0, 0:2]
        radec_end = sex_sorted[-1, 0:2]
        radec_total = (radec_end - radec_start) * 3600.
        radec_total[0] *= math.cos(math.radians(numpy.median(sex_sorted[:,1])))
        rate_total = radec_total / ((mjd_sorted[-1]-mjd_sorted[0])/24.0) 
        if (numpy.hypot(rate_total[0], rate_total[1]) < min_rate):
            logger.debug("Ignoring candidate that does not meet the global minimum rate requirement (%.3f %.3f)" % (rate_total[0], rate_total[1]))
            continue

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

        new_discovery = False
        n_likelies = numpy.sum(likely_id)
        logger.info("This is candidate # %d ..." % (cand_id+1))

        counterpart_name = 'new'
        comment = "none"
        designation = ""
        comment = ""
        if (n_likelies == 0):
            logger.info("This is a new object")
            print >>ds9_reg, '# text(%f,%f) text={NEW (%d): %.1f, %.1f}' % (radec[0], radec[1], cand_id+1, cand['rate'][0], cand['rate'][1])

            # If this object is not known, draw markers at all detected positions
            for i in range(cand['sex'].shape[0]):
                print >>ds9_reg, 'point(%f,%f) # point=circle' % (cand['sex'][i,0], cand['sex'][i,1])

            logger.info(d_total)
            logger.info(d_dtotal)
            new_discovery = True
            new_discovery_number += 1
            counterpart_name = "new_%d" % (cand_id+1)
            designation = "ODI%d" % (new_discovery_number)
            comment = "This is a new discovery"
        elif (n_likelies == 1):
            counterpart_name = (numpy.array(mpc_data['Name'])[likely_id])[0]
            comment = (numpy.array(mpc_data['comment'])[likely_id])[0] if 'comment' in mpc_data else "none"
            logger.info("Found unique counterpart: %s" % (counterpart_name))
            if (not comment == None):
                logger.info("MPC Comment: %s" % (comment))
            print >>ds9_reg, '# text(%f,%f) text={%s (%.1f %.1f)}' % (radec[0], radec[1], counterpart_name, cand['rate'][0], cand['rate'][1])
            designation = podi_mpchecker.mpc_name2id(counterpart_name)
            comment = "Incidental astrometry for %s" % (counterpart_name)
        else:
            counterpart_name = 'multiple'
            logger.info("Found more than one likely counterpart")
            print >>ds9_reg, '# text(%f,%f) text={%s}' % (radec[0], radec[1], "-?-?-")
            comment = "More than one plausible counterpart identified"

        candidates[cand_id]['new_discovery'] = new_discovery
        candidates[cand_id]['name'] = counterpart_name
        candidates[cand_id]['comment'] = comment if (not comment == "") else "none"

        time.sleep(0.01)
        logger.info("ra/dec: %s" % (str(radec)))
        logger.info("ra/dec corrected scatter: %s" % (str(radec_std)))
        logger.info("rate: %s" % (str(cand['rate'])))
        logger.info("global rate: %.3f %.3f" % (rate_total[0], rate_total[1]))
        
        # print "\n".join(cand['formatted'])
        # numpy.savetxt(sys.stdout, cand['sex'][:,0:4])
        
        # Now output the results in the MPC format
        all_fs = []
        all_fs.append("COM %s" % (comment))
        average_mag_error = numpy.median(sex_sorted[:, SXcolumn['mag_err_3.0']])
        all_fs.append("COM average photometric uncertainty: %.2f mag" % (average_mag_error))
        for i_src in range(mjd_sorted.shape[0]):
            
            fs = format_source(sex_sorted[i_src], mjd_sorted[i_src], '3.0', 
                               discovery=new_discovery, 
                               designation=designation, 
                               filter_name='r')
            new_discovery = False
            all_fs.append(fs)
        logger.info("for MPC:\n------\n%s\n------\n\n" % ("\n".join(all_fs)))
        if (not mpc == None):
            mpc.write(os.linesep.join(all_fs)+os.linesep)

    ds9_reg.close()
    if (not mpc == None):
        mpc.close()

    return candidates



if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)
    logger = logging.getLogger("Asteroidmetry: Main")

    min_count = int(cmdline_arg_set_or_default('-mincount', 5))
    min_rate = float(cmdline_arg_set_or_default('-minrate', 2))
    mpc_file = cmdline_arg_set_or_default('-mpc', None)

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    if (cmdline_arg_isset("-writeregion")):

        ref_file = get_clean_cmdline()[1]
        region_file = get_clean_cmdline()[2]

        with open("asteroidmetry.pickle", "rb") as pf:
            candidates = pickle.load(pf)

        for cand in candidates:
            print cand['rate']

        ds9_reg = open(region_file, "w")
        print >>ds9_reg, """\
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5\
        """
        
        hdulist = astropy.io.fits.open(ref_file)
        rel_ra = rel_dec = rel_mjd = 0
        rel_ra = hdulist[0].header['_NSID_RA'] if ('_NSID_RA' in hdulist[0].header) else 0.0
        rel_dec = hdulist[0].header['_NSID_DE'] if ('_NSID_DE' in hdulist[0].header) else 0.0
        rel_mjd = hdulist[0].header['_NSID_RT'] if ('_NSID_RT' in hdulist[0].header) else \
                  hdulist[0].header['MJD-OBS'] if ('MJD-OBS' in hdulist[0].header) else 0.0

        for cand in candidates:
            
            d_mjd_hours = (cand['mjd'] - rel_mjd) * 24.
            
            dra = (rel_ra * d_mjd_hours) 
            ddec = (rel_dec * d_mjd_hours)
            if (not 'name' in cand):
                continue

            if ('name' in cand): print cand['name']
            #print d_mjd_hours
            #print dra
            #print ddec

            dec_fixed = cand['sex'][:,1] - ddec / 3600.
            ra_fixed = cand['sex'][:,0] - (dra / numpy.cos(numpy.radians(dec_fixed))) / 3600.
            
            # Draw a circle for each detection
            for i in range(cand['sex'].shape[0]):
                print >>ds9_reg, 'point(%f,%f) # point=circle' % (ra_fixed[i], dec_fixed[i])
            
            # Add some label with the name
            # Use the coordinates of the first recorded point, and correct 
            # it to the MJD reference of the reference frame
            dmjd0 = (cand['mjd'][0] - rel_mjd) * 24.
            dra0 = (rel_ra - cand['rate'][0]) * dmjd0 / 3600.
            ddec0 = (rel_dec - cand['rate'][1]) * dmjd0 / 3600.
            dec0 = cand['sex'][0,1] - ddec0
            ra0 = cand['sex'][0,0] - dra0 / math.cos(math.radians(dec0))

            try:
                print >>ds9_reg, '# text(%f,%f) text={%s}' % (ra0, dec0, cand['name'])
            except:
                pass
        ds9_reg.close()

                
    else:
        sidereal_reference = get_clean_cmdline()[1]
        inputlist = get_clean_cmdline()[2:]

        candidates = find_moving_objects(sidereal_reference, inputlist, min_count, min_rate, mpc_file)

        with open("asteroidmetry.pickle", "wb") as pf:
            pickle.dump(candidates, pf)
            
        # print candidates[0]


    logger.info("Done, shutting down")
    podi_logging.shutdown_logging(options)


    sys.exit(0)
