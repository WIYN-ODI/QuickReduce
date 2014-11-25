#!/usr/bin/env python

#
# Copyright (C) 2014, Ralf Kotulla
#                     kotulla@uwm.edu
#
# All rights reserved
#

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
import multiprocessing
import copy

import podi_sitesetup as sitesetup

output_debugfiles = False

def parallel_sourceextractor(queue, dummy):

    logger = logging.getLogger("ParSEx")

    while (True):

        cmd = queue.get()
        if (cmd == None):
            queue.task_done()
            break

        sexcmd, fitsfile = cmd
        logger.info("Creating source catalog for %s" % (fitsfile))
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

        queue.task_done()

def make_catalogs(inputlist):

    jobqueue = multiprocessing.JoinableQueue()

    number_sex_runs = 0
    for fitsfile in inputlist:
        catfile = "%s.cat" % (fitsfile[:-5])

        if (os.path.isfile(catfile)):
            # Don't do anything if the catalog already exists
            continue

        logger.debug("Ordering source catalog for %s" % (fitsfile))
        sex_config_file = "%s/.config/diffphot.sex" % (sitesetup.exec_dir)
        parameters_file = "%s/.config/diffphot.sexparam" % (sitesetup.exec_dir)
        sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s" % (
            sitesetup.sextractor, sex_config_file, parameters_file, catfile, 
            fitsfile)
        
        jobqueue.put((sexcmd, fitsfile))
        number_sex_runs += 1

    logger.info("Ordered %d new source catalogs, please wait" % (number_sex_runs))

    #
    # Start worker processes
    #
    worker_args = (jobqueue, "")
    processes = []
    for i in range(sitesetup.number_cpus):
        p = multiprocessing.Process(target=parallel_sourceextractor, args=worker_args)
        p.start()
        processes.append(p)

        # also add a quit-command for each process
        jobqueue.put(None)
        
    #
    # wait until all work is done
    #
    jobqueue.join()
    logger.info("All source catalogs have been created!")

    return



def my_replace(inp, key, insert):
    out = ""+inp
    while (out.find(key) >= 0):
        begin_pos = out.find(key)
        end_pos = begin_pos + len(key)
        out = "%s%s%s" % (out[:begin_pos], insert, out[end_pos:])
    return out


RA = 0
DEC = 1
MJD = 2
DRA = 3
DDEC = 4

name_tag = "%name"

def differential_photometry(inputlist, source_coords, 
                            plot_title=None, plot_filename=None,
                            reference_star_frame=None,
                            output_catalog=None,
                            skipotas=[]):

    logger = logging.getLogger("DiffPhot")

    src_params = []
    src_names = []
    # print source_coords
    for src_coord in source_coords:
        # interpret source_coords
        try:
            sc_items = src_coord.split(",")
            src_ra = float(sc_items[0])
            src_dec = float(sc_items[1])

            src_dra, src_ddec, src_mjd = 0., 0., 0.

            src_name = None
            if (len(sc_items) >= 5):
                src_dra = float(sc_items[2])
                src_ddec = float(sc_items[3])
                if (os.path.isfile(sc_items[4])):
                    hdu = astropy.io.fits.open(sc_items[4])
                    src_mjd = hdu[0].header['MJD-OBS']
                else:
                    src_mjd = float(sc_items[4])
            if (len(sc_items) >= 6):
                src_name = sc_items[5].strip()
            elif (len(sc_items) == 3):
                src_name = sc_items[2].strip()

            if (not src_name == None):
                logger.info("Found source identifier: %s" % (src_name))

        except:
            logger.error("ERROR: Problem with interpreting coordinate string: %s" % (src_coord))
            logger.warning("Ignoring this source")
            continue
            pass

        src_params.append( [src_ra, src_dec, src_mjd, src_dra, src_ddec] )
        src_names.append(src_name)

    src_params = numpy.array(src_params)
    # print "\n"*10,src_names,"\n"*10

    # print "source coordinate data:"
    # numpy.savetxt(sys.stdout, src_params, "%14.6f")

    # print src_ra, src_dec, src_dra, src_ddec, src_mjd

    # Internally, all coordinates are corrected for cos(dec)
    cos_declination = math.cos(math.radians(src_dec))
    src_params[:,RA] *= cos_declination
    # src_ra *= cos_declination 

    # Make sure to add the reference frame to the list of frames for 
    # which we need source catalogs
    logger.info("Using reference frame: %s" % (reference_star_frame))
    if (not reference_star_frame == None and
        not reference_star_frame in inputlist):
        full_filelist = copy.copy(inputlist)
        full_filelist.append(reference_star_frame)
        make_catalogs(full_filelist)
    else:
        make_catalogs(inputlist)
    
#    return

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
            logger.debug("Read %d sources from %s" % (cat_data.shape[0], dumpfile))
        else:
            try:
                cat_data = numpy.loadtxt(catalog_filename)
                logger.debug("Read %d sources from %s" % (cat_data.shape[0], catalog_filename))
            except:
                pass
                continue
            cat_data.dump(dumpfile)

        
        # Correct for the absolute zeropoint
        magzero = 25.
        hdu = astropy.io.fits.open(fitsfile)
        magzero = hdu[0].header['MAGZERO']
        magzero = hdu[0].header['PHOTZP_X']

        # # Translate all extension numbers to proper OTAs
        # extensions = set(cat_data[:, SXcolumn['ota']])
        # for i in extensions:
        #     ota = hdu[i].header['OTA']
        #     cat_data[:, SXcolumn['ota']][cat_data[:, SXcolumn['ota']] == i] = ota

        # # Now eliminate all sources in OTAs marked for omission
        # if (not skipotas == []):
        #     for ota in skipotas:
        #         to_skip = cat_data[:, SXcolumn['ota']] == ota
        #         cat_data = cat_data[~to_skip]

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

    logger.info("All catalog files created and read...")
    logger.debug("Cos(dec) correction = %.10f" % (cos_declination))

    #############################################################################
    #
    # Open the reference star catalog and select a number of stars for 
    # the differential photometry reference
    #
    #############################################################################
    if (not reference_star_frame == None):
        ref_star_cat_file = reference_star_frame[:-5]+".cat"
        ref_stars = numpy.loadtxt(ref_star_cat_file)
        ref_stars[:, SXcolumn['ra']] *= cos_declination
        ref_stars[:, SXcolumn['mag_aper_2.0']:SXcolumn['mag_aper_12.0']+1] += 25. #(magzero - 25.0)
    else:
        ref_stars = catalog[0]

    print "top 25 reference star Ra/Dec:\n",ref_stars[:25,:2]

    # Now find sources that have good photometry in all frames that 
    # we can use as reference stars
    small_errors = (ref_stars[:, aperture_error_column] < 0.05) & \
                   (ref_stars[:, aperture_column] > 14)
    no_flags = ref_stars[:, SXcolumn['flags']] == 0
    ref_stars = ref_stars[(small_errors & no_flags)]

    # Sort by magnitude and pick the 10 brightest stars
    si = numpy.argsort(ref_stars[:, aperture_column])
    ref_sorted = ref_stars[si]
    print "magnitudes of reference stars:\n",ref_sorted[:, aperture_column]

    n_ref = 250
    phot_refs = ref_sorted[0:n_ref,0:2]
    print "Ra/Dec of 50 bright differential standards:\n", phot_refs[:,0:2]
    if (output_debugfiles): numpy.savetxt("diff_stds.cat", phot_refs)

    # Now go through all catalogs and find the source and the reference 
    # magnitudes

    # return

    output_file = open('output.test', 'w')
    match_radius = 2./3600. # 2 arcsec

    all_cats = []
    target_source = []
    target_mjd = []

    #
    # Create source catalogs for the source and n reference sources from each catalog
    #
    if (output_debugfiles): numpy.savetxt("src_params", src_params)
    for i_cat in range(len(catalog)):
        # print "\n"*5

        src_cat = catalog[i_cat]

        # turn the catalog coordinates into a KDTree
        cat_tree = scipy.spatial.cKDTree(catalog[i_cat][:,0:2])
        # print catalog[i_cat][:,0:2]

        # Search for the actual source
        # add here: fudge with coordinates to account for motion
        src_radec = src_params[:, RA:DEC+1] + (mjd_hours[i_cat]-src_params[:,MJD:MJD+1]*24.)*src_params[:, DRA:DDEC+1]/3600.
        # print "source coordinates in frame %d" % (i_cat+1)
        # numpy.savetxt(sys.stdout, src_radec)

        # src_radec = numpy.array([[src_ra, src_dec]]) + (mjd_hours[i_cat] - src_mjd*24.) * numpy.array([src_dra, src_ddec]) / 3600.
        # print src_radec * numpy.array([1./cos_declination, 1.])
        
        # print src_radec.shape, phot_refs.shape

        #########################################################################
        #
        # Determine the magnitude of the source
        #
        #########################################################################

        src_tree = scipy.spatial.cKDTree(src_radec)
        # Search for the source we are interested in
        # src_counterpart = cat_tree.query_ball_point(src_radec, r=match_radius, p=2)
#        src_counterpart = cat_tree.query_ball_tree(src_tree, r=match_radius, p=2)
        src_counterpart = src_tree.query_ball_tree(cat_tree, r=match_radius, p=2)

        src_data = numpy.empty((src_radec.shape[0], catalog[i_cat].shape[1]))
        src_data.fill(numpy.NaN)
        # print src_data

        # print src_counterpart

        for obj in range(src_radec.shape[0]):
            if (len(src_counterpart[obj]) <= 0):
                # We couldn't find any source counterpart
                continue
            elif (len(src_counterpart[obj]) > 1):
                # Found multiple counterparts
                # add some brain to pick the most likely one
                continue
            else:
                # 
                try:
                    src_data[obj,:] = (src_cat[src_counterpart[obj][0]])
                except:
                    podi_logging.log_exception()
                    continue
                #print 'src-data=',src_data.shape

        # print src_data

        # #########################################################################
        # #
        # # Determine the aperture correction from the source magnitude 
        # # to the 12'' aperture we consider for "total intensity"
        # #
        # #########################################################################

        # # Find the aperture magnitude that is valid
        # # try:
        # #     numpy.savetxt(sys.stdout, [src_data[SXcolumn['mag_aper_2.0']:SXcolumn['mag_err_12.0']+1]], " %8.4f")
        # # except:
        # #     pass
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
        if (output_debugfiles): numpy.savetxt("apercorr_%d" % (i_cat), aperture_diffs)

        aperture_diff_filtered = three_sigma_clip(aperture_diffs)
        aperture_correction = numpy.median(aperture_diff_filtered)
        aperture_correction_std = numpy.std(aperture_diff_filtered)

        # numpy.savetxt("apcorr%d.raw" % (i_cat), aperture_diffs)
        # numpy.savetxt("apcorr%d.filter" % (i_cat), aperture_diff_filtered)
        print "ap-corr/std:", inputlist[i_cat], aperture_correction, aperture_correction_std
        

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
        # print phot_data.shape
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
        # numpy.savetxt("photcat%d" % (i_cat), phot_data)
        #
        # also apply aperture correction to target source
        #
        src_data[:,SXcolumn['mag_auto']] = src_data[:,SXcolumn[src_aper_mag]] + aperture_correction
        src_data[:,SXcolumn['mag_err_auto']] = src_data[:,SXcolumn[src_aper_err]]


        #
        # Combine all catalogs for later
        #
        all_cats.append(phot_data)
        target_source.append(src_data)
        logger.debug("Added target sources from frame #%d" % (len(target_source)))
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

    print "target_source.shape", target_source.shape
    print "all_cats.shape", all_cats.shape

    print "mjds:", target_mjd
    first_mjd = numpy.min(target_mjd)

    # Dump all uncorrected source catalogs
    for src in range(target_source.shape[1]):
        if (output_debugfiles): numpy.savetxt("src.cat.%d" % src, target_source[:,src,:])


    all_cats_orig = all_cats.copy()
    #############################################################################
    #############################################################################
    ##
    ## Iteratively compute the correction
    ##
    #############################################################################
    #############################################################################
    col_name = 'mag_aper_4.0'
    mc = [SXcolumn['mag_aper_2.0'], SXcolumn['mag_aper_12.0']+1]
    mec = [SXcolumn['mag_err_2.0'], SXcolumn['mag_err_12.0']+1]
    all_mags = all_cats[:,:, mc[0]:mc[1]]

    for iteration in range(15):

        #############################################################################
        #
        # Now compute the median and filtered magnitude of each reference star
        #
        #############################################################################
        print all_cats[:,:, mc[0]:mc[1]].shape

        median_mags = bottleneck.nanmedian(all_cats[:,:, mc[0]:mc[1]], axis=0)
        std_mags = bottleneck.nanstd(all_cats[:,:, mc[0]:mc[1]], axis=0)
        print "MED/STD:", median_mags.shape, std_mags.shape
        if (output_debugfiles): numpy.savetxt("photcat_median_%d" % iteration, median_mags)
        if (output_debugfiles): numpy.savetxt("photcat_std_%d" % iteration, std_mags)

        # numpy.savetxt("median_mags.input.%d" % (iteration), all_cats[:,:, SXcolumn[col_name]])
        # numpy.savetxt("median_mags.corr.%d" % (iteration), all_cats[:,:, SXcolumn[col_name]]-median_mags[:,2])

        # all_cats[outlier, :] = numpy.NaN

        # print all_cats[:,0, mc[0]:mc[1]]
        # print median_mags[0], std_mags[0]

        #
        # Cut through the catalog and write out lightcurves 
        # for each of the reference stars
        #
        if (iteration == 14):
            for s in range(all_cats.shape[1]):
                logger.debug("Dumping reference catalog % 3d" % (s+1))

                dummy = all_cats[:,s,:].copy()
                # print dummy.shape
                dummy = numpy.append(dummy, 
                                     all_cats[:,s,0:2] - numpy.nanmean(all_cats[:,s,0:2], axis=0),
                                     axis=1)
                # print target_mjd.shape
                dummy = numpy.append(dummy, 
                                     target_mjd.reshape(target_mjd.shape[0],1)-first_mjd,
                                     axis=1)

                # compute median magnitudes in each aperture
                _mags = all_cats[:,s,SXcolumn['mag_aper_2.0']:SXcolumn['mag_err_2.0']]
                _median_mags = 0# bottleneck.nanmedian(_mags, axis=0)
                # print _mags.shape, _median_mags.shape
                _mags -= _median_mags
                dummy = numpy.append(dummy, _mags, axis=1)

                # dummy[:,:-5] = all_cats[:,s,:]
                # dummy[:,-3:-1] = all_cats[:,s,0:2] - numpy.nanmean(all_cats[:,s,0:2], axis=0)
                # dummy[:,-1] = target_mjd

                # dummy = numpy.append(dummy, diffphot_correction.reshape(target_mjd.shape[0],1), axis=1)

                dummy = numpy.append(dummy, numpy.ones((target_mjd.shape[0],1))*median_mags[s], axis=1)
                dummy = numpy.append(dummy, numpy.ones((target_mjd.shape[0],1))*std_mags[s], axis=1)

                # numpy.savetxt("lightcurve_ref%03d.cat" % (s+1), dummy)

        median_corrected_mags = all_cats[:,:, mc[0]:mc[1]] - median_mags
        if (output_debugfiles): numpy.savetxt("median_corrected_%d.ap2" % (iteration), median_corrected_mags[:,:,2])

        #
        # Now we have the differential magnitude correction for each individual star
        # Again, filter out only the corrections that seem to agree with the 
        # bulk of the other stars
        #

        diffphot_correction = bottleneck.nanmedian(median_corrected_mags, axis=1)\
                                        .reshape((median_corrected_mags.shape[0], 1, median_corrected_mags.shape[2]))
        if (output_debugfiles): numpy.savetxt("diffcorrection.%d" % (iteration), diffphot_correction[:,0,:])
        print "diffphot_correction.shape = ",diffphot_correction.shape
        print "all_mags.shape = ",all_mags.shape

        #diffphot_correction_3d = diffphot_correction
        #    diffphot_correction.reshape((all_mags.shape[0], 1, all_mags.shape[2]))
        #print "diffphot_correction_3d (pre-repeat)=",diffphot_correction_3d.shape
        #diffphot_correction_3d = diffphot_correction_3d.repeat(all_mags.shape[1], axis=1)
        #print "diffphot_correction_3d (post-repeat)=",diffphot_correction_3d.shape
        #print "median_corrected_mags.shape=", median_corrected_mags.shape

        dp_corrected = median_corrected_mags - diffphot_correction
        # dp_corrected = (median_corrected_mags.T-diffphot_correction).T

        for ap in range(dp_corrected.shape[2]):
            if (output_debugfiles): numpy.savetxt("diffcorrected.%d.ap%d" % (iteration, ap), dp_corrected[:,:,ap])
        
        #
        # Now involve some magic to mask out all outliers
        #
        dp_corr_std = bottleneck.nanstd(dp_corrected, axis=1)\
                                .reshape((dp_corrected.shape[0],1,dp_corrected.shape[2]))
        print "dp_corr_std.shape=",dp_corr_std.shape
        if (output_debugfiles): numpy.savetxt("diffcorrection-std.%d" % (iteration), dp_corr_std[:,0,:])
        
        # Mask out all datapoints that exceed the standard deviation
        # dp_corr_std_3d = dp_corr_std.reshape((all_mags.shape[0], 1, all_mags.shape[2])).repeat(all_mags.shape[1], axis=1)
        # print dp_corr_std_3d.shape
        outliers = numpy.fabs(dp_corrected) > 3*dp_corr_std
        
        all_cats[:,:, mc[0]:mc[1]][outliers] = numpy.NaN


    print diffphot_correction

    print  "\n"*5,"done with computing corrections","\n"*5

#    return

    # #############################################################################
    # #
    # # Compute the differential photometry correction
    # # 
    # #############################################################################

    # # Now we have the median magnitudes of each reference star.
    # # Additionally, all outliers have already been flagged as NaNs
    # print "median mags of all reference stars:\n", median_mags.flatten()

    # # Compute the deviation of each source from its median magnitude
    # mag_deviation = all_cats[:,:, SXcolumn[col_name]] - median_mags
    # print "mag-deviation:\n", mag_deviation.flatten()

    # # Compute the median deviation across all stars for each of the source catalogs
    # diffphot_correction = numpy.median(mag_deviation, axis=1).reshape((mag_deviation.shape[0],1))
    # print "diffphot_correction:\n", diffphot_correction.flatten()

    # print diffphot_correction.shape

    # ####
    # # write the diffphot correction for debugging
    # _x = numpy.append(target_mjd.reshape(target_mjd.shape[0],1)-first_mjd, 
    #                   diffphot_correction.reshape(target_mjd.shape[0],1),
    #                   axis=1)
    # numpy.savetxt("diffphot_corr.debug", _x)


    # #numpy.save("correction", diffphot_correction)
    # #numpy.savetxt("correction.dat", diffphot_correction)
    # #numpy.save("all_cats", all_cats)

    # # As verification, apply the correction to each of the reference star catalogs
    # print all_cats.shape
    # all_cats[:,:, SXcolumn['mag_auto']] -= diffphot_correction
    # #for i in range(all_cats.shape[0]):
    # #    numpy.savetxt("photcat%d.corr2" % (i), all_cats[i,:,SXcolumn['mag_auto']])



    #############################################################################
    #
    # Now apply the differential photometry correction to the actual source photometry
    #
    #############################################################################

    #numpy.savetxt("mjd", target_mjd)
    #numpy.savetxt("target_source.raw", target_source)
    print target_source.shape
    print target_source

    # Come up with a better way of storing the data in ascii files for 
    # post-processing with external files

    print target_source.shape
    print diffphot_correction.shape

    target_corrected = target_source.copy()

    for i in range(target_corrected.shape[1]):
        if (output_debugfiles): numpy.savetxt("target_source.prefixed.%d" % (i), target_corrected[:,i,:])

    # target_corrected[:, :, SXcolumn['mag_auto']] -= diffphot_correction
    target_corrected[:, :, mc[0]:mc[1]] -= diffphot_correction

    # Also add the scatter of the diff.phot. correction to the 
    # photometric uncertainty of the source

    print "\n"*5, "correcting uncertainties","\n"*2
    mag_errors = target_corrected[:, :, mec[0]:mec[1]]
    print mag_errors.shape, dp_corr_std.shape

    mag_errors_2 = numpy.hypot(mag_errors, dp_corr_std)
    target_corrected[:, :, mec[0]:mec[1]] = mag_errors_2

    for i in range(target_corrected.shape[1]):
        if (output_debugfiles): numpy.savetxt("target_source.postfixed.%d" % (i), target_corrected[:,i,:])


    print target_corrected.shape

    #
    # Write the outcome for each source to a text-file
    # add a header to the top of each file
    # also add MJD and some correction info to the output
    #
    output_file_line_labels = SXcolumn_descriptions
    output_file_line_labels.append("MJD")
    for i in range(len(SXapertures)):
        output_file_line_labels.append("differential photometry correction, aperture size %.1f [mag]" % (SXapertures[i]))
    for i in range(len(SXapertures)):
        output_file_line_labels.append("scatter of differential photometry correction, aperture size %.1f [mag]" % (SXapertures[i]))
    for i in range(target_corrected.shape[1]):

        if (output_catalog == None):
            filename = "target_source.postfixed.%d" % (i)
        else:
            filename = my_replace(output_catalog, name_tag, src_names[i])
        print "putput catalog:", output_catalog, filename
        filename = my_replace(filename, " ", "_")

        fn = open(filename, "w")
        for line in range(len(output_file_line_labels)):
            print >>fn, "Line %02d: %s" % (line+1, output_file_line_labels[line])
        _dum = target_corrected[:,i,:].copy()
        print _dum.shape, target_mjd.shape, target_mjd.reshape((-1,1)).shape
        _dum = numpy.append(_dum, target_mjd.reshape((-1,1))/24., axis=1)
        _dum = numpy.append(_dum, diffphot_correction[:,0,:], axis=1)
        _dum = numpy.append(_dum, dp_corr_std[:,0,:], axis=1)
        numpy.savetxt(fn, _dum)
        fn.close()

    #############################################################################
    #
    # Create some plots for the light curve before & after correction
    #
    #############################################################################

    print src_names
    for i_source in range(src_params.shape[0]):

        try:
            fig = matplotlib.pyplot.figure()
            #ax_final = fig.add_subplot(311)
            #ax_raw = fig.add_subplot(312)
            #ax_corr = fig.add_subplot(313)

            time = target_mjd  /24.
            #print time.shape
            #print target_source.shape
            #print target_source[0,:,0].shape
            # break

            width = 0.83
            ax_final = fig.add_axes([0.15,0.3,width,0.65])
            ax_corr = fig.add_axes([0.15,0.1,width,0.2])

            if (not plot_title == None):
                this_plot_title = plot_title

                if (not src_names[i_source] == None):
                    print src_names[i_source] 
                    this_plot_title = my_replace(plot_title, name_tag, src_names[i_source])

                    # begin_pos = plot_title.find(name_tag)
                    # end_pos = begin_pos + len(name_tag)
                    # this_plot_title = "%s%s%s" % (plot_title[:begin_pos], src_names[i_source], plot_title[end_pos+1:])

                ax_final.set_title(this_plot_title)

            from matplotlib.ticker import NullFormatter
            nullfmt   = NullFormatter()         # no labels
            # no labels
            ax_final.xaxis.set_major_formatter(nullfmt)
            ax_corr.set_xlabel("modified julian date")

            #
            # Also plot the uncorrected light curve as crossed connected with a thin line
            #
            c = '#A75858'
            ax_final.scatter(time, target_source[:, i_source, SXcolumn['mag_auto']],
                             color=c, marker='x', label='aperture-corr.')
            ax_final.plot(time, target_source[:, i_source, SXcolumn['mag_auto']],
                             color=c, linestyle='-')

            #
            # Also show the actual 4'' aperture data
            #
            c = '#4AA2A5'
            ax_final.scatter(time, target_source[:, i_source, SXcolumn['mag_aper_4.0']],
                             color=c, marker='x', label='raw 4"')
            ax_final.plot(time, target_source[:, i_source, SXcolumn['mag_aper_4.0']],
                             color=c, linestyle='-')

            #
            # Plot the light curve as data points with uncertainties
            # 
            ax_final.errorbar(time, target_corrected[:, i_source, SXcolumn[col_name]], 
                        yerr=target_corrected[:, i_source, SXcolumn['mag_err_4.0']], 
                        fmt='o', c='blue')
            ax_final.scatter(time, target_corrected[:, i_source, SXcolumn[col_name]],
                             c='blue', label='final')
            ax_final.set_ylabel("app. mag.\n[mag (12'')]")
            ax_final.legend(loc='best', borderaxespad=0.5, prop={'size':9})

            #
            # As bottom plot show correction amplitude
            #
            # Determine the plotting range
            max_corr = numpy.nanmax(numpy.fabs(diffphot_correction[:,0,:])) * 1.35
            ax_corr.plot(time, diffphot_correction[:, 0,:], marker='o')
            # ax_corr.set_ylim((-0.09,+0.09))
            ax_corr.set_ylim((-1*max_corr,+1*max_corr))
            ax_corr.axhline(y=0, linewidth=1, color='grey')
            ax_corr.set_ylabel("correction\n[mag]")

            if (plot_filename == None):
                matplotlib.pyplot.show()
                fig.show()
            else:
                if (not src_names[i_source] == None and plot_filename.find('%name') >= 0):
                    names_underscore = my_replace(src_names[i_source], " ", "_")
                    this_plot = my_replace(plot_filename, "%name", names_underscore)
                elif (plot_filename.find('%name') >= 0):
                    name_id = "src_%d" % (i_source+1)
                    this_plot = my_replace(plot_filename, "%name", name_id)
                else:
                    this_plot = "%s.%d.%s" % (plot_filename[:-4], i_source+1, plot_filename[-3:])

    #            fig.savefig(plot_filename)
                fig.savefig(this_plot)
                logger.info("Writing plot: %s"  % (this_plot))
        except:
            podi_logging.log_exception()
            pass

    return



if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)
    logger = logging.getLogger("DiffPhot: Main")

    print sys.argv
    print " ".join(sys.argv)

    clean_cmdline = get_clean_cmdline()[1:]
    coord_end = 0
    for i in range(len(clean_cmdline)):
        if (clean_cmdline[i] == ":"):
            coord_end = i
            break

    skipota = []
    _skipota = cmdline_arg_set_or_default("-skipota", None)
    if (not _skipota == None):
        items = _skipota.split(",")
        for i in items:
            skipota.append(int(i))
    logger.info("Excluding all sources from these OTAs: %s" % (
        ", ".join(skipota)))
    source_data = clean_cmdline[:coord_end]
    inputlist = clean_cmdline[coord_end+1:]

    # source_data = get_clean_cmdline()[1]
    # catalog_dir = None #get_clean_cmdline()[2]
    # inputlist = get_clean_cmdline()[2:]

    print "source_data:"
    print "\n".join(source_data)
    print "\n\ninput-list:"
    print "\n".join(inputlist)

    # print "source data", source_data
    # print "filelist",inputlist

    title = None
    if (cmdline_arg_isset('-title')):
        title=get_cmdline_arg('-title')
    plot_filename = cmdline_arg_set_or_default('-plotfile', None)

    reference_frame = cmdline_arg_set_or_default('-refframe', None)

    output_catalog = cmdline_arg_set_or_default('-catout', None)

    data = differential_photometry(inputlist, source_coords=source_data,
                                   plot_title=title, plot_filename=plot_filename,
                                   reference_star_frame=reference_frame,
                                   output_catalog=output_catalog,
                                   skipotas=skipota,
    )

    logger.info("Done, shutting down")
    podi_logging.shutdown_logging(options)


    sys.exit(0)
