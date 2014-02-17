#! /usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestions on how to improve the code or its 
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
This module handles all functionality related to saturation and persistency
effects. Most functions are called during reduction from within collectcells.


Standalone functions and command-line flags
-------------------------------------------
* **-makecat**

  ``podi_persistency -makecat (-persistency=dir/) file1.fits file2.fits``

  Create saturation catalogs for a number of files and write results to the
  directory given with the -persistency flag.


* **-masksattrails**

  ``podi_persistency -masksattrails input.fits catalog.fits output.fits``

  Apply the persistency masking to the specified input file, using the
  saturation table from file catalog.fits and write the resulting file into
  output.fits. This assumes that the input.fits file is a valid file created
  with collectcells.

* **-findclosemjds**

  ``podi_persistency -findclosemjds (-persistency=/dir) input.fits``

  Test-routine to find saturation catalogs within a fixed range of 
  [-1,600] seconds around the MJD of the specified input frame.


* **-fixpersistency**

  ``podi_persistency -fixpersistency (-persistency=/dir) input.fits output.fits``

  Similar to the -masksettrails functionality, but using all files within a
  fixed MJD range ([-1,1800] seconds) around the MJD of the input frame. Results
  are written to output.fits. As above it is assumed that input.fits is a valid
  frame created with collectcells.

Modules
-------

"""

import sys
import os
import pyfits
import numpy
import scipy
import pywcs
from astLib import astWCS
import jdcal

#import podi_plotting
import matplotlib
import matplotlib.pyplot

import time
import multiprocessing
import Queue
import itertools

from podi_definitions import *
import podi_sitesetup as sitesetup
import podi_collectcells
import podi_logging
import logging

try:
    import cPickle as pickle
except:
    import pickle


SXFocusColumn = {
    "ra": 0,
    "dec": 1,
    "x": 2, 
    "y": 3,
    "extension": 4,
    "mag_auto": 6,
    "fwhm_image": 7,
    "fwhm_world": 8
}

def mp_measure_focus(queue_in, queue_ret, verbose=False):
    """
    This is a small helper routine for the process of creating the saturation
    catalogs.  It reads filenames from job queue, creates the arrays of pixel
    coordinates, and posts the results to a return queue. Actually creating the
    fits tables is then handled by the main process.

    Parameters
    ----------
    queue_in : Queue

        Holds all the input files

    queue_ret : Queue

        Queue to report results back to main process

    """

    logger = logging.getLogger("MakeSetCat")
    logger.debug("Starting worker process")

    while (True):
        task = queue_in.get()

        if (task == None):
            queue_in.task_done()
            logger.debug("Received shutdown command, terminating")
            return

        filename, n_stars = task

        # print "\n"*10,filename,"\n"*10
        cat_name = measure_focus_ota(filename, n_stars)

        queue_ret.put( cat_name )

        queue_in.task_done()

    return





def measure_focus_ota(filename, n_stars=5):
    """
    Create a saturation table for a given OTA exposure.

    """


    # print"\n\n\nworking on file ",filename,"\n\n\n"

    try:
        hdulist = pyfits.open(filename)
    except IOError:
        logger.debug("Can't open file %s" % (filename))
        return None
    except:
        podi_logging.log_exception()
        return None

    obsid = hdulist[0].header['OBSID']
    ota = hdulist[1].header['WN_OTAX'] * 10 + hdulist[1].header['WN_OTAY']

    logger = logging.getLogger("MeasureFocusOTA: %s(%02d)" % (obsid, ota))
    logger.info("Starting work ...")

    obsid = hdulist[0].header['OBSID']
    ota = int(hdulist[0].header['FPPOS'][2:4])

    # Run SourceExtractor on the file
    basedir, _ = os.path.split(os.path.abspath(sys.argv[0]))
    sex_config = "%s/.config/focus.sexconf" % (basedir)
    sex_param = "%s/.config/focus.sexparam" % (basedir)
    catfile = "%s/tmp.%s_OTA%02d.cat" % (sitesetup.scratch_dir, obsid, ota)
    sex_cmd = "%(sexcmd)s -c %(sex_config)s -PARAMETERS_NAME %(sex_param)s -CATALOG_NAME %(catfile)s %(filename)s" % {
        "sexcmd": sitesetup.sextractor,
        "sex_config": sex_config,
        "sex_param": sex_param,
        "catfile": catfile,
        "filename": filename,
        "redirect": sitesetup.sex_redirect,
    }
    # print sex_cmd
    # Run source extractor
    # catfile = "/tmp//tmp.pid4383.20121008T221836.0_OTA33.cat"

    if (not os.path.isfile(catfile)):
        logger.debug("Running source extractor to search for stars")
        start_time = time.time()
        try:
            ret = subprocess.Popen(sex_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
            (sex_stdout, sex_stderr) = ret.communicate()
            if (ret.returncode != 0):
                logger.warning("Sextractor might have a problem, check the log")
                logger.info("Stdout=\n"+sex_stdout)
                logger.info("Stderr=\n"+sex_stderr)
        except OSError as e:
            podi_logging.log_exception()
            print >>sys.stderr, "Execution failed:", e
        end_time = time.time()
        logger.debug("SourceExtractor finished after %.2f seconds" % (end_time - start_time))
    else:
        logger.debug("Source catalog already exists, re-using the old file")

    #
    # delete the tmp catalog
    #

    #
    # Load the source catalog.
    # Handle cases of non-existing or empty catalogs
    # 
    logger.debug("loading the source catalog from %s" % (catfile))
    try:
        source_cat = numpy.loadtxt(catfile)
    except IOError:
        logger.warning("The Sextractor catalog is empty, ignoring this OTA")
        source_cat = None
        return None
    if (source_cat.shape[0] <= 0):
        # no sources found
        return None

    # print "\n\n total sources in raw file",source_cat.shape,"\n\n"
    logger.debug("Found %d sources in raw SourceExtractor catalog %s" % (source_cat.shape[0], catfile))

    #
    # Now convert all X/Y values to proper OTA X/Y coordinates based on their 
    # extension number
    #
    corr_cat = None
    # print "extensions:", source_cat[:,SXFocusColumn['extension']]

    for i in range(len(hdulist)):
        if (not is_image_extension(hdulist[i])):
            # print "skipping extension",i,", this is not an image extension"
            continue

        # get cell_x, cell_y from this extension
        cell_x = hdulist[i].header['WN_CELLX']
        cell_y = hdulist[i].header['WN_CELLY']
        x1, c2, y1, y2 = cell2ota__get_target_region(cell_x, cell_y, 1)
        
        # Create a mask for all sources in this cell
        in_this_cell = (source_cat[:,SXFocusColumn['extension']] == i)
        if (numpy.sum(in_this_cell) <= 0):
            logger.debug("Couldn't find any sources in cell %d,%d" % (cell_x, cell_y))
            continue

        # print "found",numpy.sum(in_this_cell),"sources for cell",cell_x, cell_y, "  adding", x1, y1
        
        cell_cat = source_cat[in_this_cell]
        cell_cat[:,SXFocusColumn['x']] += x1
        cell_cat[:,SXFocusColumn['y']] += y1

        corr_cat = cell_cat if corr_cat == None else numpy.append(corr_cat, cell_cat, axis=0)


    # print "\n\n\n\ntotal corrected catalog:",corr_cat.shape
    # save the source catalog
    #numpy.savetxt("focus_cat.ota%02d" % (ota), source_cat)
    #numpy.savetxt("focus_cat2.ota%02d" % (ota), corr_cat)
    logger.debug("done fixing the pixel coordinates")

    # only select bright enough sources
    bright_enough = corr_cat[:,SXFocusColumn['mag_auto']] < -10
    corr_cat = corr_cat[bright_enough]
    #numpy.savetxt("focus_cat3.ota%02d" % (ota), corr_cat)

    #dummy_test = open("dummy.test", "w")
    # Now try to match up stars in a sequence
    all_angles, all_distances = [], []
    for s1 in range(corr_cat.shape[0]):
        # Assume this is the middle star in the sequence

        # Find all stars above and below it in a cone
        # compute the angle to all other stars in the catalog
        dx = corr_cat[:,SXFocusColumn['x']] - corr_cat[s1,2]
        dy = corr_cat[:,SXFocusColumn['y']] - corr_cat[s1,3]
        d_total  = numpy.hypot(dx, dy)

        in_cone = numpy.fabs(dx/dy) < 0.1

        # Need to be at most n_stars * 10'' and at least 5'' 
        close_enough = (d_total < (n_stars * 10. / 0.11)) & (d_total > 5 / 0.11)
        
        candidates = corr_cat[in_cone & close_enough]
        # if (numpy.sum(candidates) < n_stars):
        #     # Only use full sequences
        #     continue

        # are the magnitudes comparable
        similar_brightness = numpy.fabs(candidates[:,SXFocusColumn['mag_auto']] \
                                            - corr_cat[s1,SXFocusColumn['mag_auto']]) < 1
        #print >>dummy_test, "#", candidates.shape[0], numpy.sum(similar_brightness)
        good_candidates = candidates[similar_brightness]


        # Now sort the data with increasing y values
        si = numpy.argsort(good_candidates[:,SXFocusColumn['y']])
        sorted_candidates = good_candidates[si]

        #numpy.savetxt(dummy_test, sorted_candidates)
        #print >>dummy_test, "\n\n\n"

        # Now compute the slope and distance between each point and each point 
        # above it
        for p1, p2 in itertools.combinations(range(sorted_candidates.shape[0]), 2):
            angle = numpy.arctan2(sorted_candidates[p1,2] - sorted_candidates[p2,2],
                                  sorted_candidates[p1,3] - sorted_candidates[p2,3])
            distance = numpy.sqrt( (sorted_candidates[p1,2] - sorted_candidates[p2,2])**2
                                   + (sorted_candidates[p1,3] - sorted_candidates[p2,3])**2 )
            all_angles.append(angle)
            all_distances.append(distance)

    #dummy_test.close()

    # Once we are through with the first iteration find the best-fitting angle

    all_angles = numpy.array(all_angles)
    all_distances = numpy.array(all_distances)

    # Find the best or rather most frequently occuring angle
    all_angles[all_angles < 0] += 2*math.pi

    #numpy.savetxt("dummy.angles", all_angles)
    #numpy.savetxt("dummy.distances", all_distances)
 
    filtered_angles = three_sigma_clip(all_angles)
    if (filtered_angles == None or
        filtered_angles.ndim < 1 or
        filtered_angles.shape[0] <= 0):
        return None

    angle_width = scipy.stats.scoreatpercentile(filtered_angles, [16,84])
    angle_width = scipy.stats.scoreatpercentile(filtered_angles, [5,95])

    logger.debug("Found median angle %f [%f ...%f]" % (
            numpy.degrees(numpy.median(filtered_angles)), 
            numpy.degrees(angle_width[0]), numpy.degrees(angle_width[1])
            )
    )

    #
    # Now we can do another proper search for all stars
    # This time, only search for complete series (#stars as specified)
    #
    #focus_stars = open("focus_stars", "w")
    all_candidates = []
    for s1 in range(corr_cat.shape[0]):

        # Find all stars above and below it in a cone
        # compute the angle to all other stars in the catalog
        dx = corr_cat[:,SXFocusColumn['x']] - corr_cat[s1,2]
        dy = corr_cat[:,SXFocusColumn['y']] - corr_cat[s1,3]

        angles = numpy.arctan2(dx, dy)
        angles[angles < 0] += 2*math.pi

        d_total = numpy.hypot(dx, dy)
        #print numpy.degrees(angle_width), numpy.degrees(angles)

        #print angle_width[0], angle_width[1]
        in_cone1 = (angles > angle_width[0]) & (angles < angle_width[1])
        in_cone2 = (angles+math.pi > angle_width[0]) & (angles+math.pi < angle_width[1])
        in_cone = in_cone1 | in_cone2
        # print angles[in_cone]
        #print angles[in_cone][0], angles[in_cone][0] > angle_width[0], angles[in_cone][0] < angle_width[1]
        close_enough = (d_total < ((n_stars+1) * 10. / 0.11)) & (d_total > 5 / 0.11)
        similar_brightness = numpy.fabs(corr_cat[:,SXFocusColumn['mag_auto']] \
                                            - corr_cat[s1,SXFocusColumn['mag_auto']]) < 1
        good = in_cone & close_enough & similar_brightness
        good[s1] = True
        # print s1, ":", numpy.sum(in_cone), numpy.sum(close_enough), numpy.sum(similar_brightness), numpy.sum(good)

        if (numpy.sum(good) <= 1):
            continue

        candidates = corr_cat[good]
        # print "# canddates =", candidates.shape[0]

        if (not candidates.shape[0] == n_stars): 
             # Only use full sequences
            continue
            pass

        # print "found match:",s1

        # Now we have a set with the right number of stars, matching the overall 
        # angle, and with similar brightnesses
        # sort them top to bottom
        si = numpy.argsort(candidates[:,3])
        #numpy.savetxt(focus_stars, candidates[:,3])
        #numpy.savetxt(focus_stars, si)
        sorted_candidates = candidates[si]

        sorted_candidates[:,0] = numpy.arange(sorted_candidates.shape[0])[::-1]+1.
        #numpy.savetxt(focus_stars, sorted_candidates)
        #numpy.savetxt(focus_stars, numpy.degrees(angles[good]))
        #numpy.savetxt(focus_stars, in_cone[good])
        #numpy.savetxt(focus_stars, d_total[good])
        #numpy.savetxt(focus_stars, (corr_cat[:,10] - corr_cat[s1,10])[good])
        #print >>focus_stars, "\n\n\n\n"

        all_candidates.append(sorted_candidates)
    
    #focus_stars.close()

    all_candidates = numpy.array(all_candidates)
    logger.debug(str(all_candidates.shape))

    #xxx = open("steps", "w")
    # Now compute the distances from each star to the previous

    step_vectors = []

    for i in range(1, n_stars):
        logger.debug("Candidates: %d %s\n%s" % (all_candidates.ndim, str(all_candidates.shape), str(all_candidates)))
        if (all_candidates.ndim < 1 or
            all_candidates.shape[0] <= 0):
            # We ran out of candidates
            logger.debug("We ran out of viable candidates after %s stars" % (i))
            return None

        steps = all_candidates[:,i,2:4] - all_candidates[:,i-1,2:4]
        #numpy.savetxt(xxx, steps)
        #print >>xxx, "\n\n\n\n"

        logger.debug("Computing average step size, star %d" % (i))
        logger.debug("Steps-X:\n%s" % (str(steps[:,0])))
        logger.debug("Steps-y:\n%s" % (str(steps[:,1])))
        clean_dx = three_sigma_clip(steps[:,0])
        clean_dy = three_sigma_clip(steps[:,1])

        # Check if both clean_dx and clean_dy are not empty
        logger.debug("clean-dx=%s" % (str(clean_dx)))
        logger.debug("clean-dy=%s" % (str(clean_dy)))
        if (clean_dx.ndim < 1 or clean_dx.shape[0] <= 0 or
            clean_dy.ndim < 1 or clean_dy.shape[0] <= 0):
            logger.debug("Can't find a clean dx/dy shift in iteration %d" % (i))
            return None

        dx = numpy.median(clean_dx)
        dy = numpy.median(clean_dy)

        distx = scipy.stats.scoreatpercentile(clean_dx, [16,84])
        disty = scipy.stats.scoreatpercentile(clean_dy, [16,84])
        sigma_x = 0.5 * (distx[1] - distx[0])
        sigma_y = 0.5 * (disty[1] - disty[0])

        step_vectors.append([dx, dy, sigma_x, sigma_y])

        good_steps = (steps[:,0] > (dx - 3*sigma_x)) & (steps[:,0] < (dx + 3*sigma_x)) \
            & (steps[:,1] > (dy - 3*sigma_y)) & (steps[:,1] < (dy + 3*sigma_y))

        logger.debug("before step-matching #%d: %s" % (i, str(all_candidates.shape)))
        all_candidates = all_candidates[good_steps]
        logger.debug("after step-matching: #%d: %s" % (i, str(all_candidates.shape)))

    logger.debug("%s: %s" % (filename, str(step_vectors)))

    # final_focus = open("final_focus", "w")
    # for i in range(all_candidates.shape[0]):
    #     numpy.savetxt(final_focus, all_candidates[i])
    #     print >>final_focus, "\n\n\n\n\n"
    # final_focus.close()

    logger.debug("Found %d focus stars" % (all_candidates.shape[0]))
    return all_candidates

    logger.debug("Returning final FITS table catalog")
    return None
    



def get_focus_measurement(filename, n_stars=5, mp=False):

    """

    Parameters
    ----------
    filename : string
    
        One file of the exposure. This file is mainly used to obtain the
        necessary information to create all the other filenames for this
        exposure.

    output_dir : string

        Directory to hold all the saturation catalogs. This is the directory
        that will be fed into collectcells via the -persistency command line
        flag.

    mp : bool - not used

    redo : bool

        Recreate the saturation catalog if it already exists

    Returns
    -------

    """

    logger = logging.getLogger("MeasureFocus")
    logger.info("Starting focus measurement for %s (%d *)..." % (filename, n_stars))

    if (os.path.isfile(filename)):
        # This is one of the OTA fits files
        # extract the necessary information to generate the 
        # names of all the other filenames
        try:
            hdulist = pyfits.open(filename)
        except IOError:
            logger.warning("\rProblem opening file %s...\n" % (filename))
            return
        except:
            podi_logging.log_exception()


        hdr_filename = hdulist[0].header['FILENAME']
        hdr_items = hdr_filename.split('.')
        basename = "%s.%s" % (hdr_items[0], hdr_items[1])
        hdulist.close()

        # Split the input filename to extract the directory part
        directory, dummy = os.path.split(filename)

    elif (os.path.isdir(filename)):
        # As a safety precaution, if the first parameter is the directory containing 
        # the files, extract just the ID string to be used for this script
        if (filename[-1] == "/"):
            filename = filename[:-1]

        basedir, basename = os.path.split(filename)
        directory = filename


    # Setup parallel processing
    queue        = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.JoinableQueue()
        
    number_jobs_queued = 0

    for (ota_x, ota_y) in available_ota_coords:
        ota = ota_x * 10 + ota_y

        filename = "%s/%s.%02d.fits" % (directory, basename, ota)

        if (not os.path.isfile(filename)):
            filename = "%s/%s.%02d.fits.fz" % (directory, basename, ota)
            if (not os.path.isfile(filename)):
                continue

        logger.debug("Adding file %s to task list" % (filename))

        queue.put( (filename, n_stars) )
        number_jobs_queued += 1
        # break

    # Now start all the workers
    logger.debug("Starting worker processes")
    processes = []
    for i in range(sitesetup.number_cpus):
        p = multiprocessing.Process(target=mp_measure_focus, args=(queue, return_queue, False))
        p.start()
        processes.append(p)
        time.sleep(0.01)
        
        # Tell all workers to shut down when no more data is left to work on
        queue.put( (None) )

    logger.info("Collecting catalogs for each OTA")
    all_foci = None
    for i in range(number_jobs_queued):
        focus_positions = return_queue.get()
        if (focus_positions == None):
            continue

        logger.debug("Received %d focus positions" % (focus_positions.shape[0]))

        all_foci = focus_positions if all_foci == None else numpy.append(all_foci, focus_positions, axis=0)
         #print cat_name

        return_queue.task_done()

    # Join each process to make thre they terminate(d) correctly
    logger.debug("Joining process to ensure proper termination")
    for p in processes:
        p.join()

    #print all_foci, all_foci.ndim, all_foci.shape
    if (all_foci == None or
        all_foci.ndim < 2 or
        all_foci.shape[0] <= 0):
        logger.error("Couldn't find any star patterns!")
        return

    logger.info("Found a grand total of %d focus positions" % (all_foci.shape[0]))

    stats = get_mean_focuscurve(all_foci)
    pfit, uncert, fwhm_median, fwhm_std, fwhm_cleaned, best_focus_position, best_focus = stats

    plotfilename = "%s_focus.png" % (basename)
    create_focus_plot(all_foci, stats, basename, plotfilename)

    logger.debug("all done!")
    return

def poly_fit(p, x):
    return p[2] * x**2 + p[1] * x + p[0]
def poly_err(p,x,y,err):
    fit = poly_fit(p,x)
    return (fit-y)/err

def get_mean_focuscurve(foci):

    nstars = foci.shape[1]
    print "using patterns with %d stars" % (nstars)
    positions = foci[0,:,0]

    median_focus = numpy.median(foci, axis=0)
    # print median_focus
    # print "fwhm-column:",SXFocusColumn['fwhm_world']
    fwhms = [foci[:,a,SXFocusColumn['fwhm_world']]*3600. for a in range(foci.shape[1])]
    # print fwhms

    # good_value = numpy.isfinite(fwhms)
    # for i in range(3):

    #     fwhm_sigmas = scipy.stats.scoreatpercentile(fwhms[good_value], [16,84], axis=1)
    #     fwhm_median = numpy.median(fwhms[good_value], axis=1)

    #     good_value = (fwhms > fwhm_median - 3 * (fwhm_median - fwhm_sigmas[0])) & \
    #         (fwhms < fwhm_median + 3 * (fwhm_sigmas[1] - fwhm_median))

    fwhm_cleaned = [ three_sigma_clip(fwhms[a]) for a in range(len(fwhms)) ]

    # print fwhm_cleaned

    fwhm_median = [numpy.median(fwhm_cleaned[a]) for a in range(len(fwhms)) ]
    fwhm_std = [numpy.std(fwhm_cleaned[a]) for a in range(len(fwhms)) ]

    print fwhm_median
    print fwhm_std

    # Fit a polynomial to the data, using the uncertainties 
    # in each data point as error 

    pinit = [0,0,0]
    args = (positions, fwhm_median, fwhm_std)
    fit = scipy.optimize.leastsq(poly_err, pinit, args=args, full_output=1)
    pfit = fit[0]
    uncert = numpy.sqrt(numpy.diag(fit[1]))
    print pfit, uncert

    best_focus_position = -pfit[1] / (2 * pfit[2])
    best_focus = poly_fit(pfit, best_focus_position)
    print best_focus_position,"-->",best_focus
    
    print "found", "maximum" if pfit[2] < 0 else "minimum"

    return pfit, uncert, fwhm_median, fwhm_std, fwhm_cleaned, best_focus_position, best_focus

#    all_fwhms = foci[:,:,6]
#    print all_fwhms.shape

    
def create_focus_plot(plotdata, stats, obsid, plotfile):

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    pfit, uncert, fwhm_median, fwhm_std, fwhm_cleaned, best_focus_position, best_focus = stats

    fwhm_median = numpy.array(fwhm_median)
    fwhm_std = numpy.array(fwhm_std)

    # Compute the minimum and maximum
    min_y = numpy.min(fwhm_median - 2*fwhm_std) - 0.1
    min_y = 0 if min_y < 0 else min_y
    max_y = numpy.max(fwhm_median + fwhm_std) + 0.5

    focus_step = math.fabs(plotdata[0,0,0] - plotdata[0,1,0])
    min_x = numpy.min(plotdata[0,:,0]) - 0.5*focus_step
    max_x = numpy.max(plotdata[0,:,0]) + 0.5*focus_step

    if (pfit[2] > 0):
        # This means we found a minimum
        min_x = numpy.min([best_focus_position-0.3*focus_step, min_x])
        max_x = numpy.max([best_focus_position+0.3*focus_step, max_x])


    ax.set_xlim((min_x, max_x))
    ax.set_ylim((min_y, max_y))
    ax.set_xlabel("Focus position")
    ax.set_ylabel("Image quality - FWHM [arcsec]")

    x_curve = numpy.linspace(min_x,max_x,200)
    y_curve = poly_fit(pfit, x_curve)

    # Determine focus step width  
    scatterwidth = 0.3

    # Plot all individual focus data points
    for i in range(plotdata.shape[0]):
        focus_x = plotdata[i,:,0] - 0.5*scatterwidth*focus_step + \
                  numpy.random.random(plotdata.shape[1])*focus_step*scatterwidth
        ax.scatter(focus_x, plotdata[i,:,SXFocusColumn['fwhm_world']]*3600., marker="+", c="#a0a0a0")

    # Plot the median values and the uncertainties
    ax.errorbar(x=plotdata[0,:,0], y=fwhm_median, 
                xerr=focus_step*0.1, yerr=fwhm_std,
                c='green')
    ax.scatter(plotdata[0,:,0], fwhm_median, c='green')
                

    ax.plot(x_curve, y_curve, "-", linewidth=2, c='blue')

    # Set title
    if (pfit[2] > 0):
        # This means we found a minimum
        info = "Best focus: %f at position %f" % (best_focus, best_focus_position)
    else:
        median_focus_pos = numpy.median(plotdata[0,:,0])
        slope = 2*pfit[2]*median_focus_pos + pfit[1]
        info = "No best focus found, try %s focus positions" % (
            "larger" if slope < 0 else "smaller")

    title = "Focus %(obsid)s\n%(info)s" % {"obsid": obsid,
                                         "info": info,}
    ax.set_title(title)

    # Add a little arrow pointing at the minimum
    arrow_height = 0.1 * (max_y - min_y)
    arrow_pos = best_focus_position
    ax.arrow(x=best_focus_position, y=best_focus+arrow_height, 
             dx=0., dy=-1*arrow_height,
             linewidth=1.5, color="blue",
             head_starts_at_zero=False, 
             head_width=0.02*(max_x-min_x),
             head_length=0.03*(max_y-min_y),
             length_includes_head=True,
             )
    # ax.annotate(" ", (best_focus_position, best_focus), 
    #          xytext=None, xycoords='data',
    #          textcoords='data', arrowprops=None)

    # ax.text(min_x + 0.05*(max_x-min_x), 
    #         max_y - 0.05*(max_y-min_y),
    #         "best-focus: %f at position %f" % (best_focus, best_focus_position),
    #         horizontalalignment='left',
    #         verticalalignment='bottom',
    #         fontsize=10, backgroundcolor='white')

    fig.tight_layout()

    print "Saving plot to",plotfile
    fig.savefig(plotfile)

    return


if __name__ == "__main__":

    options = podi_collectcells.read_options_from_commandline()
    podi_logging.setup_logging(options)

    n_stars = int(cmdline_arg_set_or_default('-nstars', 5))

    #output_dir = cmdline_arg_set_or_default('-persistency', '.')
    #verbose = cmdline_arg_isset("-verbose")
    for filename in get_clean_cmdline()[1:]:
        get_focus_measurement(filename, n_stars)



    podi_logging.shutdown_logging(options)
