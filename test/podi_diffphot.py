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
import multiprocessing
import copy

import podi_sitesetup as sitesetup


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


RA = 0
DEC = 1
MJD = 2
DRA = 3
DDEC = 4

def differential_photometry(inputlist, source_coords, 
                            plot_title=None, plot_filename=None,
                            reference_star_frame=None):

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
                src_name = sc_items[5]
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
    print "\n"*10,src_names,"\n"*10

    print "source coordinate data:"
    numpy.savetxt(sys.stdout, src_params, "%14.6f")

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

    print ref_stars[:25,:2]

    # Now find sources that have good photometry in all frames that 
    # we can use as reference stars
    small_errors = (ref_stars[:, aperture_error_column] < 0.05) & \
                   (ref_stars[:, aperture_column] > 14)
    no_flags = ref_stars[:, SXcolumn['flags']] == 0
    ref_stars = ref_stars[(small_errors & no_flags)]

    # Sort by magnitude and pick the 10 brightest stars
    si = numpy.argsort(ref_stars[:, aperture_column])
    ref_sorted = ref_stars[si]
    print ref_sorted[:, aperture_column]

    n_ref = 50
    phot_refs = ref_sorted[0:n_ref,0:2]
    print phot_refs[:,0:2]

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
    numpy.savetxt("src_params", src_params)
    for i_cat in range(len(catalog)):
        print "\n"*5

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

        print src_counterpart

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
                    src_data[obj,:] = (catalog[i_cat][src_counterpart[obj][0]])
                except:
                    podi_logging.log_exception()
                    continue
                #print 'src-data=',src_data.shape

        # print src_data

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

    print target_source.shape
    print all_cats.shape

    #############################################################################
    #
    # Now compute the median and filtered magnitude of each reference star
    #
    #############################################################################

    for i in range(all_cats.shape[1]):
        numpy.savetxt("allcats_%d.cat" % (i+1), all_cats[:,i,:])

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
    diffphot_correction = numpy.median(mag_deviation, axis=1).reshape((mag_deviation.shape[0],1))
    print diffphot_correction
    print diffphot_correction.shape

    numpy.save("correction", diffphot_correction)
    numpy.savetxt("correction.dat", diffphot_correction)
    numpy.save("all_cats", all_cats)

    # As verification, apply the correction to each of the reference star catalogs
    print all_cats.shape
    all_cats[:,:, SXcolumn['mag_auto']] -= diffphot_correction
    for i in range(all_cats.shape[0]):
        numpy.savetxt("photcat%d.corr2" % (i), all_cats[i,:,SXcolumn['mag_auto']])

    #############################################################################
    #
    # Now apply the differential photometry correction to the actual source photometry
    #
    #############################################################################

    numpy.savetxt("mjd", target_mjd)
    #numpy.savetxt("target_source.raw", target_source)
    print target_source.shape
    print target_source

    # Come up with a better way of storing the data in ascii files for 
    # post-processing with external files

    print target_source.shape
    print diffphot_correction.shape

    target_corrected = target_source.copy()
    target_corrected[:, :, SXcolumn['mag_auto']] -= diffphot_correction
    # numpy.savetxt("target_source.fixed", target_corrected)

    #############################################################################
    #
    # Create some plots for the light curve before & after correction
    #
    #############################################################################

    def my_replace(inp, key, insert):
        out = ""+inp
        while (out.find(key) >= 0):
            begin_pos = out.find(key)
            end_pos = begin_pos + len(key)
            out = "%s%s%s" % (out[:begin_pos], insert, out[end_pos:])
        return out

    print src_names
    for i_source in range(src_params.shape[0]):

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

            name_tag = "%name"
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
        ax_final.errorbar(time, target_corrected[:, i_source, SXcolumn['mag_auto']], 
                    yerr=target_corrected[:, i_source, SXcolumn['mag_err_auto']], 
                    fmt='o', c='blue')
        ax_final.scatter(time, target_corrected[:, i_source, SXcolumn['mag_auto']],
                         c='blue', label='final')
        ax_final.set_ylabel("app. mag.\n[mag (12'')]")
        ax_final.legend(loc='best', borderaxespad=0.5, prop={'size':9})

        #
        # As bottom plot show correction amplitude
        #
        # Determine the plotting range
        max_corr = numpy.nanmax(numpy.fabs(diffphot_correction)) * 1.35
        ax_corr.plot(time, diffphot_correction, marker='o')
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

    return



if __name__ == "__main__":

    
    options = set_default_options()
    podi_logging.setup_logging(options)
    options = read_options_from_commandline(options)
    logger = logging.getLogger("DiffPhot: Main")


    clean_cmdline = get_clean_cmdline()[1:]
    coord_end = 0
    for i in range(len(clean_cmdline)):
        if (clean_cmdline[i] == ":"):
            coord_end = i
            break

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

    data = differential_photometry(inputlist, source_coords=source_data,
                                   plot_title=title, plot_filename=plot_filename,
                                   reference_star_frame=reference_frame
    )

    logger.info("Done, shutting down")
    podi_logging.shutdown_logging(options)


    sys.exit(0)
