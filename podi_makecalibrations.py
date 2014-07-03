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
podi_makecalibrations is the main tool to create typical calibration products
(bias, dark, flat-field frames) from a list of raw frames. 

To do so, podi_makecalibrations reads a list of input files, sorts them by the 
type of frame (e.g. bias) and creates each calibration product in sequence, 
using the newly generated bias frame to correct the dark-frame, etc. Internally,
all work is done by podi_collectcells.

It transparently handles binned and unbinned data. 

Additional external calibration products to be used during the reduction
(pupilghost templates, non-linearity coefficients) are specified in the same was
as in podi_collectcells.

Temporary and intermediate calibration products (such as individual, normalized 
flat-field frames, are stored in a tmp sub-directory. By default, these files 
are deleted if no longer needed, but can be retained with the -keeptemps flag.


How to run podi_makecalibrations
--------------------------------
``podi_makecalibrations input.list calib_directory (options)``

The typical way of execution is to create a calibration directory, and run 
podi_makecalibrations from within this directory, specifying . as calib_dir.


Command line options
--------------------
* **-keeptemps**

  Do not delete the intermediate files once they are no longer needed. This 
  option is handy when testing, e.g., the quality of a pupilghost subtraction, 
  as the individual files do not need to be recomputed if podi_makecalibrations
  is called for the second time.

* **-only=(bias|dark|flat)**

  Only compute the specified calibration products. This is most commonly used to 
  inspect products as they are being created, e.g. making sure the bias-frame 
  look good before using it to reduce the dark-frames.

* **-keepprepg**
 
  Keep a copy of the flat-fields BEFORE the pupilghost has been reduced.




All other command-line options are handled by podi_collectcells, so see the full
list of command line options for this task.

Typical additional options are -nonlinearity and -pupilghost.


Methods
-------
"""

import sys
import os
import pyfits
import numpy
import scipy

gain_correct_frames = False
from podi_definitions import *
from podi_commandline import *
from podi_collectcells import *
from podi_imcombine import *
from podi_makeflatfield import *
import podi_matchpupilghost
import logging
import podi_sitesetup as sitesetup


def strip_fits_extension_from_filename(filename):
    """
    Remove the trailing .fits or .fits.fz from the given filename.
    """

    if (filename[-5:] == ".fits"):
        return filename[:-5]
    elif (filename[-8:] == ".fits.fz"):
        return filename[:-8]
    return filename

def compute_readnoise(biases, binning):

    if (len(gain_readnoise_bias[binning]) < 1):
        # Without frame there's really nothing we can do
        # But this case should never happen
        logger.warning("Can't compute readnoise, zero frames available")
        return

    readnoise = {}

    bordermargin = bm = 25 if binning == 1 else 50

    if (len(biases) == 1):
        # This is not the best way of doing it, but we can still 
        # compute a estimate of the readnoise from a single frame
        logger.info("Estimating read-out noise from a single frame")
    
        # Loop over all extensions and each cell in each extension
        # print biases[0]
        for ext in range(len(biases[0])):
            if (not is_image_extension(biases[0][ext])):
                continue
            extname = biases[0][ext].header['EXTNAME']

            readnoise_ota = numpy.ones(shape=(8,8))
            for cx, cy in itertools.product(range(8), repeat=2):
                x1,x2,y1,y2 = cell2ota__get_target_region(cx, cy, binning)
                cell = biases[0][ext].data[y1:y2, x1:x2][bm:-bm, bm:-bm]
                cell = cell[numpy.isfinite(cell)]
                
                try:
                    _ron = scipy.stats.scoreatpercentile(cell, [16, 84])
                    ron = 0.5 * (_ron[1] - _ron[0])
                except ValueError:
                    # This means there were likely not enough valid pixels 
                    # set the readnoise to some negative default value that
                    # is easy to ignore in all other cases
                    ron = -10

                # logger.debug("RON for %s, cell %d,%d: %f" % (extname, cx, cy, ron))
                readnoise_ota[cx,cy] = ron

            readnoise[extname] = readnoise_ota
            # logger.debug("All rons for ext %s: %s" % (extname, str(readnoise_ota)))
            # print "readnoise-%s:\n" % (extname),readnoise_ota

    else:
        #
        # We start here if we have more than a single frame to work on
        #
        logger.info("Computing read-out noise from %d frames" % (len(biases)))

        for a,b in itertools.combinations(range(len(biases)), 2):
            logger.debug("Working on bias pair %d and %d" % (a,b))

            for ext in range(len(biases[a])):
                if (not is_image_extension(biases[0][ext])):
                    continue
                extname = biases[a][ext].header['EXTNAME']

                if (not extname in readnoise):
                    readnoise[extname] = []

                readnoise_ota = numpy.ones(shape=(8,8))

                frame_a = biases[a][extname].data
                frame_b = biases[b][extname].data

                diff_ab = frame_a - frame_b

                for cx, cy in itertools.product(range(8), repeat=2):
                    x1,x2,y1,y2 = cell2ota__get_target_region(cx, cy, binning)

                    cell = diff_ab[y1:y2, x1:x2][bm:-bm, bm:-bm]
                    cell = cell[numpy.isfinite(cell)]
                
                    try:
                        _ron = scipy.stats.scoreatpercentile(cell, [16, 84])
                        ron = 0.5 * (_ron[1] - _ron[0])
                    except ValueError:
                        # This means there were likely not enough valid pixels 
                        # set the readnoise to some negative default value that
                        # is easy to ignore in all other cases
                        ron = -10

                    # logger.debug("RON for %s, cell %d,%d: %f" % (extname, cx, cy, ron))
                    readnoise_ota[cx,cy] = ron

                readnoise[extname].append(readnoise_ota)
                # logger.debug("All rons for ext %s: %s" % (extname, str(readnoise_ota)))
                # print "readnoise-%s:\n" % (extname),readnoise_ota


        # Now we have readnoise measurements for each combination of frames
        for ext in readnoise:
            # Convert the list of arrays to 3-d array
            all_data = numpy.array(readnoise[ext])

            # mask out all negative (=invalid) numbers so we can compute the average
            all_data[all_data < 0] = numpy.NaN

            # average the readnoise values from each pair of frames
            combined = bottleneck.nanmean(all_data, axis=0)

            # Correct for the fact we are working on difference frames
            combined /= math.sqrt(2)

            # replace the NaNs back to negative values to not cause illegal FITS 
            # header values
            combined[numpy.isnan(combined)] = -10.
            
            # and prepare the resulting Readnoise values for return
            readnoise[ext] = combined
        
    logger.debug("All rons: %s" % (str(readnoise)))
    return readnoise






def compute_overscan_levels(hdulist):

    overscan_levels = {}
    extname2id = {}

    binning = get_binning(hdulist[1].header)

    for i in range(len(hdulist)):
        if (not is_image_extension(hdulist[i])):
            continue
        extname = hdulist[i].header['EXTNAME']
        cx, cy = hdulist[i].header['WN_CELLX'], hdulist[i].header['WN_CELLY']

        overscan_region = extract_biassec_from_cell(hdulist[i].data, binning)
        overscan_level = numpy.median(overscan_region)

        overscan_levels[extname.lower()] = overscan_level
        extname2id[extname.lower()] = i

    return overscan_levels, extname2id





def compute_techdata_from_bias_flat(flatlist, biaslist, ota):
    
    logger = logging.getLogger("ComputeTechData-OTA%02d" % (ota))

    if (len(flatlist) < 2):
        # not enough flats to compute gain
        logger.warning("Only found %d flats, needs >=2!" % (len(flatlist)))
        return None, None

    if (len(flatlist) != len(biaslist)):
        logger.warning("Need equal number of flats and biases (%d vs %d)!" % (len(flatlist), len(biaslist)))
        return None, None


    gain = {}
    readnoise = {}

    binning = 1
    bordermargin = bm = 25 if binning == 1 else 50
    # logger.info("Computing gain and readnoise from %d flats and biases" % (len(flatlist)))


    #
    # Open all files and compute all overscan levels
    #
    logger.debug("computing all overscan levels")
    overscan_levels = {}
    extname2id = {}

    flat_hdus = [None] * len(flatlist)
    for i in range(len(flatlist)):
        flat_hdus[i] = pyfits.open(flatlist[i])
        #print flatlist[i]
        overscan_levels[flatlist[i]], extname2id[flatlist[i]] = compute_overscan_levels(flat_hdus[i])

    bias_hdus = [None] * len(biaslist)
    for i in range(len(biaslist)):
        #print biaslist[i]
        bias_hdus[i] = pyfits.open(biaslist[i])
        overscan_levels[biaslist[i]], extname2id[biaslist[i]] = compute_overscan_levels(bias_hdus[i])
        
    #print overscan_levels

    # Prepare some arrays for all gain, readnoise, readnoise_electron data
    all_gain_readnoise = []
    sqrt_two = math.sqrt(2)

    # Find combinations of two flats and two biases
    for a,b in itertools.combinations(range(len(flat_hdus)), 2):
        logger.debug("Computing gain & readnoise from frame pair %d and %d" % (a,b))

        gain_ron = numpy.ones(shape=(3,8,8)) * -99

        # if requested, find the non-linearity correction file for each of the frames
        if (options['nonlinearity-set']):
            def load_nonlinearity_correction(options, hdu):
                ota = hdu[1].header['WN_OTAX'] * 10 + hdu[1].header['WN_OTAY']
                mjd = hdu[0].header['MJD-OBS']
                nonlinearity_file = options['nonlinearity']
                if (options['nonlinearity'] == None or 
                    options['nonlinearity'] == "" or
                    not os.path.isfile(nonlinearity_file)):
                    nonlinearity_file = podi_nonlinearity.find_nonlinearity_coefficient_file(mjd, options)
                nonlin_data = podi_nonlinearity.load_nonlinearity_correction_table(nonlinearity_file, ota)
                # logger.info("NLCorr: %s, OTA %02d" % (nonlinearity_file, ota))
                return nonlin_data

            nonlin_flat_a = load_nonlinearity_correction(options, flat_hdus[a])
            nonlin_flat_b = load_nonlinearity_correction(options, flat_hdus[b])
            nonlin_bias_a = load_nonlinearity_correction(options, bias_hdus[a])
            nonlin_bias_b = load_nonlinearity_correction(options, bias_hdus[b])


        # loop over all cells
        for cx, cy in itertools.product( range(8), range(8) ):

            cellname = "xy%d%d" % (cx, cy)

            if (not ((cellname in extname2id[flatlist[a]]) and
                     (cellname in extname2id[flatlist[b]]) and
                     (cellname in extname2id[biaslist[a]]) and
                     (cellname in extname2id[biaslist[b]])) ):
                # This extension does not exist in one of the involved frames, 
                # skip this case (but it shouldn't happen unless the file is corrupted)
                continue

            flat_a = extract_datasec_from_cell( flat_hdus[a][extname2id[flatlist[a]][cellname]].data, binning )\
                     - overscan_levels[flatlist[a]][cellname]
            flat_b = extract_datasec_from_cell( flat_hdus[b][extname2id[flatlist[b]][cellname]].data, binning )\
                     - overscan_levels[flatlist[b]][cellname]

            bias_a = extract_datasec_from_cell( bias_hdus[a][extname2id[biaslist[a]][cellname]].data, binning )\
                     - overscan_levels[biaslist[a]][cellname]
            bias_b = extract_datasec_from_cell( bias_hdus[b][extname2id[biaslist[b]][cellname]].data, binning )\
                     - overscan_levels[biaslist[b]][cellname]

            if (options['nonlinearity-set']):
                flat_a += podi_nonlinearity.compute_cell_nonlinearity_correction(flat_a, cx, cy, nonlin_flat_a)
                flat_b += podi_nonlinearity.compute_cell_nonlinearity_correction(flat_b, cx, cy, nonlin_flat_b)
                bias_a += podi_nonlinearity.compute_cell_nonlinearity_correction(bias_a, cx, cy, nonlin_bias_a)
                bias_b += podi_nonlinearity.compute_cell_nonlinearity_correction(bias_b, cx, cy, nonlin_bias_b)

            f1 = bottleneck.nanmean(flat_a[bm:-bm, bm:-bm])
            f2 = bottleneck.nanmean(flat_b[bm:-bm, bm:-bm])

            b1 = bottleneck.nanmean(bias_a[bm:-bm, bm:-bm])
            b2 = bottleneck.nanmean(bias_b[bm:-bm, bm:-bm])

            diff_bias = bias_a[bm:-bm, bm:-bm] - bias_b[bm:-bm, bm:-bm]
            diff_flat = flat_a[bm:-bm, bm:-bm] - flat_b[bm:-bm, bm:-bm]

            sigma_bias = numpy.std(diff_bias)
            sigma_flat = numpy.std(diff_flat)


            gain_ron[0, cx, cy] = ((f1 + f2) - (b1 + b2)) / (sigma_flat**2 - sigma_bias**2)
            gain_ron[1, cx, cy] = sigma_bias / sqrt_two
            gain_ron[2, cx, cy] = gain_ron[0, cx,cy] * gain_ron[1,cx,cy]

        all_gain_readnoise.append(gain_ron)

    
    # Convert the list of arrays to 4-d array (cellx, celly, value, framepair) 
    # where value = (gain, readnoise, readnoise_electrons)

    all_data = numpy.array(all_gain_readnoise)
    # print all_data.shape

    # mask out all negative (=invalid) numbers so we can compute the average
    all_data[all_data < 0] = numpy.NaN

    # average the readnoise values from each pair of frames
    avg = bottleneck.nanmean(all_data, axis=0)
    std = bottleneck.nanstd(all_data, axis=0)

    # replace the NaNs back to negative values to not cause illegal FITS 
    # header values
    avg[numpy.isnan(avg)] = -10.
    std[numpy.isnan(std)] = -10.
            
    return avg, std



def compute_techdata_mpwrapper(input_queue, return_queue):
    
    while (True):
        lists = input_queue.get()
        if (lists == None):
            input_queue.task_done()
            break
        
        biaslist, flatlist, ota = lists

        bordermargin = bm = 25 if binning == 1 else 50
        logger.debug("Computing gain and readnoise from %d flats and biases, OTA %02d" % (len(flatlist), ota))
        
        avg, std = compute_techdata_from_bias_flat(flatlist, biaslist, ota)
        logger.debug("Returning results for OTA %02d" % (ota))

        return_queue.put((ota, avg, std))

        input_queue.task_done()

    return


def compute_techdata(calib_biaslist, calib_flatlist, output_dir, options, n_frames=2):

    for binning in calib_biaslist:
        if (not binning in calib_flatlist):
            continue

        bin_flatlist = calib_flatlist[binning]
        bin_biaslist = calib_biaslist[binning]

        for filtername in bin_flatlist:
            #print filtername

            #print bin_flatlist
            #print bin_flatlist[filtername]

            out_filename = "%s/techdata_%s_bin%d.fits" % (output_dir, filtername, binning)

            if (os.path.isfile(out_filename)):
                logger.info("Techdata for %s, bin %d already exists, skipping" % (
                        filtername, binning))
                continue

            if (len(bin_flatlist[filtername]) < n_frames or
                len(bin_biaslist) < n_frames):
                logger.warning("Insufficient data to create techdata for filter %s, binning %d (%d biases, %d flats, need >= %d)" % (
                    filtername, binning, len(bin_biaslist), len(bin_flatlist[filtername]), n_frames))
                continue

            association_table = {}

            biaslist = bin_biaslist
            flatlist = bin_flatlist[filtername]

            #print biaslist
            #print flatlist
            
            #
            # Run all otas in parallel as usual
            #
            work_queue = multiprocessing.JoinableQueue()
            result_queue = multiprocessing.Queue()

            # Figure out what files we need
            bias_construct = []
            flat_construct = []
            for filename in biaslist:
                fitspos = filename.find(".fits")
                base = filename[:fitspos-2]
                ext = filename[fitspos:]
                bias_construct.append( (base, ext) )
            for filename in flatlist:
                fitspos = filename.find(".fits")
                base = filename[:fitspos-2]
                ext = filename[fitspos:]
                flat_construct.append( (base, ext) )
                
            #
            # Construct the filenames all biases and flats for each of the 
            # extensions and send the files off for processing
            #
            #print bias_construct
            #print flat_construct
            list_of_otas_to_collect = available_ota_coords
            number_parallel_jobs = 0
            for (otax, otay) in list_of_otas_to_collect:
                ota = otax * 10 + otay
                #print otax, otay

                ota_biases = []
                ota_flats = []

                for (base, ext) in bias_construct:
                    biasfile = "%s%d%d%s" % (base, otax, otay, ext)
                    if (os.path.isfile(biasfile)):
                        ota_biases.append(biasfile)
                        collect_reduction_files_used(association_table, {"raw": biasfile})
                for (base, ext) in flat_construct:
                    flatfile = "%s%d%d%s" % (base, otax, otay, ext)
                    if (os.path.isfile(flatfile)):
                        ota_flats.append(flatfile)
                        collect_reduction_files_used(association_table, {"raw": flatfile})

                n_use = numpy.min([n_frames, 
                                   len(ota_biases), 
                                   len(ota_flats)
                            ])
                
                # Add this list of files to the work queue
                work_queue.put( (ota_biases[:n_use], ota_flats[:n_use], ota) )
                number_parallel_jobs += 1

            # Start worker processes
            worker_args = (work_queue, result_queue)
            processes = []
            for i in range(sitesetup.number_cpus):
                p = multiprocessing.Process(target=compute_techdata_mpwrapper, args=worker_args)
                p.start()
                processes.append(p)
                time.sleep(0.01)

                # Also send on termination note per process
                work_queue.put(None)

            # Prepare the Tech-HDU
            techhdu = pyfits.ImageHDU(name='TECHDATA')

            # Receive all results 
            for i in range(number_parallel_jobs):
                results = result_queue.get()
                if (results == None):
                    continue
                ota, avg, std = results
                logger.debug("Adding results for OTA %02d to tech-hdu" % (ota))

                for cx, cy in itertools.product(range(8), range(8)):

                    ids4 = "%02d%d%d" % (ota, cx, cy)
                    label = "OTA%02d,cell%d%d" % (ota, cx, cy)

                    techhdu.header["GN__%s" % ids4] = (avg[0, cx,cy], "gain, %s" % (label))
                    techhdu.header["GN_E%s" % ids4] = (std[0, cx,cy], "gain std.dev., %s" % (label))

                    techhdu.header["RN__%s" % ids4] = (avg[1, cx,cy], "readnoise [counts], %s" % (label))
                    techhdu.header["RN_E%s" % ids4] = (std[1, cx,cy], "readnoise std.dev. [counts], %s" % (label))

                    techhdu.header["RNE_%s" % ids4] = (avg[2, cx,cy], "readnoise [e-], %s" % (label))
                    techhdu.header["RNEE%s" % ids4] = (std[2, cx,cy], "readnoise std.dev. [e-], %s" % (label))

                # next cell
            # next OTA

            
            #
            # Construct a filename and write the techdata hdu to file
            #
            logger.info("Assembling output file")
            prim_header = pyfits.PrimaryHDU()

            # open any of the files to get some info about filter, etc.

            assoc_hdu = create_association_table(association_table, verbose=False)

            out_filename = "%s/techdata_%s_bin%d.fits" % (output_dir, filtername, binning)
            out_hdulist = pyfits.HDUList([prim_header, techhdu, assoc_hdu])

            clobberfile(out_filename)
            out_hdulist.writeto(out_filename, clobber=True)
            logger.info("wrote output to file: %s" % (out_filename))

        # next filter

    # next binning
    return

  
def gain_readnoise_to_tech_hdu(hdulist, gain, readnoise):

    if (gain == None and readnoise== None):
        # nothing to do
        return 

    try:
        techhdu = hdulist['TECHDATA']
    except:
        techhdu = pyfits.ImageHDU(name='TECHDATA')
        hdulist.append(techhdu)
        techhdu = hdulist['TECHDATA']

    # print gain
    # print readnoise

    #
    # Add all gain information to tech-hdu
    #
    if (not gain == None):
        for extname in gain:

            ota = int(extname[3:5]) # extname has to be of form OTAxy.SCI

            for cx, cy in itertools.product(range(gain[extname].shape[0]), range(gain[extname].shape[1])):
                keyword = "GN__%02d%d%d" % (ota, cx, cy)
                techhdu.header[keyword] = (gain[extname][cx,cy], "gain, OTA %02d, cell %d,%d" % (ota, cx, cy))
                # print "adding to header:", keyword, readnoise[extname][cx,cy]

    #
    # Add all read-noise information to tech-hdu
    #
    if (not readnoise == None):
        for extname in readnoise:
            #print extname

            ota = int(extname[3:5]) # extname has to be of form OTAxy.SCI
            #print ota

            #print readnoise[extname].shape[0], readnoise[extname].shape[1]
            for cx, cy in itertools.product(range(readnoise[extname].shape[0]), range(readnoise[extname].shape[1])):
                keyword = "RN__%02d%d%d" % (ota, cx, cy)
                techhdu.header[keyword] = (readnoise[extname][cx,cy], "readnoise, OTA %02d, cell %d,%d" % (ota, cx, cy))
                # print "adding to header:", keyword, readnoise[extname][cx,cy]
                
    #
    # If we have both readnoise AND gain, also remember the readnoise in electrons
    # 
    if (not readnoise == None and not gain == None):
        for extname in gain:
            if (not extname in readnoise):
                continue
            ota = int(extname[3:5])
            readnoise_electrons = readnoise[extname] * gain[extname]
            for cx, cy in itertools.product(range(readnoise_electrons.shape[0]), range(readnoise_electrons.shape[1])):
                keyword = "RNE_%02d%d%d" % (ota, cx, cy)
                techhdu.header[keyword] = (readnoise_electrons[cx,cy], "readnoise in e-, OTA %02d, cell %d,%d" % (ota, cx, cy))


    # print techhdu.header

    return 


if __name__ == "__main__":

    stdout_write("""\

    **********************************************************************
    * This is podi_makecalibrations                                      *
    * (c) 2012-2013: Ralf Kotulla, kotulla@uwm.edu                       *
    *                University of Wisconsin (Milwaukee & Madison)       *
    *                WIYN Observatory, Inc                               *
    *                                                                    *
    * Please acknowledge the author when using any products generated    *
    * with this tool. For comments, questions or ideas for improvement   *
    * please send an email to kotulla@uwm.edu. Thank you!                *
    **********************************************************************

""")

    # Set the options for collectcells to some reasonable start values
    options = set_default_options()
    # Then read the actual given parameters from the command line
    options = read_options_from_commandline()
    # Setup everything we need for logging
    podi_logging.setup_logging(options)
    logger = logging.getLogger("MakeCalibration_Init")
    
    verbose = cmdline_arg_isset("-verbose")

    # Read the input file that has the list of files
    filelist_filename = get_clean_cmdline()[1]

    # Assign a fallback output filename if none is given 
    output_directory = get_clean_cmdline()[2]
    if (not os.path.isdir(output_directory)):
        logger.info("Output directory (%s) does not exist, creating it" % (output_directory))
        os.mkdir(output_directory)
        if (not os.path.isdir(output_directory)):
            logger.error("Failed to create output-directory, aborting!")
            sys.exit(-1)

    tmp_directory = cmdline_arg_set_or_default("-tmpdir", output_directory + "/tmp")

    if (not os.path.isfile(filelist_filename)):
        logger.critical("Unable to open input filelist %s" % (filelist_filename))
        podi_logging.shutdown_logging(options)
        sys.exit(-1)

    if (not os.path.isdir(output_directory)):
        logger.critical("Specified output directory does not exists..." % (output_directory))
        podi_logging.shutdown_logging(options)
        sys.exit(-2)

    # We also need a "tmp" subdirectory  directory in the output directory
    if (not os.path.exists(tmp_directory)):
        logger.debug("Creating tmp directory %s" % (tmp_directory))
        os.makedirs(tmp_directory)

    compute_gain_readnoise = True #cmdline_arg_isset("-gainreadnoise")
    nframes_gain_readnoise = -1
    if (compute_gain_readnoise):
        nframes_gain_readnoise = int(cmdline_arg_set_or_default("-gainreadnoise", 2))
        logger.info("Read techdata specs: use %d frames" % (nframes_gain_readnoise))
    if (nframes_gain_readnoise <= 0):
        compute_gain_readnoise = False

    if (options['nonlinearity-set']):
        logger.info("Applying non-linearity correction")

    #
    # Read the list of files
    #

    dark_list = []
    bias_list = []

    filters = []
    flat_list = []

    stdout_write("####################\n#\n# Sighting input data\n#\n####################\n")
    _list = open(filelist_filename, "r")
    calib_file_list = []
    binning_list = []
    filter_list = []

    calib_bias_list = {}
    calib_dark_list = {}
    calib_flat_list = {}

    for full_filename in _list.readlines():
        if (len(full_filename)<=1):
            continue
        if (full_filename[0] == "#"):
            continue

        ota00 = full_filename.strip().split()[0]
        #print ota00

        directory, filename = os.path.split(ota00)
        
        hdulist = pyfits.open(ota00)
        binning = get_binning(hdulist[1].header)
        obstype = hdulist[0].header['OBSTYPE']

        if (not binning in calib_bias_list):
            calib_bias_list[binning] = []
            calib_dark_list[binning] = []
            calib_flat_list[binning] = {}

        logger.info("   %s --> %s BIN=%d" % (directory, obstype, binning))

        filter = hdulist[0].header['FILTER']
        if (obstype in ["DFLAT", "TFLAT"]):
            filter_list.append(filter)
            if (not filter in calib_flat_list[binning]):
                calib_flat_list[binning][filter] = []
            calib_flat_list[binning][filter].append(ota00)
        elif (obstype == "DARK"):
            filter = None
            calib_dark_list[binning].append(ota00)
        elif (obstype == "BIAS"):
            filter = None
            calib_bias_list[binning].append(ota00)
        else:
            logger.warning("%s is not a calibration frame" % (directory))
            hdulist.close()
            continue

        hdulist.close()
        del hdulist

        calib_entry = (ota00, obstype, filter, binning)
        calib_file_list.append(calib_entry)
        binning_list.append(binning)

    # Determine all binning values encountered
    binning_set = set(binning_list)

    # print calib_bias_list
    # print calib_dark_list
    # print calib_flat_list

    # Allocate some storage so we can save the temporary data from bias and 
    # flats that we need to compute the gain and read noise
    gain_readnoise_bias = {}
    gain_readnoise_flat = {}

    # Also create a unique set of filters. This, for now, ignores 
    # the fact that not all filters have all binning factors. These
    # cases are handled below
    filter_set = set(filter_list)

    #
    # First of all, let's combine all bias frames
    #
    logger = logging.getLogger("MakeCalibration_Bias")

    for binning in binning_set:
        gain_readnoise_bias[binning] = []

        bias_frame = "%s/bias_bin%d.fits" % (output_directory, binning)
        
        # From the full filelist, extract only the bias frames with the right bias
        bias_list = []
        for (filename, obstype, filter, bin) in calib_file_list:
            if (obstype == "BIAS" and binning == bin):
                bias_list.append(filename)

        if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "bias"): 
            stdout_write("####################\n#\n# Creating bias-frame (binning %d)\n#\n####################\n" % binning)
            bias_to_stack = []
            if (not os.path.isfile(bias_frame) or cmdline_arg_isset("-redo")):
                for cur_bias in bias_list:
                    logger.debug("Running collectcells for bias-frame %s" % (cur_bias))
                    # if (verbose): print "Collecting cells for bias",cur_bias
                    # First run collectcells
                    dummy, basename = os.path.split(cur_bias)
                    bias_outfile = "%s/bias.b%d.%s.fits" % (tmp_directory, binning, strip_fits_extension_from_filename(basename))
                    if (not os.path.isfile(bias_outfile) or cmdline_arg_isset("-redo")):
                        
                        bias_hdu = collectcells(cur_bias, bias_outfile,
                                     options=options,
                                     process_tracker=None,
                                     batchmode=False,
                                     showsplash=False)

                    # # Save the HDU if we need it later to compute gain and read-noise
                    # if (compute_gain_readnoise 
                    #     and len(gain_readnoise_bias[binning]) < nframes_gain_readnoise):
                    #     gain_readnoise_bias[binning].append(bias_hdu)

                    bias_to_stack.append(bias_outfile)
                #print bias_list

                logger.info("Stacking %d frames into %s ..." % (len(bias_to_stack), bias_frame))
                bias_hdu = imcombine(bias_to_stack, bias_frame, "sigmaclipmean", return_hdu=True)

                #
                # Compute the read-noise for each cell
                #
                if (compute_gain_readnoise):
                    readnoise = compute_readnoise(gain_readnoise_bias[binning], binning)
                    gain_readnoise_to_tech_hdu(bias_hdu, None, readnoise)
                    # print readnoise
                    if (not readnoise == None):
                        # Compute average readnoise numbers and add to extensions
                        for extname in readnoise:
                            all_ron = readnoise[extname]
                            avg_readnoise = numpy.average(all_ron[all_ron > 0])
                            bias_hdu[extname].header['RDNOISE'] = (avg_readnoise, "average readnoise [counts]")
                            std_readnoise = numpy.std(all_ron[all_ron > 0])
                            bias_hdu[extname].header['RON_STD'] = (std_readnoise, "readnoise std.dev. [counts]")
                            logger.debug("Found readnoise (%s): %.4f +/- %.4f" % (extname, avg_readnoise, std_readnoise))

                
                # Relabel the file as 'master-bias" and save to disk
                if (not bias_hdu == None):
                    bias_hdu[0].header['OBJECT'] = "master-bias"
                    bias_hdu.writeto(bias_frame, clobber=True)

                logger.debug("Stacking %s done!" % (bias_frame))
            else:
                logger.info("Bias-frame already exists, nothing to do!\n")
            if (not cmdline_arg_isset("-keeptemps")):
                for file in bias_to_stack:
                    logger.debug("Deleting tmp-file %s" % (file))
                    clobberfile(file)

    #
    # Now that we have the master bias frame, go ahead and reduce the darks
    #
    logger = logging.getLogger("MakeCalibration_Dark")
    options['normalize'] = "EXPMEAS"
    # For now set all darks to detector-glow "yes"
    for binning in binning_set:
        
        dark_frame = "%s/dark_yes_bin%d.fits" % (output_directory, binning)

        # From the full filelist, extract only the dark frames with the right binning
        dark_list = []
        for (filename, obstype, filter, bin) in calib_file_list:
            if (obstype == "DARK" and binning == bin):
                dark_list.append(filename)

        if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "dark"): 
            cmdline_opts = read_options_from_commandline()
            options['bias_dir'] = output_directory if (cmdline_opts['bias_dir'] == None) else cmdline_opts['bias_dir']
            stdout_write("####################\n#\n# Creating dark-frame (binning %d)\n#\n####################\n" % (binning))
            darks_to_stack = []
            if (not os.path.isfile(dark_frame) or cmdline_arg_isset("-redo")):
                for cur_dark in dark_list:
                    logger.debug("Running collectcells for dark-frame %s" % (cur_dark))
                    # if (verbose): print "Collecting cells for dark",cur_dark
                    # First run collectcells
                    dummy, basename = os.path.split(cur_dark)
                    dark_outfile = "%s/dark.b%d.%s.fits" % (tmp_directory, binning, strip_fits_extension_from_filename(basename))
                    if (not os.path.isfile(dark_outfile) or cmdline_arg_isset("-redo")):
                        collectcells(cur_dark, dark_outfile,
                                     process_tracker=None,
                                     options=options,
                                     batchmode=False, showsplash=False)
                    darks_to_stack.append(dark_outfile)
                #print darks_to_stack

                logger.info("Stacking %d frames into %s ..." % (len(darks_to_stack), dark_frame))
                dark_hdu = imcombine(darks_to_stack, dark_frame, "sigmaclipmean", return_hdu=True)
                
                # Relabel the file as 'master-dark" and save to disk
                if (not dark_hdu == None):
                    dark_hdu[0].header['OBJECT'] = "master-dark"
                    dark_hdu.writeto(dark_frame, clobber=True)

                logger.debug("Stacking %s done!" % (dark_frame))
            else:
                logger.info("Dark-frame already exists, nothing to do!\n")
            if (not cmdline_arg_isset("-keeptemps")):
                for file in darks_to_stack:
                    logger.debug("Deleting tmp-file %s" % (file))
                    clobberfile(file)
    if ('normalize' in options):
        del options['normalize']

    #
    # And finally, reduce the flats using the biases and darks.
    #
    logger = logging.getLogger("MakeCalibration_Flat")
    logger.debug("Flat-field filter set: %s" % (filter_set))

    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "flat"): 

        cmdline_opts = read_options_from_commandline()
        options['bias_dir'] = output_directory if (cmdline_opts['bias_dir'] == None) else cmdline_opts['bias_dir']
        options['dark_dir'] = output_directory if (cmdline_opts['dark_dir'] == None) else cmdline_opts['dark_dir']

        pupilghost_dir = options['pupilghost_dir']

        for binning in binning_set:
            for filter in filter_set:
                gain_readnoise_flat[binning] = []

                flat_frame = "%s/flat_%s_bin%d.fits" % (output_directory, filter, binning)

                # From the full filelist, extract only the dark frames with the right binning
                logger.info("Preparing files for flat-field %s (bin=%d)" % (filter, binning))

                flat_list = []
                for (filename, obstype, _filter, bin) in calib_file_list:
                    if (obstype == "DFLAT" and binning == bin and filter == _filter):
                        flat_list.append(filename)

                if (len(flat_list) <= 0):
                    continue

                # Overwrite the pupil ghost correction so we don't do it twice
                options['pupilghost_dir'] = None
                logger.debug("overwriting (for now) pupilghost dir=%s" % (pupilghost_dir))

                flats_to_stack = []
                if (not os.path.isfile(flat_frame) or cmdline_arg_isset("-redo")):
                    stdout_write("####################\n#\n# Reducing flat-field %s (binning=%d)\n#\n####################\n" % (filter, binning))
                    for cur_flat in flat_list:
                        logger.debug("Running collectcells for flat-frame %s" % (cur_flat))
                        # if (verbose): print "Collecting cells for flat",cur_flat
                        # First run collectcells
                        dummy, basename = os.path.split(cur_flat)
                        flat_outfile = "%s/nflat.b%d.%s.%s.fits" % (tmp_directory, binning, filter, strip_fits_extension_from_filename(basename))
                        if (not os.path.isfile(flat_outfile) or cmdline_arg_isset("-redo")):
                            #wcs_solution = os.path.split(os.path.abspath(sys.argv[0]))[0]+"/wcs_distort2.fits"
                            #wcs_solution = cmdline_arg_set_or_default("-wcs", wcs_solution)
                            hdu_list = collectcells(cur_flat, flat_outfile,
                                                    process_tracker=None,
                                                    options=options,
                                                    batchmode=True, showsplash=False)

                            # # Save the HDU if we need it later to compute gain and read-noise
                            # if (compute_gain_readnoise 
                            #     and len(gain_readnoise_flat[binning]) < nframes_gain_readnoise):
                            #     gain_readnoise_flat[binning].append(hdu_list)

                            normalize_flatfield(None, flat_outfile, binning_x=8, binning_y=8, repeats=3, batchmode_hdu=hdu_list)

                        flats_to_stack.append(flat_outfile)
                    #print flats_to_stack

                    logger.info("Stacking %d frames into %s ..." % (len(flats_to_stack), flat_frame))
                    flat_hdus = imcombine(flats_to_stack, flat_frame, "sigmaclipmean", return_hdu=True)

                    if (flat_hdus == None):
                        # There was a problem with the stacking, so skip to the next flat-field
                        continue

                    # Relabel the file
                    flat_hdus[0].header['OBJECT'] = "master-flat %s" % (filter)

                    logger.debug("Stacking %s done!" % (flat_frame))

                    #
                    # Now apply the pupil ghost correction 
                    # Only do this if requested via keyword -pupilghost=(dirname)
                    #
                    if (not pupilghost_dir == None): #options['pupilghost_dir'] != None):
                        # Reset the pupil ghost option to enable it here
                        options['pupilghost_dir'] = pupilghost_dir

                        logger.info("Performing pupil ghost correction ...")
                        # Get level os active filter and determine what the template filename is
                        filter_level = get_filter_level(flat_hdus[0].header)

                        pg_filename = "pupilghost_template___level_%d__bin%d.fits" % (filter_level, binning)
                        pg_template = check_filename_directory(options['pupilghost_dir'], pg_filename)
                        logger.debug("Using template file %s" % (pg_template))

                        # If we have a template for this level
                        if (os.path.isfile(pg_template)):
                            logger.debug("Using pupilghost template in  %s ... " % (pg_template))
                            pg_hdu = pyfits.open(pg_template)
                            if (filter in podi_matchpupilghost.scaling_factors and
                                podi_matchpupilghost.scaling_factors[filter] > 0):

                                scaling = podi_matchpupilghost.scaling_factors[filter]

                                # Also save a copy before the pupil ghost correction.
                                if (cmdline_arg_isset("-keepprepg")):
                                    logger.debug("Writing flat-field before pupil ghost correction ...")
                                    flat_hdus.writeto(flat_frame[:-5]+".prepg.fits", clobber=True)

                                podi_matchpupilghost.subtract_pupilghost(flat_hdus, pg_hdu, scaling)
                                flat_hdus[0].header["PUPLGOST"] = (pg_template, "p.g. template")
                                flat_hdus[0].header["PUPLGFAC"] = (scaling, "pupilghost scaling")
                                logger.debug("Pupilghost subtraction complete")
                        else:
                            logger.info("Couldn't find the pupilghost template for level %d" % (filter_level))
                            logger.debug("Missing pg-template file: %s" % (pg_template))

                    # And finally write the (maybe pupilghost-corrected) flat-field to disk
                    flat_hdus.writeto(flat_frame, clobber=True)
                else:
                    logger.info("Flatfield (%s) already exists, nothing to do!\n" % (filter))
                if (not cmdline_arg_isset("-keeptemps")):
                    for file in flats_to_stack:
                        logger.debug("Deleting tmp-file %s" % (file))
                        clobberfile(file)

    #            options['pupilghost_dir'] = pupilghost_dir

    #
    # Insert here: Compute the GAIN for each cell
    #
    logger = logging.getLogger("MakeCalibration_TechData")
    if ((not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "techdata") 
        and compute_gain_readnoise):
        logger.info("Computing gain and readnoise for each cell")
        techdatafile = "%s/techdata_bin%d.fits" % (output_directory, binning)
        compute_techdata(calib_bias_list, calib_flat_list, 
                         output_directory, options, 
                         n_frames=nframes_gain_readnoise)
                                       

    logger.debug("All calibrations done successfully!")
    podi_logging.shutdown_logging(options)

    stdout_write("\nAll done, yippie :-)\n\n")
