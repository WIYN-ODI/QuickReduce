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
from podi_collectcells import *
from podi_imcombine import *
from podi_makeflatfield import *
import podi_matchpupilghost
import logging


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
                cell = biases[0][ext].data[x1:x2, y1:y2][bm:-bm, bm:-bm]
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

                    cell = diff_ab[x1:x2, y1:y2][bm:-bm, bm:-bm]
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


def compute_gain(flats, biases, _readnoise, binning):
    
    if (len(flats) < 2):
        # not enough flats to compute gain
        logger.warning("Only found %d flats, needs >=2!" % (len(flats)))
        return None, None

    if (len(flats) != len(biases)):
        logger.warning("Need equal number of flats and biases (%d vs %d)!" % (len(flats), len(biases)))
        return None, None

    # for i in range(len(flats)):
    #     for ext in range(len(flats[i])):
    #         if ("EXTNAME" in flats[i][ext].header):
    #             print i, ext, "-->", flats[i][ext].header['EXTNAME']
    #         else:
    #             print i, ext, "--> NO NAME"

    # for i in range(len(biases)):
    #     for ext in range(len(biases[i])):
    #         if ("EXTNAME" in biases[i][ext].header):
    #             print i, ext, "-->", biases[i][ext].header['EXTNAME']
    #         else:
    #             print i, ext, "--> NO NAME"

            
    gain = {}
    readnoise = {}

    bordermargin = bm = 25 if binning == 1 else 50
    logger.info("Computing gain and readnoise from %d flats and biases" % (len(flats)))

    for a,b in itertools.combinations(range(len(flats)), 2):
        logger.debug("Computing gain & readnoise from frame pair %d and %d" % (a,b))

        for ext in range(len(flats[a])):
            if (not is_image_extension(biases[0][ext])):
                    continue
            extname = flats[a][ext].header['EXTNAME']
            # print "starting gain for ext",extname

            if (not extname in gain):
                gain[extname] = []
            if (not extname in readnoise):
                readnoise[extname] = []

            gain_ota = numpy.ones(shape=(8,8))
            readnoise_ota = numpy.ones(shape=(8,8))

            # Make sure the current extension exists for all flats and biases
            try:
                _ = flats[a][extname]
                _ = flats[b][extname]
                _ = biases[a][extname]
                _ = biases[b][extname]
            except KeyError:
                # This extension wasn't found in one or more of the frames
                logger.warning("Couldn't find extension %s in all involved files" % (extname))
                pass
                continue

            diff_flat = flats[a][extname].data - flats[b][extname].data
            diff_bias = biases[a][extname].data - biases[b][extname].data

            for cx, cy in itertools.product(range(8), repeat=2):
                x1,x2,y1,y2 = cell2ota__get_target_region(cx, cy, binning)

                dflat = diff_flat[x1:x2, y1:y2][bm:-bm, bm:-bm]
                dflat = dflat[numpy.isfinite(dflat)]

                f1 = bottleneck.nanmean(flats[a][extname].data[x1:x2, y1:y2][bm:-bm, bm:-bm])
                f2 = bottleneck.nanmean(flats[b][extname].data[x1:x2, y1:y2][bm:-bm, bm:-bm])
                
                dbias = diff_bias[x1:x2, y1:y2][bm:-bm, bm:-bm]
                dbias = dbias[numpy.isfinite(dbias)]
                
                b1 = bottleneck.nanmean(biases[a][extname].data[x1:x2, y1:y2][bm:-bm, bm:-bm])
                b2 = bottleneck.nanmean(biases[b][extname].data[x1:x2, y1:y2][bm:-bm, bm:-bm])
                
                try:
                    _sigma_f = scipy.stats.scoreatpercentile(dflat, [16, 84])
                    sigma_f = 0.5 * (_sigma_f[1] - _sigma_f[0])
                except ValueError:
                    # This means there were likely not enough valid pixels 
                    # set the readnoise to some negative default value that
                    # is easy to ignore in all other cases
                    sigma_f = numpy.NaN


                try:
                    _sigma_b = scipy.stats.scoreatpercentile(dbias, [16, 84])
                    sigma_b = 0.5 * (_sigma_b[1] - _sigma_b[0])
                except ValueError:
                    sigma_b = numpy.NaN

                gain_ota[cx,cy] = ((f1 + f2) - (b1 + b2)) / (sigma_f**2 - sigma_b**2)
                readnoise_ota[cx, cy] = sigma_b

                logger.debug("gain: ota %s, cell %d,%d: f1/2=%f %f  B1/2=%f %f sigma_f=%f, sigma_b=%f" % (
                    extname, cx, cy,
                    f1, f2, b1, b2, sigma_f, sigma_b)
                )

            # print "appending gain_ota", gain_ota.shape, "\n",gain_ota 
            gain[extname].append(gain_ota)
            # print "appending readnoise_ota", readnoise_ota.shape, "\n",readnoise_ota 
            readnoise[extname].append(readnoise_ota)
            # print "gain-%s:\n" % (extname),gain_ota

    # Now we have gain measurements for each combination of frames
    logger.debug("Computing average gain from pair data")
    for ext in gain:
        # Convert the list of arrays to 3-d array
        all_data = numpy.array(gain[ext])
        # print all_data.shape

        # mask out all negative (=invalid) numbers so we can compute the average
        all_data[all_data < 0] = numpy.NaN

        # average the readnoise values from each pair of frames
        combined = bottleneck.nanmean(all_data, axis=0)
        # print combined.shape

        # replace the NaNs back to negative values to not cause illegal FITS 
        # header values
        combined[numpy.isnan(combined)] = -10.
            
        # and prepare the resulting Readnoise values for return
        gain[ext] = combined

    # Repeat the above for the readnoise calculation
    logger.debug("Computing average gain and readnoise from pair data")
    for ext in readnoise:
        all_data = numpy.array(readnoise[ext])
        all_data[all_data < 0] = numpy.NaN
        combined = bottleneck.nanmean(all_data, axis=0)
        combined /= math.sqrt(2)
        combined[numpy.isnan(combined)] = -10.
        readnoise[ext] = combined
 
    return gain, readnoise
  
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

            for cx, cy in range(gain[extname].shape[0], gain[extname].shape[1]):
                keyword = "GN_%02d_%d%d" % (ota, cx, cy)
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
                keyword = "RN_%02d_%d%d" % (ota, cx, cy)
                techhdu.header[keyword] = (readnoise[extname][cx,cy], "readnoise, OTA %02d, cell %d,%d" % (ota, cx, cy))
                # print "adding to header:", keyword, readnoise[extname][cx,cy]
                
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

    verbose = cmdline_arg_isset("-verbose")

    # Read the input file that has the list of files
    filelist_filename = get_clean_cmdline()[1]

    # Assign a fallback output filename if none is given 
    output_directory = get_clean_cmdline()[2]

    tmp_directory = cmdline_arg_set_or_default("-tmpdir", output_directory + "/tmp")

    logger = logging.getLogger("MakeCalibration_Init")
    
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

    compute_gain_readnoise = cmdline_arg_isset("-gainreadnoise")
    nframes_gain_readnoise = -1
    if (compute_gain_readnoise):
        nframes_gain_readnoise = int(cmdline_arg_set_or_default("-gainreadnoise", 2))
    if (nframes_gain_readnoise <= 0):
        compute_gain_readnoise = False

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

        logger.info("   %s --> %s BIN=%d" % (directory, obstype, binning))

        filter = hdulist[0].header['FILTER']
        if (obstype == "DFLAT"):
            filter_list.append(filter)
        elif (obstype == "DARK" or obstype == "BIAS"):
            filter = None
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
                                     batchmode=True,
                                     showsplash=False)
                        # With batchmode enabled, we have th save the file to disk here
                        clobberfile(bias_outfile)
                        bias_hdu.writeto(bias_outfile, clobber=True)
                    else:
                        bias_hdu = pyfits.open(bias_outfile)

                    # Save the HDU if we need it later to compute gain and read-noise
                    if (compute_gain_readnoise 
                        and len(gain_readnoise_bias[binning]) < nframes_gain_readnoise):
                        gain_readnoise_bias[binning].append(bias_hdu)

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
                    logger.debug("Running collectcells for bias-frame %s" % (cur_bias))
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

                            # Save the HDU if we need it later to compute gain and read-noise
                            if (compute_gain_readnoise 
                                and len(gain_readnoise_flat[binning]) < nframes_gain_readnoise):
                                gain_readnoise_flat[binning].append(hdu_list)

                            normalize_flatfield(None, flat_outfile, binning_x=8, binning_y=8, repeats=3, batchmode_hdu=hdu_list)

                        flats_to_stack.append(flat_outfile)
                    #print flats_to_stack

                    logger.info("Stacking %d frames into %s ..." % (len(flats_to_stack), flat_frame))
                    flat_hdus = imcombine(flats_to_stack, flat_frame, "sigmaclipmean", return_hdu=True)

                    # Relabel the file
                    flat_hdus[0].header['OBJECT'] = "master-flat %s" % (filter)

                    logger.debug("Stacking %s done!" % (flat_frame))

                    #
                    # Insert here: Compute the GAIN for each cell
                    #
                    if (compute_gain_readnoise):
                        gain, readnoise = compute_gain(gain_readnoise_flat[binning], 
                                                       gain_readnoise_bias[binning],
                                                       readnoise, binning)

                        techhdu = gain_readnoise_to_tech_hdu(flat_hdus, gain, readnoise)

                        # print gain
                        if (not gain == None):
                             # Compute average readnoise numbers and add to extensions
                             for extname in gain:
                                 all_gain = gain[extname]
                                 avg_gain = numpy.average(all_gain[all_gain > 0])
                                 std_gain = numpy.std(all_gain[all_gain > 0])

                                 flat_hdus[extname].header['GAIN'] = (avg_gain, "average OTA gain [e-/ct]")
                                 flat_hdus[extname].header['GAIN_STD'] = (std_gain, "std.dev. of OTA gain")

                                 if (not readnoise == None and extname in readnoise):
                                     all_ron = readnoise[extname] * gain[extname]
                                     valid_ron = (readnoise[extname] > 0) & (gain[extname] > 0)
                                     avg_ron = numpy.average(all_ron[valid_ron])
                                     std_ron = numpy.std(all_ron[valid_ron])
                                     
                                     avg_ron_cts = numpy.average(readnoise[extname][readnoise[extname] > 0])
                                     std_ron_cts = numpy.std(readnoise[extname][readnoise[extname] > 0])
                                     flat_hdus[extname].header['RDNOISEE'] = (avg_ron, "average readnoise in e-")
                                     flat_hdus[extname].header['RON_STDE'] = (std_ron, "std.dev. of readnoise in e-")

                                     flat_hdus[extname].header['RDNOISE'] = (avg_ron_cts, "average readnoise in counts")
                                     flat_hdus[extname].header['RON_STD'] = (std_ron_cts, "std.dev. of readnoise in counts")

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

    stdout_write("\nAll done, yippie :-)\n\n")
    logger.debug("All calibrations done successfully!")
    podi_logging.shutdown_logging(options)
