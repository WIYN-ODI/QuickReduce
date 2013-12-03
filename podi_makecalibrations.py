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
How-to use:

./podi_collectcells.py 


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

def strip_fits_extension_from_filename(filename):
    #print filename[-5:]
    #print filename[-8:]
    if (filename[-5:] == ".fits"):
        return filename[:-5]
    elif (filename[-8:] == ".fits.fz"):
        return filename[:-8]
    return filename

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

    verbose = cmdline_arg_isset("-verbose")

    # Read the input file that has the list of files
    filelist_filename = get_clean_cmdline()[1]

    # Assign a fallback output filename if none is given 
    output_directory = get_clean_cmdline()[2]

    tmp_directory = cmdline_arg_set_or_default("-tmpdir", output_directory + "/tmp")

    if (not os.path.isfile(filelist_filename)):
        stdout_write("Unable to open input filelist %s" % (filelist_filename))
        sys.exit(-1)

    if (not os.path.isdir(output_directory)):
        stdout_write("Specified output directory does not exists..." % (output_directory))
        sys.exit(-2)

    # We also need a "tmp" subdirectory  directory in the output directory
    if (not os.path.exists(tmp_directory)):
        os.makedirs(tmp_directory)

    #
    # Read the list of files
    #

    dark_list = []
    bias_list = []

    filters = []
    flat_list = []

    options = read_options_from_commandline()

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

        print "   %s --> %s BIN=%d" % (directory, obstype, binning)

        filter = hdulist[0].header['FILTER']
        if (obstype == "DFLAT"):
            filter_list.append(filter)
        elif (obstype == "DARK" or obstype == "BIAS"):
            filter = None
        else:
            stdout_write("%s is not a calibration frame" % (directory))
            hdulist.close()
            continue

        hdulist.close()
        del hdulist

        calib_entry = (ota00, obstype, filter, binning)
        calib_file_list.append(calib_entry)
        binning_list.append(binning)

    # Determine all binning values encountered
    binning_set = set(binning_list)

    # Also create a unique set of filters. This, for now, ignores 
    # the fact that not all filters have all binning factors. These
    # cases are handled below
    filter_set = set(filter_list)

    #
    # First of all, let's combine all bias frames
    #
    for binning in binning_set:
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
                    if (verbose): print "Collecting cells for bias",cur_bias
                    # First run collectcells
                    dummy, basename = os.path.split(cur_bias)
                    bias_outfile = "%s/bias.b%d.%s.fits" % (tmp_directory, binning, strip_fits_extension_from_filename(basename))
                    if (not os.path.isfile(bias_outfile) or cmdline_arg_isset("-redo")):
                        collectcells(cur_bias, bias_outfile,
                                     options=options,
                                     process_tracker=None,
                                     batchmode=False,
                                     showsplash=False)
                    bias_to_stack.append(bias_outfile)
                #print bias_list

                stdout_write("Stacking %d frames into %s ..." % (len(bias_to_stack), bias_frame))
                imcombine(bias_to_stack, bias_frame, "median")
                stdout_write("done!\n")
            else:
                stdout_write("Bias-frame already exists!\n")       
            if (not cmdline_arg_isset("-keeptemps")):
                for file in bias_to_stack:
                    clobberfile(file)

    #
    # Now that we have the master bias frame, go ahead and reduce the darks
    #
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
                    if (verbose): print "Collecting cells for dark",cur_dark
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

                stdout_write("Stacking %d frames into %s ..." % (len(darks_to_stack), dark_frame))
                imcombine(darks_to_stack, dark_frame, "median")
                stdout_write("done!\n")
            else:
                stdout_write("Darkframe already exists!\n")       
            if (not cmdline_arg_isset("-keeptemps")):
                for file in darks_to_stack:
                    clobberfile(file)

    #
    # And finally, reduce the flats using the biases and darks.
    #
    print "filter set", filter_set

    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "flat"): 

        cmdline_opts = read_options_from_commandline()
        options['bias_dir'] = output_directory if (cmdline_opts['bias_dir'] == None) else cmdline_opts['bias_dir']
        options['dark_dir'] = output_directory if (cmdline_opts['dark_dir'] == None) else cmdline_opts['dark_dir']

        pupilghost_dir = options['pupilghost_dir']

        for binning in binning_set:
            for filter in filter_set:

                flat_frame = "%s/flat_%s_bin%d.fits" % (output_directory, filter, binning)

                # From the full filelist, extract only the dark frames with the right binning
                print "Workiing on", filter, binning

                flat_list = []
                for (filename, obstype, _filter, bin) in calib_file_list:
                    if (obstype == "DFLAT" and binning == bin and filter == _filter):
                        flat_list.append(filename)

                if (len(flat_list) <= 0):
                    continue

                # Overwrite the pupil ghost correction so we don't do it twice
                options['pupilghost_dir'] = None
                print "pupilghost dir=",pupilghost_dir

                flats_to_stack = []
                if (not os.path.isfile(flat_frame) or cmdline_arg_isset("-redo")):
                    stdout_write("####################\n#\n# Reducing flat-field %s (binning=%d)\n#\n####################\n" % (filter, binning))
                    for cur_flat in flat_list:
                        if (verbose): print "Collecting cells for flat",cur_flat
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
                            normalize_flatfield(None, flat_outfile, binning_x=8, binning_y=8, repeats=3, batchmode_hdu=hdu_list)
                        flats_to_stack.append(flat_outfile)
                    #print flats_to_stack

                    stdout_write("Stacking %d frames into %s ... " % (len(flats_to_stack), flat_frame))
                    flat_hdus = imcombine(flats_to_stack, flat_frame, "median", return_hdu=True)
                    stdout_write(" done!\n")

                    #
                    # Now apply the pupil ghost correction 
                    # Only do this if requested via keyword -pupilghost=(dirname)
                    #
                    if (not pupilghost_dir == None): #options['pupilghost_dir'] != None):
                        # Also save a copy before the pupil ghost correction.
                        print "Writing flat-field before pupil ghost correction ..."
                        flat_hdus.writeto(flat_frame[:-5]+".prepg.fits", clobber=True)

                        # Reset the pupil ghost option to enable it here
                        options['pupilghost_dir'] = pupilghost_dir

                        stdout_write("Performing pupil ghost correction ...")
                        # Get level os active filter and determine what the template filename is
                        filter_level = get_filter_level(flat_hdus[0].header)

                        pg_filename = "pupilghost_template___level_%d__bin%d.fits" % (filter_level, binning)
                        pg_template = check_filename_directory(options['pupilghost_dir'], pg_filename)
                        print pg_template

                        # If we have a template for this level
                        if (os.path.isfile(pg_template)):
                            stdout_write("\n   Using file %s ... " % (pg_template))
                            pg_hdu = pyfits.open(pg_template)
                            scaling = podi_matchpupilghost.scaling_factors[filter]
                            podi_matchpupilghost.subtract_pupilghost(flat_hdus, pg_hdu, scaling)
                            flat_hdus[0].header["PUPLGOST"] = (pg_template, "p.g. template")
                            flat_hdus[0].header["PUPLGFAC"] = (scaling, "pupilghost scaling")
                            stdout_write(" done!\n")
                        else:
                            stdout_write(" --> problem:\n")
                            stdout_write("Couldn't find the pupilghost template for level %d\n" % (filter_level))
                            stdout_write("   I was looking for file %s\n" % (pg_template))

                    # And finally write the (maybe pupilghost-corrected) flat-field to disk
                    flat_hdus.writeto(flat_frame, clobber=True)
                else:
                    stdout_write("Flatfield (%s) already exists!\n" % filter)       
                if (not cmdline_arg_isset("-keeptemps")):
                    for file in flats_to_stack:
                        clobberfile(file)

    #            options['pupilghost_dir'] = pupilghost_dir

    stdout_write("\nAll done, yippie :-)\n\n")
