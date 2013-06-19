#!/usr/bin/env python

#
# (c) Ralf Kotulla for WIYN/pODI
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

if __name__ == "__main__":

    verbose = cmdline_arg_isset("-verbose")

    # Read the input file that has the list of files
    filelist_filename = get_clean_cmdline()[1]

    # Assign a fallback output filename if none is given 
    output_directory = get_clean_cmdline()[2]

    tmp_directory = cmdline_arg_set_or_default("-tmpdir", output_directory + "/tmp")
    pgtempdir = cmdline_arg_set_or_default("-pupilghost", None)


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

    stdout_write("####################\n#\n# Sighting input data\n#\n####################\n")
    _list = open(filelist_filename, "r")
    for full_filename in _list.readlines():
        if (len(full_filename)<=1):
            continue
        if (full_filename[0] == "#"):
            continue

        ota00 = full_filename.strip().split()[0]
        #print ota00

        directory, filename = os.path.split(ota00)
        
        hdulist = pyfits.open(ota00)
        obstype = hdulist[0].header['OBSTYPE']
        print "   %s --> %s" % (directory, obstype)

        if (obstype == "DFLAT"):
            filter = hdulist[0].header['FILTER']
            if (not filter in filters):
                # Ok, this is a new filter
                pos = len(filters)
                filters.append(filter)
                flat_list.append([])
                #print "Found new filter", filter
            else:
                pos = filters.index(filter)
            #print "Adding frame to filter #",pos
            flat_list[pos].append(ota00)
        elif (obstype == "DARK"):
            dark_list.append(ota00)
        elif (obstype == "BIAS"):
            bias_list.append(ota00)
        else:
            stdout_write("%s is not a calibration frame" % (directory))
        hdulist.close()
        del hdulist

    
    #
    # First of all, let's combine all bias frames
    #
    bias_frame = "%s/bias.fits" % (output_directory)
    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "bias"): 
        stdout_write("####################\n#\n# Creating bias-frame\n#\n####################\n")
        bias_to_stack = []
        if (not os.path.isfile(bias_frame) or cmdline_arg_isset("-redo")):
            for cur_bias in bias_list:
                if (verbose): print "Collecting cells for bias",cur_bias
                # First run collectcells
                dummy, basename = os.path.split(cur_bias)
                bias_outfile = "%s/bias.%s.fits" % (tmp_directory, basename)
                if (not os.path.isfile(bias_outfile) or cmdline_arg_isset("-redo")):
                    collectcells(cur_bias, bias_outfile,
                                 bias_dir=None, dark_dir=None, flatfield_dir=None, bpm_dir=None, 
                                 batchmode=False)
                bias_to_stack.append(bias_outfile)
            #print bias_list

            stdout_write("Stacking %d frames into %s ..." % (len(bias_to_stack), bias_frame))
            imcombine(bias_to_stack, bias_frame, "median")
        else:
            stdout_write("Bias-frame already exists!\n")       
        if (not cmdline_arg_isset("-keeptemps")):
            for file in bias_to_stack:
                clobberfile(file)


    #
    # Now that we have the master bias frame, go ahead and reduce the darks
    #
    # For now set all darks to detector-glow "yes"
    dark_frame = "%s/dark_yes.fits" % (output_directory)
    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "dark"): 
        stdout_write("####################\n#\n# Creating dark-frame\n#\n####################\n")
        darks_to_stack = []
        if (not os.path.isfile(dark_frame) or cmdline_arg_isset("-redo")):
            for cur_dark in dark_list:
                if (verbose): print "Collecting cells for dark",cur_dark
                # First run collectcells
                dummy, basename = os.path.split(cur_dark)
                dark_outfile = "%s/dark.%s.fits" % (tmp_directory, basename)
                if (not os.path.isfile(dark_outfile) or cmdline_arg_isset("-redo")):
                    collectcells(cur_dark, dark_outfile,
                                 bias_dir=output_directory, dark_dir=None, flatfield_dir=None, bpm_dir=None, 
                                 batchmode=False)
                darks_to_stack.append(dark_outfile)
            #print darks_to_stack

            stdout_write("Stacking %d frames into %s ..." % (len(darks_to_stack), dark_frame))
            imcombine(darks_to_stack, dark_frame, "median")
        else:
            stdout_write("Darkframe already exists!\n")       
        if (not cmdline_arg_isset("-keeptemps")):
            for file in darks_to_stack:
                clobberfile(file)




    #
    # And finally, reduce the flats using the biases and darks.
    #
    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "flat"): 
        for cur_filter_id in range(len(filters)):
            filter = filters[cur_filter_id]
            flat_frame = "%s/flat_%s.fits" % (output_directory, filter)
            flats_to_stack = []
            if (not os.path.isfile(flat_frame) or cmdline_arg_isset("-redo")):
                stdout_write("####################\n#\n# Reducing flat-field %s\n#\n####################\n" % filter)
                for cur_flat in flat_list[cur_filter_id]:
                    if (verbose): print "Collecting cells for flat",cur_flat
                    # First run collectcells
                    dummy, basename = os.path.split(cur_flat)
                    flat_outfile = "%s/nflat.%s.%s.fits" % (tmp_directory, filter, basename)
                    if (not os.path.isfile(flat_outfile) or cmdline_arg_isset("-redo")):
                        wcs_solution = os.path.split(os.path.abspath(sys.argv[0]))[0]+"/wcs_distort2.fits"
                        hdu_list = collectcells(cur_flat, flat_outfile,
                                     bias_dir=output_directory, dark_dir=output_directory, flatfield_dir=None, bpm_dir=None, 
                                     wcs_solution=wcs_solution,
                                     batchmode=True)
                        normalize_flatfield(None, flat_outfile, binning_x=8, binning_y=8, repeats=3, batchmode_hdu=hdu_list)
                    flats_to_stack.append(flat_outfile)
                #print flats_to_stack

                stdout_write("Stacking %d frames into %s ..." % (len(flats_to_stack), flat_frame))
                flat_hdus = imcombine(flats_to_stack, flat_frame, "median", return_hdu=True)

                #
                # Now apply the pupil ghost correction 
                # Only do this if requested via keyword -pupilghost=(dirname)
                #
                if (pgtempdir != None):
                    # Get level os active filter and determine what the template filename is
                    filter_level = get_filter_level(flat_hdus[0].header)
                    pg_template = "%s/pupilghost_template___level_%d.fits" % (pgtempdir, filter_level)

                    # If we have a template for this level
                    if (os.path.isfile(pg_template)):
                        pg_hdu = pyfits.open(pg_template)
                        scaling = podi_matchpupilghost.scaling_factors[filter]
                        podi_matchpupilghost.subtract_pupilghost(flat_hdus, pg_hdu, scaling)
                        flat_hdus[0].header.update("PUPLGOST", pg_template, "pupilghost template")
                        flat_hdus[0].header.update("PUPLGFAC", scaling, "pupilghost scaling")

                # And finally write the (maybe pupilghost-corrected) flat-field to disk
                flat_hdus.writeto(flat_frame, clobber=True)
            else:
                stdout_write("Flatfield (%s) already exists!\n" % filter)       
            if (not cmdline_arg_isset("-keeptemps")):
                for file in flats_to_stack:
                    clobberfile(file)

