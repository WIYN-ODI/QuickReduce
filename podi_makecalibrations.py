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

if __name__ == "__main__":

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

    stdout_write("####################\n#\n# Sighting input data\n#\n####################\n")
    _list = open(filelist_filename, "r")
    for full_filename in _list.readlines():
        ota00 = full_filename.strip()
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
            flat_list[pos].append(directory)
        elif (obstype == "DARK"):
            dark_list.append(directory)
        elif (obstype == "BIAS"):
            bias_list.append(directory)
        else:
            stdout_write("%s is not a calibration frame" % (directory))
        hdulist.close()
        del hdulist

    
    #
    # First of all, let's combine all bias frames
    #
    stdout_write("####################\n#\n# Creating bias-frame\n#\n####################\n")
    bias_frame = "%s/bias.fits" % (output_directory)
    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "bias"): 
        bias_to_stack = []
        for cur_bias in bias_list:
            # First run collectcells
            dummy, basename = os.path.split(cur_bias)
            bias_outfile = "%s/bias.%s.fits" % (tmp_directory, basename)
            if (not os.path.isfile(bias_outfile) and not cmdline_arg_isset("-redo")):
                collectcells(cur_bias, bias_outfile,
                             bias_dir=None, dark_dir=None, flatfield_dir=None, bpm_dir=None, 
                             batchmode=False)
            bias_to_stack.append(bias_outfile)
        #print bias_list

        if (not os.path.isfile(bias_frame) and not cmdline_arg_isset("-redo")):
            stdout_write("Stacking %d frames into %s ..." % (len(bias_to_stack), bias_frame))
            imcombine(bias_to_stack, bias_frame, "median")
        if (not cmdline_arg_isset("-keeptemps")):
            for file in bias_to_stack:
                clobberfile(file)


    #
    # Now that we have the master bias frame, go ahead and reduce the darks
    #
    stdout_write("####################\n#\n# Creating dark-frame\n#\n####################\n")
    # For now set all darks to detector-glow "yes"
    dark_frame = "%s/dark_yes.fits" % (output_directory)
    if (not cmdline_arg_isset("-only") or get_cmdline_arg("-only") == "dark"): 
        darks_to_stack = []
        for cur_dark in dark_list:
            # First run collectcells
            dummy, basename = os.path.split(cur_dark)
            dark_outfile = "%s/dark.%s.fits" % (tmp_directory, basename)
            if (not os.path.isfile(dark_outfile) and not cmdline_arg_isset("-redo")):
                collectcells(cur_dark, dark_outfile,
                             bias_dir=output_directory, dark_dir=None, flatfield_dir=None, bpm_dir=None, 
                             batchmode=False)
            darks_to_stack.append(dark_outfile)
        #print darks_to_stack

        if (not os.path.isfile(dark_frame) and not cmdline_arg_isset("-redo")):
            stdout_write("Stacking %d frames into %s ..." % (len(darks_to_stack), dark_frame))
            imcombine(darks_to_stack, dark_frame, "median")
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
            stdout_write("####################\n#\n# Reducing flat-field %s\n#\n####################\n" % filter)
            flats_to_stack = []
            for cur_flat in flat_list[cur_filter_id]:
                # First run collectcells
                dummy, basename = os.path.split(cur_flat)
                flat_outfile = "%s/nflat.%s.%s.fits" % (tmp_directory, filter, basename)
                if (not os.path.isfile(flat_outfile) and not cmdline_arg_isset("-redo")):
                    hdu_list = collectcells(cur_flat, flat_outfile,
                                 bias_dir=output_directory, dark_dir=output_directory, flatfield_dir=None, bpm_dir=None, 
                                 batchmode=True)
                    normalize_flatfield(None, flat_outfile, binning_x=8, binning_y=8, repeats=3, batchmode_hdu=hdu_list)
                flats_to_stack.append(flat_outfile)
            #print flats_to_stack

            if (not os.path.isfile(flat_frame) and not cmdline_arg_isset("-redo")):
                stdout_write("Stacking %d frames into %s ..." % (len(flats_to_stack), flat_frame))
                imcombine(flats_to_stack, flat_frame, "median")
            if (not cmdline_arg_isset("-keeptemps")):
                for file in flats_to_stack:
                    clobberfile(file)

