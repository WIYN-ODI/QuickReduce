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
import traceback

gain_correct_frames = False
from podi_collectcells import *
from podi_definitions import *


if __name__ == "__main__":

    # Loop over all files and collect all cells. Name output files after the basename of the input file.
    
    start = 1
    output_directory = None
    keeppath = False
    
    output_directory = cmdline_arg_set_or_default("-output", None)
    if (output_directory != None):
        if (not os.path.isdir(output_directory)):
            os.mkdir(output_directory)
    
    # Handle all reduction flags from command line
    bias_dir, dark_dir, flatfield_dir, bpm_dir, start = read_reduction_directories(start=start, warn=False)

    if (cmdline_arg_isset("-fromfile")):
        filename = get_cmdline_arg("-fromfile")
        file = open(filename, "r")
        lines = file.readlines()
        filelist = []
        for line in lines:
            if (len(line) <=1):
                continue
            if (line[0] == "#"):
                continue
            fitsfilename = line.strip().split()[0]
            filelist.append(fitsfilename)
            #print filelist
        #sys.exit(0)
    else:
        filelist = get_clean_cmdline()[1:]
        
    # For now assume that the WCS template file is located in the same directory as the executable
    root_dir, py_exe = os.path.split(os.path.abspath(sys.argv[0]))
    wcs_solution = root_dir + "/wcs_distort2.fits"
    wcs_solution = cmdline_arg_set_or_default("-wcs", wcs_solution)
    stdout_write("Using canned WCS solution from %s...\n" % (wcs_solution))

    fixwcs = cmdline_arg_isset("-fixwcs")
    
    user_wcs_offset = None
    if (cmdline_arg_isset("-wcsoffset")):
        tmp = get_cmdline_arg("-wcsoffset")
        items = tmp.split(',')
        user_wcs_offset = [float(items[0]), float(items[1])]
        stdout_write("Applying a user-defined WCS offset of %.3f, %.3f degrees\n" % (user_wcs_offset[0], user_wcs_offset[1]))
    

    # Now loop over all files and run collectcells
    for folder in filelist:
        if (folder[-1] == "/"):
            folder = folder[:-1]

        directory,basename = os.path.split(folder)
        if (directory == ""):
            directory = "."

        # Figure out into what directory we should put the output file
        if (output_directory != None):
            outputfile_dir = output_directory
        elif (cmdline_arg_isset("-keeppath")):
            outputfile_dir = directory
        else:
            outputfile_dir = "."

        # And also figure out what the filename is supposed to be
        if (cmdline_arg_isset("-formatout")):
            outputfile = get_cmdline_arg("-formatout")
        else:
            # Then assemble the actual filename
            outputfile_name = "%s.fits" % (basename)
            outputfile = "%s/%s" % (outputfile_dir, outputfile_name)


        stdout_write("Collecting cells from %s ==> %s\n" % (folder,outputfile))

        clobber_mode = not cmdline_arg_isset("-noclobber")

        # Collect all cells, perform reduction and write result file
        try:
            collectcells(folder, outputfile,
                         bias_dir, dark_dir, flatfield_dir, bpm_dir,
                         wcs_solution=wcs_solution,
                         fixwcs=fixwcs,
                         clobber_mode=clobber_mode,
                         user_wcs_offset=user_wcs_offset)
        except:
            stdout_write("\n\n##############################\n#\n# Something terrible happened!\n")
            etype, error, stackpos = sys.exc_info()
            stdout_write("# Exception report:")
            stdout_write("#  ==> %s\n" % (error))
            print traceback.format_exc()
            stdout_write("##############################\n")

        stdout_write("\n")
        
    stdout_write("Yippie, completely done!\n")
    
