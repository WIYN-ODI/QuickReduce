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

    # Now loop over all files and run collectcells
    for folder in get_clean_cmdline()[1:]:
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
            outputfile_name = get_cmdline_arg("-formatout")
        else:
            outputfile_name = "%s.fits" % (basename)

        # Then assemble the actual filename
        outputfile = "%s/%s" % (outputfile_dir, outputfile_name)

        stdout_write("Collecting cells from %s ==> %s\n" % (folder,outputfile))

        # Collect all cells, perform reduction and write result file
        collectcells(folder, outputfile,
                     bias_dir, dark_dir, flatfield_dir, bpm_dir)

        stdout_write("\n")
        
    stdout_write("Yippie, completely done!\n")
    
