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
    if (sys.argv[1] == "-output"):
        output_directory = sys.argv[2]
        start = 3
        if (not os.path.isdir(output_directory)):
            os.mkdir(output_directory)
    elif (sys.argv[1] == "-keeppath"):
        keeppath = True
        start = 2
        
    # Handle all reduction flags from command line
    bias_dir, dark_dir, flatfield_dir, bpm_dir, start = read_reduction_directories(start=start, warn=False)

    # Now loop over all files and run collectcells
    for folder in get_clean_cmdline()[1:]:
        if (folder[-1] == "/"):
            folder = folder[:-1]

        directory,basename = os.path.split(folder)
        if (directory == ""):
            directory = "."
            
        if (output_directory != None):
            outputfile = "%s/%s.fits" % (output_directory, basename)
        elif (keeppath):
            outputfile = "%s/%s.fits" % (directory, basename)
        else:
            outputfile = "%s.fits" % (basename)

        stdout_write("Collecting cells from %s ==> %s\n" % (folder,outputfile))

        # Collect all cells, perform reduction and write result file
        collectcells(folder, outputfile,
                     bias_dir, dark_dir, flatfield_dir, bpm_dir)

        stdout_write("\n")
        
    stdout_write("Yippie, completely done!\n")
    
