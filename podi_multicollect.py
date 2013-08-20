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

    # Read all options from command line. This is 100% compatible 
    # with an individual call to podi_collectcells.
    options = read_options_from_commandline(None)

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

        # Collect all cells, perform reduction and write result file
        try:
            collectcells(folder, outputfile, options=options, batchmode=False)
        except:
            stdout_write("\n\n##############################\n#\n# Something terrible happened!\n")
            etype, error, stackpos = sys.exc_info()
            stdout_write("# Exception report:")
            stdout_write("#  ==> %s\n" % (error))
            print traceback.format_exc()
            stdout_write("##############################\n")

        stdout_write("\n")
        
    stdout_write("Yippie, completely done!\n")
    
