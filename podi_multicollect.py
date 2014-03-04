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

"""podi_multicollect is a wrapper around podi_collectcells that allows to
batch-reduce a number of files consecutively.

Usage:
-----------

    ``podi_multicollect.py -fromfile=file.list -formatout=(format) (options)``

file.list is list of either directories named after the OBSID of the frames it
contains (e.g. ``/some/dir/o20131212T012345.1/``) or a list of files, with one file
per exposure (e.g. ``/some/dir/o20131212T012345.1/o20131212T012345.1.33.fits``). In
both cases, multicollect or collectcells automatically determines the name of
all other FITS files for this exposure.

Warning
-------

    If the input file-list contains two FITS files of the same exposure 
    (e.g. xxx.33.fits and xxx.00.fits), this results in this frame being 
    reduced twice. No harm done, just a waste of time.

**-formatout=format** 

    This allows to put information from the FITS header into the filename to 
    form informative, human-readable filenames.

    Valid format string include

        * %OBJECT
        * %OBSID
        * %FILTER
        * %EXPTIME

    These can be assembled in any order.  for example into
    ``-formatout=%OBSID___%OBJECT___%FILTER___%EXPTIME.fits``. This would be
    translated during execution into:
    ``o20131112T012345.0___M51_pointing1___odi_r___1200.fits``

    **Note** If the formatout-string contains directories, these will be created
    during execution.


**(options)**
   
    All other options are handled by collectcells, thus all collectcells 
    options can also be given to podi_multicollect.

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
from podi_commandline import *


if __name__ == "__main__":

    # Loop over all files and collect all cells. Name output files after the basename of the input file.
    
    start = 1
    output_directory = None
    keeppath = False
    
    # Read all options from command line. This is 100% compatible 
    # with an individual call to podi_collectcells.
    options = read_options_from_commandline(None)

    # Setup everything we need for logging
    podi_logging.setup_logging(options)

    output_directory = cmdline_arg_set_or_default("-output", None)
    if (output_directory != None):
        if (not os.path.isdir(output_directory)):
            os.makedirs(output_directory)
    
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
            collectcells(folder, outputfile, process_tracker=None, options=options, batchmode=False)
        except:
            stdout_write("\n\n##############################\n#\n# Something terrible happened!\n")
            etype, error, stackpos = sys.exc_info()
            stdout_write("# Exception report:")
            stdout_write("#  ==> %s\n" % (error))
            print traceback.format_exc()
            stdout_write("##############################\n")

        stdout_write("\n")
        
    # Shutdown logging to shutdown cleanly
    podi_logging.shutdown_logging(options)

    stdout_write("Yippie, completely done!\n")
    
