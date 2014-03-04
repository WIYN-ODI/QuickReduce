#! /usr/bin/env python
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

Small tool to display the associations table in a collectcells-reduced frame.


Run as
-----------

``podi_showassociations.py (-full) input.fits``

By default, only the filenames of the associated input files are shown. With the
-full option, show_associations displays the full filenames of these files.

"""

import sys
import os
import pyfits
import numpy
from podi_definitions import *
from podi_commandline import *

if __name__ == "__main__":

    inputfile = get_clean_cmdline()[1]
    show_full_file = cmdline_arg_isset("-full")

    try:
        hdulist = pyfits.open(inputfile)
    except IOError:
        print "Can't open the input file: %s" % (inputfile)
        sys.exit(0)

    
    try:
        assoctable = hdulist['ASSOCIATIONS']
        # print assoctable
    except:
        stdout_write("Couldn't find association table, quitting\n\n")
        sys.exit(0)

    reduction_steps = assoctable.data.field('correction')
    if (show_full_file):
        filenames = assoctable.data.field('filename_full')
    else:
        filenames = assoctable.data.field('filename')

    # print reduction_steps

    stdout_write("\nAssociations for file %s ...\n" % (inputfile))
    stdout_write("-------------------------------\n")
    for step in set(reduction_steps):

        assoc_filenames = filenames[reduction_steps == step]

        for i in range(len(assoc_filenames)):
            if (i==0):
                stdout_write("% 15s : %s\n" % (step, assoc_filenames[i]))
            else:
                stdout_write("% 15s : %s\n" % (' ', assoc_filenames[i]))

    stdout_write("-------------------------------\n\n")

