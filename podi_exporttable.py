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

Extract a single extension from the input FITS file.

Usage:
---------

``podi_extractextension input.fits EXT.NAME output.fits``

"""

import sys
import os
import pyfits
import numpy
from podi_definitions import *

if __name__ == "__main__":

    inputfile = get_clean_cmdline()[1]
    tablename = get_clean_cmdline()[2]
    outputfile = get_clean_cmdline()[3]

    if (not os.path.isfile(inputfile)):
        print "Unable to open file %s" % (inputfile)
        sys.exit(0)

    try:
        hdulist = pyfits.open(inputfile)
    except:
        pass
        sys.exit(0)

    try:
        table = hdulist[tablename]
    except:
        print "File %s does not contain table %s" % (inputfile, tablename)
        pass
        sys.exit(0)

    with open(outputfile, "w") as of:
    
        # print header
        stdout_write("(%s) writing header ..." % (inputfile))
        for i in range(len(table.columns)):
            col = table.columns[i]
            print >>of, "# Column %02d: %-20s %-50s [%-10s]" % (i+1, col.name, col.disp, col.unit)


        n_rows = table.data.field(0).shape[0]
        n_fields = len(table.columns)

        for row in range(n_rows):
            for field in range(n_fields):
                print >>of, table.data.field(field)[row],
            print >>of
            stdout_write("\r(%s) Writing row %d of %d ..." % (inputfile, row+1, n_rows))
        stdout_write("done!\n")


    
