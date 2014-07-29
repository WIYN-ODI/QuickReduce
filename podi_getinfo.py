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

Create a night-log for a given list of input frame. 

This log lists, for each frame 
* filename
* type of observation: bias/dark/flat/science 
* binning
* filter name
* exposure time
* Object/target name as specified during the observation (this is not necessarily
  the real name of the target, rather what the user entered)
* pointing coordinates: Ra/Dec

The output of this file is compatible with the file-list requirements of, for 
example, podi_makecalibrations.

"""

import sys
import os
import pyfits
from podi_definitions import get_binning


if __name__ == "__main__":

    for filename in sys.argv[1:]:

        if (not os.path.isfile(filename)):
            continue

        try:
            hdulist = pyfits.open(filename)
            hdr = hdulist[0].header
        except:
            continue

        
        obstype, object, exptime, filter, dateobs = "???", "???", -999, "???", "???"
        try:
            obstype = hdr['OBSTYPE']
            object = hdr['OBJECT']
            exptime = hdr['EXPTIME']
            filter = hdr['FILTER']
            dateobs = hdr['DATE-OBS']
            binning = get_binning(hdulist[1].header)

        except:
            pass

        ra = hdr['RA']
        dec = hdr['DEC']

        print "%s %s %6s %d %15s %6.1f %60s %14s %14s" % (filename, dateobs, obstype, binning, filter, exptime, object, ra, dec)
