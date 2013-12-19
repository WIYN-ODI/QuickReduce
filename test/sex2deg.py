#!/usr/local/bin/python
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

import sys
import os
import ephem
import numpy

if __name__ == "__main__":

    print "Press Ctrl-D to abort"

    while (True):
        
        inp = raw_input("enter sexagesimal string: ")
        ra_dec = inp.split()
        ra, dec = ra_dec[0], ra_dec[1]
        #print ra, dec

        coord_j2000 = ephem.Equatorial(ra, dec, epoch=ephem.J2000)
        print numpy.degrees(coord_j2000.ra), numpy.degrees(coord_j2000.dec) 
