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

"""
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy
from astLib import astWCS



if __name__ == "__main__":


    catfile = sys.argv[1]
    data = numpy.loadtxt(catfile)

    fitsfile = sys.argv[2]
    hdu = pyfits.open(fitsfile)

    for i in range(1, len(hdu)):
        wcs = astWCS.WCS(hdu[i].header, mode="pyfits")

        sel = (data[:,4] == i)

        data_sel = data[sel]

        for j in range(data_sel.shape[0]):
            x, y = data_sel[j,0], data_sel[j,1]
            _ra, _dec = wcs.pix2wcs(x-1, y-1)

            print x, y, i, data_sel[j,2], data_sel[j,3], _ra, _dec
